{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE Rank2Types #-}
{- | Error propagation by using "Numeric.AD.Mode.Forward"

This differs from <https://github.com/leftaroundabout/uncertainly-haskell> in that the sources of
uncertainty are named, so. For example, if our calculation involves a formula:

>>> f = b * sin a

And we've measured @b@ to within 0.02 units and @a@ to within 0.01 units, then:

>>> summary $ let
>>>     b = 5 ± ("b",0.02)
>>>     a = 5 ± ("a", 0.01)
>>>  in b*sin a

This prints:

> -4.794621373315692
>     RSD in %
>       aggregate   0.49749902612566316
>       a           0.2958129155327455
>       b           -0.39999999999999997

This is the same as:

>    aggregate = sqrt(diff (\x -> x * sin 5 ) 5 ^2 * 0.02^2 + diff (\x -> 5 * sin x) 5 ^ 2 * 0.01^2)
>                        / abs(5 * sin 5)


-}
module ErrorProp.Linearized where

import Data.List
import Numeric.AD.Mode.Forward
import qualified Data.Map as M
import Text.PrettyPrint
import Text.XFormat.Show
import ErrorProp.Common

-- | linearized
data LV a = LV
      { -- | The expected value
        lv_mean :: !a,
        {- | by sensitivity it is meant:

           > (d (lv_mean m) / d "v")  == (lv_sens m ! "v")

           this is close to a proper definition of sensitivity (of Y with
           respect to x) which is

           >    (dY/dx) * (x/Y)

           which can be interpreted as the percentage change in Y resulting
           from a percentage change in x

        -}
        lv_sens :: ! (M.Map String a),
        {- | variances of uncertain input variables -}
        lv_var :: ! (M.Map String a) }
    deriving (Show, Eq)

summary :: LV Double -> Doc
summary x = hang (double (lv_mean x)) 2
            $ hang (text "RSD in %") 2
                $ vcat $ map (\[a,b] -> a $$ b) (

        [text "aggregate",
         nest gap $ double (100 * sqrt (sum
                                            [ (lv_sens x M.! k)^2 * v
                                                | (k,v) <- M.assocs (lv_var x) ])
                    / abs(lv_mean x))]
        :
        [ [text k ,
            nest gap $ double (100 * lv_sens x M.! k * sqrt v / abs(lv_mean x)) ]
                                        | (k,v) <- M.assocs (lv_var x) ])

        where (<+>) = (Text.PrettyPrint.<+>)
              gap = maximum $ 12: map length (M.keys (lv_var x))


-- * Creating 'LV'

instance ErrorProp LV where
    type ErrorPropLabel LV a = (String, a)
    a +/- (s,n) = LV a (M.singleton s 1) (M.singleton s (n*n))
    certain a = LV a M.empty M.empty


-- * Lifting normal functions to LV
llv1 :: ShowNum a => (forall s. AD s (Forward a) -> AD s (Forward a)) -> LV a -> LV a
llv1 f x = llvn (\[a] -> f $! a) [x]

llv2 :: ShowNum a => (forall s. AD s (Forward a) -> AD s (Forward a) -> AD s (Forward a)) -> LV a -> LV a -> LV a
llv2 f x y = llvn (\[a,b] -> f a b) [x,y]

llvn :: ShowNum a => (forall s. [AD s (Forward a)] -> AD s (Forward a)) -> [LV a] -> LV a
llvn f lvs = case llvnn (\x -> [f x]) lvs of [a] -> a

llvnn :: ShowNum a => (forall s. [AD s (Forward a)] -> [AD s (Forward a)]) -> [LV a] -> [LV a]
llvnn f lvs = do
    (fab,dfs) <- jacobian' f (map lv_mean lvs)
    return $ let
       m1 = foldl' (M.unionWith (+)) M.empty $ zipWith3 (\a dfa -> fmap (*dfa))
                (map lv_mean lvs)
                dfs
                (map lv_sens lvs)
       m2 = foldl' (M.unionWithKey
            (\k v1 v2 ->
                -- should be an approximate comparison?
                if v1 /= v2 then error $
                        showf ("llvnn: variance "%String%
                               " duplicated with 2 different values "%
                               String%String)
                              k (showNum v1) (showNum v2)
                else v1))
            M.empty
            (map lv_var lvs)
     in fab `seq` M.foldr seq () m1 `seq` M.foldr seq () m2 `seq` LV fab m1 m2

instance ShowNum a => Num (LV a) where
    (+) = llv2 (+)
    (-) = llv2 (-)
    (*) = llv2 (*)
    abs = llv1 abs
    signum = llv1 signum
    fromInteger = certain . fromInteger

instance Ord a => Ord (LV a) where
    compare a b = compare (lv_mean a) (lv_mean b)
instance (Real a, ShowNum a) => Real (LV a) where
    toRational = toRational . lv_mean
instance (RealFrac a, ShowNum a) => RealFrac (LV a) where
    properFraction (LV a ss nn) = let (c,d) = properFraction a
        in (c, LV d (fmap (\x -> x*a/d) ss) nn)
instance (Fractional a, ShowNum a) => Fractional (LV a) where
    fromRational = certain . fromRational
    (/) = llv2 (/)
    recip = llv1 recip

instance (Floating a, ShowNum a) => Floating (LV a) where
  pi = certain pi
  (**) = llv2 (**)
  logBase = llv2 logBase
  exp   = llv1 exp
  sqrt  = llv1 sqrt
  log   = llv1 log
  sin   = llv1 sin
  tan   = llv1 tan
  cos   = llv1 cos
  asin  = llv1 asin
  atan  = llv1 atan
  acos  = llv1 acos
  sinh  = llv1 sinh
  tanh  = llv1 tanh
  cosh  = llv1 cosh
  asinh = llv1 asinh
  atanh = llv1 atanh
  acosh = llv1 acosh
