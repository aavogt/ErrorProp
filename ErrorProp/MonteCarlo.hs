{-# LANGUAGE TypeFamilies #-}
{- | Monte Carlo error propagation

A well known example, where the mean of the samples (4.9) is not the same
as what you get by just squaring the original values (2^2)

>>> g <- newMTGen (Just 1)
>>> let n = 50000 in (/fromIntegral n) . V.sum . snd . V.unzip <$> sample g (2 ± ("a", 1)^2 ) n
4.9904548142954965


-}
module ErrorProp.MonteCarlo where


import ErrorProp.Common
import System.Random.Mersenne
import Control.Monad
import Data.Maybe
import Control.Applicative
import qualified Data.Map as M
import qualified Data.Vector as V
import qualified Data.Traversable as T


data MV a = MV
    { mv_mean :: a
    , mv_vars :: M.Map String (Maybe a) -- ^ variance
    , mv_sample :: M.Map String a -> a -- ^ function to get the value of the sampled variable
    }


{- | > sample g mv n

 generates a 'V.Vector' of length @n@ which contains the samples


 [todo?]

 would it make sense to integrate this approach with the Linearized?
 Second-order approaches have been taken. Implementing that could provide
 some guidance for where to sample. For example if @d2f(x,y)/dxdy@ is zero,
 while the other partial derivatives (@f,xx@ and @f,yy@) are not, then perhaps
 it makes sense for the sampling to focus on the x and y variables independently.

 latin hypercube design? other designs

 make sure / show that this plays nicely with dimensional

 benchmarks / use something better that @Map String@


 variance reduction:

   antithetic pairs:

   Control variable: monte-carlo the difference between the linear approximation (kind of like LV)
   and actual 


 if we have f(x, g(y)), it may make sense to keep g(y) constant while varying x

-}
sample ::  (MTRandom a, Floating a) => MTGen -> MV a -> Int -> IO (V.Vector (M.Map String a, a))
sample g mv n = V.forM (V.replicate n ()) $ \_ -> (\s -> (s, mv_sample mv s)) <$> T.mapM
    (\x -> scale x . toUnitVar <$> random g) (mv_vars mv)
    where
    toUnitVar x = sqrt 12 * (x-0.5)
    -- 1/12 is the variance of the uniform distribution
    scale x y = maybe 1 sqrt x * y


-- * Creating 'MV'

instance ErrorProp MV where
    type ErrorPropLabel MV a = (String, a)
    a +/- (s,n) = MV a (M.singleton s (Just (n*n))) (\m -> a + (m M.! s))
    certain x = MV x M.empty (const x)

instance (Eq a, ShowNum a) => Num (MV a) where
    (+) = liftMV2 (+)
    (-) = liftMV2 (-)
    (*) = liftMV2 (*)
    abs = liftMV1 abs
    signum = liftMV1 signum
    fromInteger = certain . fromInteger

instance Eq a => Eq (MV a) where
    MV v vs _ == MV v' vs' _ = v == v' && vs == vs'

liftMV2 :: (ShowNum a, Eq a) => (a -> a -> a) -> MV a -> MV a -> MV a
liftMV2 f (MV m1 v1 s1) (MV m2 v2 s2) =
    MV (f m1 m2)
       (M.unionWithKey mergeFn v1 v2)
       (\s -> s1 s `f` s2 s)
    where
    mergeFn k (Just a') (Just b') | a' /= b' = error $
                        unwords ["liftMV2: ", k ,
                                    " has two different variances " ,
                                    showNum a', " and " , showNum b']
    mergeFn _ a b = mplus a b

liftMV1 :: (a -> a) -> MV a -> MV a
liftMV1 f (MV m1 v1 s1) = MV (f m1) v1 (\s -> f (s1 s))

instance Ord a => Ord (MV a) where
    compare a b = compare (mv_mean a) (mv_mean b)

instance (Real a, ShowNum a) => Real (MV a) where
    toRational = toRational . mv_mean

-- | can't really be implemented correctly in the sense that
-- the remainder produced is for the mean, not the actual sample.
-- For example, if the mean is 0.99, some of the sampled values
-- will be 1.01. The value of the remainder cannot depend on the
-- value of the sample (requiring ss to escape the MV constructor),
-- so you can end up with cases where @mv_sample ((a `mod` y) > y)@
instance (RealFrac a, ShowNum a, Floating a) => RealFrac (MV a) where
    properFraction (MV a v s) = let (c,d) = properFraction a
        in (c, MV d v (\ss -> s ss *a/d))

instance (Eq a, ShowNum a, Fractional a) => Fractional (MV a) where
    fromRational = certain . fromRational
    recip = liftMV1 recip
    (/) = liftMV2 (/)

instance (Floating a, ShowNum a, RealFrac a) => Floating (MV a) where
  pi = certain pi
  (**) = liftMV2 (**)
  logBase = liftMV2 logBase
  exp   = liftMV1 exp
  sqrt  = liftMV1 sqrt
  log   = liftMV1 log
  sin   = liftMV1 sin
  tan   = liftMV1 tan
  cos   = liftMV1 cos
  asin  = liftMV1 asin
  atan  = liftMV1 atan
  acos  = liftMV1 acos
  sinh  = liftMV1 sinh
  tanh  = liftMV1 tanh
  cosh  = liftMV1 cosh
  asinh = liftMV1 asinh
  atanh = liftMV1 atanh
  acosh = liftMV1 acosh

