{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE BangPatterns #-}
-- | Fourier amplitude sensitivity test, a quasi-monte carlo method
module ErrorProp.QuasiMC where

import Math.Polynomial
import Math.Polynomial.Bernoulli
import qualified Data.Vector.Unboxed as V
import qualified Data.Vector.Generic as VG
import Data.List

import Data.Numbers.Primes
import Data.VectorSpace

import ErrorProp.QuasiMC.Cukier (cukierOmegas)

import Statistics.Sample

f1 alpha2 = constPoly 1 ^+^ ((-1)^(1 + alpha2 `mod` 2) * (2*pi)^(2*alpha2)) *^ bx
  where bx = bernoulliPoly !! (2*alpha2)

-- the explicit equations for F2,F4,F6 (equation 4.14,15,16) are correct
prop_f11 x = f11 x ~= f1 1 `evalPoly` x
  where f11 x = 1 + 2 * pi^2 * (x^2 - x + 1/6)

prop_f12 x = f12 x ~= f1 2 `evalPoly` x
  where f12 x = 1 + pi^4/45 * (1 - 30* x^2 *  (1-x)^2)

prop_f13 x = f13 x ~= f1 3 `evalPoly` x
  where f13 x = 1 + (2*pi^6/945) * (1 - 21*x^2 + 105*x^4 -126*x^5 + 42*x^6)

x ~= y = (2 * abs(x-y) / (1 + abs x + abs y)) < (1e-8 :: Double)
infix 4 ~=

-- eqn 4.10. Is the function integrated to generate the "best"
-- lattice points
f2 alpha2 xs = V.product $ V.map (f1 alpha2 `evalPoly`) xs

-- | Korbov's construction for good lattice points. Is an LCG
--
-- 4.32
korbov l n s = V.iterateN s (\p -> (p * l) `mod` n) 1

-- | try all values of L up to n/2
korbovLExhaustive alpha2 n s =
    head $
    sort [ (abs (r1-1) :: Double, l)
            | l <- [1 .. n `div` 2 ],
              let r1 = rank1 (korbov l n s) n (V.replicate s 0) (f2 alpha2) ]

-- | assume L has k prime factors. Provides candidate L values,
-- and the error in the integral used to decide that L is the best.
-- The error is absolute and relative, since the target function
-- integrates to 1.
korbovLCompositeK
  :: Int -- ^ alpha/2, roughness of the function being integrated
  -> Int -- ^ k
  -> Int -- ^ n, number of points in the quadrature
  -> Int -- ^ s, dimension
  -> [(Double, Int)] -- ^ (error, best l)
korbovLCompositeK alpha2 k n s =
    sort [ (abs (rankL l - 1) :: Double, l)
            | l <- go k ps 1 ]
  where 
    rankL l = rank1 (korbov l n s) n (V.replicate s 0) (f2 alpha2)
    go 0 _ n = [n]
    go k ps1 n = do
      q : ps1' @ (_:_) <- tails ps1
      go (k-1) ps1' (n*q)
    
    ps = takeWhile ( < round ( (fromIntegral n / 2) ** (1/fromIntegral k)))
                primes

{- | rank1 shifted lattice rule

Given

> f :: [0,1)^s -> R

calculate

> Q(Z,N,C)f = (1/N) sum_j=0^N f { j z/N + c }

is an approximation to the integral of @f@ (@If@) over that region.

equation 4.42 Sloan&Joe 1994 pg 89


If C is multivariant uniform [0,1)^s, we can do a CI for the actual
integral

P( |Q - I| < sigma v) >= 1 - (1/v^2)

with sigma^2 estimated by sum_k=1^q (Q - Qbar)^2 / ( q (q-1))

is unbiased. Quantiles are biased (Lemieux 2008). Bootstrap/jacknife resample to estimate bias.


trapezoidal / rectangle rule is good for the periodic integral here,
but possibly try "GA Evans, JR Webster 1999"


-}
rank1 :: (Integral i, RealFrac a, VG.Vector v a, VG.Vector vi i,
            VG.Vector vi a)
  => vi i -- ^ z
  -> Int -- ^ N
  -> v a -- ^ c
  -> (v a -> a) -- ^ f
  -> a
rank1 z n c f = go 0 0
  where
    zn = VG.convert $ VG.map (\zi -> fromIntegral zi / fromIntegral n) z

    go !j !accum
      | j == VG.length zn = accum
      | otherwise = go (j+1) $ f $ VG.zipWith (+) c
              $ VG.map (\zni -> fracPart $ fromIntegral j*zni) zn 

{- |

>>> map fracPart [12.3, -2.5]
[0.3, 0.5]

-}
fracPart x = x - fromInteger (floor x)


-- | fourier amplitude sensitivity test
fast :: 
  (RealFrac b, VG.Vector V.Vector b) =>
    V.Vector Int -- ^ N incommensurate frequencies (see 'korbov', 'cukierOmegas', or others)
  -> (V.Vector b -> Double) -- ^ a function of N variables which will be sampled on @[0, 1)^N@
  -> V.Vector Double -- ^ fraction of the total variance attributed to the i'th variable
fast omegas f =
    let omegaMax = VG.maximum omegas
        q = omegaMax

        xj k = VG.map (\w -> fracPart $ fromIntegral (w*k) / fromIntegral (2*q+1)) omegas
        sks = V.enumFromTo (-q) q
        fks = VG.map (f . xj) sks

        withK k = VG.map (\w -> mean
                              $ VG.zipWith (\s f -> f * k (w*s)) sks fks)
                omegas

        nq = fromIntegral (2*q+1)

        a i = withK (cos . (*(2*i*pi/nq)) . fromIntegral)
        b i = withK (sin . (*(2*i*pi/nq)) . fromIntegral)

        vjs = foldl1 (V.zipWith (+))
              $ map (V.map (^2))
              [  if odd j then b fij else a fij
                  | j <- [1 .. 2],
                    let fij = fromIntegral j ]

        normalizeFac = 2 / variance fks

    in V.map (*normalizeFac) vjs

fastKi s np k =
  let n = primes !! np
      (_, l) = korbovLExhaustive 1 n s
      omegas = korbov l n s

  in fast omegas (\v -> v V.! k)

fastKis np s = [ (k,fastKi s np k)  | k <- [0 .. s-1]]
{- properties.

> ki k = fast omegas (\x -> x V.! k)

> ki k V.! j == 1 if k==j else 0

approximately satisfied: normalization is off somewhere however

-}
