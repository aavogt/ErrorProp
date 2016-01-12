module ErorrProp.QuasiMC.Cukier where
-- omega relatively prime to N (Cheng Suzukawa Wolfsberg pg 3994),
-- NAG D01GYF.1
--
-- a1 = 1
--
-- ai = a^(i-1) mod p
--
-- where a is adjustible ... which is too difficult...


-- table 6 from R. I. Cukier, J. H. Schaibly, and K. E. Shuler
-- J Chem Phys vol 63, pg 1140 (1975)
cukierO = V.fromList [0, 0, 1, 5, 11, 1, 17, 23, 19, 25, 41, 31,
  23, 87, 67, 73, 85, 143, 149, 99, 119, 237,
  267, 283, 151, 385, 157, 215, 449, 163, 337,
  253, 375, 441, 673, 773, 875, 873, 587, 849,
  623, 637, 891, 943, 1171, 1225, 1335, 1725, 1663, 2019]

cukierD = V.fromList [4, 8, 6, 10, 20, 22, 32, 40, 38, 26,
  56, 62, 46, 76, 96, 60, 86, 126, 134,
  112, 92, 128, 154, 196, 34, 416, 106, 208,
  328, 198, 382, 88, 348, 186, 140, 170, 284,
  568, 302, 438, 410, 248, 448, 388, 596, 216,
  100, 488, 166, 0]

cukierN = V.fromList [0, 0, 38, 78, 142, 182, 334, 486,
  630, 806, 974, 1158, 1374, 1814, 2038, 2446, 2734, 3310,
  3838, 4174, 4702, 5542, 6174, 6854, 7110, 8182, 8934, 9590,
  11358, 11526, 13014, 14206, 15046,
  16702, 18374, 19334, 20422, 21550,
  22678, 14934, 25782, 27478, 29486,
  31486, 33950, 36550, 37854, 39814,
  41518, 43606]

-- http://en.wikipedia.org/wiki/Fourier_amplitude_sensitivity_testing
-- amj :: Vector Double
-- amj =  

cukierOmegas n | n > 50 = error "cukierOmegas: n <= 50"
cukierOmegas n =
    V.constructN n
         (\ v -> if V.length v == 0
                   then cukierO V.! (n-1)
                   else V.last v + cukierD V.! (n - V.length v))

