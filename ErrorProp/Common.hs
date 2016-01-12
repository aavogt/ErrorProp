{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE DefaultSignatures #-}
module ErrorProp.Common where

import Data.List
import Foreign.C
import Numeric.Compensated

-- | the x in @LV x@ must be in this class. An empty instance can work. For
-- example:
--
-- > instance (Compensable a, Show a) => ShowNum (Compensated a)
class (Eq a, Num a) => ShowNum a where
    -- | used to show the variance, when a source of error is given
    -- inconsistent values such as:
    --
    -- > 1 ± ("x",2) + 2 ± ("x",1)
    showNum :: a -> String
    default showNum :: Show a => a -> String
    showNum = show

instance ShowNum Double
instance ShowNum Float
instance ShowNum CDouble
instance ShowNum CFloat
instance (Compensable a, Show a) => ShowNum (Compensated a)


class ErrorProp c where
    type ErrorPropLabel c a
    (+/-) :: ShowNum a => a -> ErrorPropLabel c a -> c a
    certain :: a -> c a

-- | unicode version of +/-
a ± b = a +/- b

-- | add an additional uncertainty to a value
a +- sn = a + (0 ± sn)


