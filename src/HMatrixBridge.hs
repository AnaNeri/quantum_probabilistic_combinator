module HMatrixBridge (
    matMulHM,
    toHMatrix,
    fromHMatrix,
    tensorListHM
) where

import Numeric.LinearAlgebra (Matrix, fromLists, toLists, (<>))
import qualified Numeric.LinearAlgebra as LA
import Data.Complex

-- Kronecker and list-tensor helpers
kroneckerHM :: Matrix (Complex Double) -> Matrix (Complex Double) -> Matrix (Complex Double)
kroneckerHM = LA.kronecker

-- Build a tensor product from a list of small matrices using hmatrix kronecker,
-- returning a list-of-rows representation.
tensorListHM :: [[[Complex Double]]] -> [[Complex Double]]
tensorListHM mats =
    case mats of
        [] -> []
        (x:xs) ->
            let mx = toHMatrix x
                mres = foldl (\r m -> kroneckerHM r (toHMatrix m)) mx xs
            in fromHMatrix mres

-- Convert from list-of-rows to hmatrix Matrix
toHMatrix :: [[Complex Double]] -> Matrix (Complex Double)
toHMatrix = fromLists

-- Convert back
fromHMatrix :: Matrix (Complex Double) -> [[Complex Double]]
fromHMatrix = toLists

-- Multiply using hmatrix
matMulHM :: [[Complex Double]] -> [[Complex Double]] -> [[Complex Double]]
matMulHM a b =
    let ma = toHMatrix a
        mb = toHMatrix b
    in fromHMatrix (ma LA.<> mb)
