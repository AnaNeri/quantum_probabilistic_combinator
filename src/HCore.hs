-- Optimized functions for linear algebra operations using hmatrix, to be used in critical paths of the quantum combinator implementation. 
-- This module serves as a bridge between our list-of-lists matrix representation and the efficient hmatrix library, 
-- allowing us to leverage hmatrix's optimized routines for matrix multiplication and Kronecker products 
-- while maintaining our existing data structures in other parts of the codebase.

module HCore (
    HMatrix,
    toHMatrix,
    fromHMatrix,
    matMulH,
    matMulList,
    kroneckerH,
    tensorListFromLists,
    tensorListH
) where

import Numeric.LinearAlgebra (Matrix)
import qualified Numeric.LinearAlgebra as LA
import Data.Complex

-- HMatrix alias for convenience
type HMatrix = Matrix (Complex Double)

-- Convert list-of-rows to hmatrix
toHMatrix :: [[Complex Double]] -> HMatrix
toHMatrix = LA.fromLists

-- Convert back
fromHMatrix :: HMatrix -> [[Complex Double]]
fromHMatrix = LA.toLists

-- Multiply two hmatrix matrices
matMulH :: HMatrix -> HMatrix -> HMatrix
matMulH = (LA.<>)

-- Multiply lists by converting to hmatrix (helper during refactor)
matMulList :: [[Complex Double]] -> [[Complex Double]] -> [[Complex Double]]
matMulList a b = fromHMatrix $ (toHMatrix a) LA.<> (toHMatrix b)

-- Kronecker product
kroneckerH :: HMatrix -> HMatrix -> HMatrix
kroneckerH = LA.kronecker

-- Build tensor product from list-of-lists returning an HMatrix
tensorListFromLists :: [[[Complex Double]]] -> HMatrix
tensorListFromLists [] = LA.fromLists [[]]
tensorListFromLists (x:xs) = foldl (\r m -> kroneckerH r (toHMatrix m)) (toHMatrix x) xs

-- Build tensor product from a list of HMatrix
tensorListH :: [HMatrix] -> HMatrix
tensorListH [] = LA.fromLists [[]]
tensorListH (x:xs) = foldl kroneckerH x xs
