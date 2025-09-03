module Matrices (
    Matrix,
    matMul,
    scalarMul,
    matAdd,
    dagger,
    tensor_prod,
    printMatrix,
    h, x, y, z, id_m, cx, ccx, cu3, swap,
    test_matrices
) where

import Data.Complex
import Data.List (transpose)
import Data.Matrix (identity, toLists)

type Matrix = [[Complex Double]]

------------------------------------------------------
-- Gates 
------------------------------------------------------
-- Hadamard gate
h :: [[Complex Double]]
h = [[ (1/sqrt 2) :+ 0,  (1/sqrt 2) :+ 0]
    ,[ (1/sqrt 2) :+ 0, (-1/sqrt 2) :+ 0]
    ]
-- X gate
x :: Num a => [[a]]
x = [[0, 1], 
     [1, 0]]
-- Z gate
z :: [[Complex Double]]
z = [[1 :+ 0, 0 :+ 0], 
     [0 :+ 0, (-1) :+ 0]]
-- Y gate
y :: [[Complex Double]]
y = [[0 :+ 0, 0 :+ (-1)], 
     [0 :+ 1, 0 :+ 0]]
-- Id gate
id_m :: Num a => [[a]]
id_m = toLists $ identity 2
-- CX gate
cx :: Num a => [[a]]
cx = [[1,0,0,0],
      [0,1,0,0],
      [0,0,0,1],
      [0,0,1,0]]
-- CCX gate 
ccx :: Num a => [[a]]
ccx = [[1,0,0,0,0,0,0,0],
       [0,1,0,0,0,0,0,0],
       [0,0,1,0,0,0,0,0],
       [0,0,0,1,0,0,0,0],
       [0,0,0,0,1,0,0,0],
       [0,0,0,0,0,1,0,0],
       [0,0,0,0,0,0,0,1],
       [0,0,0,0,0,0,1,0]]
-- CU3 gate 
cu3 :: Double -> Double -> Double -> [[Complex Double]]
cu3 theta phi lambda =
    [ [1, 0, 0, 0]
    , [0, 1, 0, 0]
    , [0, 0, (cos (theta/2) :+ 0), - (exp (0 :+ lambda)) * (sin (theta/2) :+ 0)]
    , [0, 0, (exp (0 :+ phi)) * (sin (theta/2) :+ 0), (exp (0 :+ (phi + lambda))) * (cos (theta/2) :+ 0)]
    ]
-- swap gate 
swap :: Num a => [[a]]
swap = [[1,0,0,0],
        [0,0,1,0],
        [0,1,0,0],
        [0,0,0,1]]

------------------------------------------------------
-- Operations 
------------------------------------------------------
-- matrix multiplications
matMul :: Matrix -> Matrix -> Matrix
matMul a b = [[sum [x*y | (x,y) <- zip row col] | col <- transpose b] | row <- a]

-- multiplication between a scalar and a matrix
scalarMul :: Double -> Matrix -> Matrix
scalarMul s = map (map (* (s :+ 0)))

-- matrix addition
matAdd :: Matrix -> Matrix -> Matrix
matAdd a b = [[x + y | (x, y) <- zip rowA rowB] | (rowA, rowB) <- zip a b]

-- dagger
dagger :: Matrix -> Matrix
dagger m = map (map conjugate) (transpose m)

-- tensor product
tensor_prod :: Matrix -> Matrix -> Matrix
tensor_prod a b = concatMap (\rowA -> map (\rowB -> concatMap (\aElem -> map (aElem *) rowB) rowA) b) a

-- Function to print the matrix
printMatrix :: Show a => [[a]] -> IO ()
printMatrix matrix = mapM_ (putStrLn . unwords . map show) matrix

------------------------------------------------------
------------------------------------------------------
-- Main function to run and print the gates
test_matrices :: IO ()
test_matrices = do
    putStrLn "Hadamard Gate:"
    printMatrix h
    putStrLn "\nX Gate:"
    printMatrix x
    putStrLn "\nZ Gate:"
    printMatrix z
    putStrLn "\nY Gate:"
    printMatrix y
    putStrLn "\nId Gate:"
    printMatrix id_m
    putStrLn "\nMultiplication of matrices (hzh):"
    let t1 = matMul h (matMul z h)
    printMatrix t1
    putStrLn "\nMultiplication of matrix and scalar (2+id):"
    let t2 = scalarMul 2 id_m
    printMatrix t2
    putStrLn "\nAddition of matrices (id+id):"
    let t3 = matAdd id_m id_m
    printMatrix t3
    putStrLn "\nMatrix conjugate transpose:"
    let m4_start = [[1:+0, (-2) :+ (-1), 5:+ 0],  [1:+1, 0:+1, 4:+(-2)]]
    let m4_end = [[1 :+ 0, 1 :+ (-1)],  [(-2) :+ 1, 0 :+ (-1)],  [5 :+ 0, 4 :+ 2]]
    let t4 = dagger m4_start
    putStrLn $ "should be true:" ++ show (m4_end == t4)
    putStrLn "\nTensor Product (x by h)"
    let t5 = tensor_prod x h
    printMatrix t5

