import Data.Complex
import Data.Matrix (identity, toLists)

-- Define the Hadamard gate as a 2x2 matrix
h :: Floating a => [[a]]
h = [[ 1/sqrt 2,  1/sqrt 2]
    ,[ 1/sqrt 2, -1/sqrt 2]
    ]

-- X gate
x :: Num a => [[a]]
x = [[0, 1], 
     [1, 0]]

-- Z gate
z :: [[Complex Double]]
z = [[1 :+ 0, 0 :+ 0], 
     [0 :+ 0, 0 :+ (-1)]]

-- Y gate
y :: [[Complex Double]]
y = [[0 :+ 0, 0 :+ (-1)], 
     [0 :+ 1, 0 :+ 0]]

-- Id gate (using Data.Matrix)
id_m :: [[Int]]
id_m = toLists $ identity 2

-- Function to print the matrix
printMatrix :: Show a => [[a]] -> IO ()
printMatrix matrix = mapM_ (putStrLn . unwords . map show) matrix

-- Main function to run and print the gates
main :: IO ()
main = do
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

