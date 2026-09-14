module ProbabilisticCombinator (
    choice,
    quantumChoice,
    quantumChoiceMix,
    roundMatrix,
    test_prob_comb
) where

import Matrices
import System.Random (randomRIO)
import Data.Map (Map)
import qualified Data.Map as Map
import Data.Complex (Complex(..))

import DensityMatrices 

-- Probabilistic combinator: applies f with probability p, g with probability 1-p
choice :: (a -> a) -> (a -> a) -> Double -> a -> IO a
choice f g p x = do
    prob <- randomRIO (0, 1) :: IO Double
    if prob <= p
        then return (f x)
        else return (g x)

-- Quantum channel: randomly applies U' with probability p and U with probability 1-p
quantumChoice :: Matrix -> Matrix -> Double -> Matrix -> IO Matrix
quantumChoice u' u p v = do
    prob <- randomRIO (0, 1) :: IO Double
    let term0 = matMul (matMul u' v) (dagger u')
        term1 = matMul (matMul u v) (dagger u)
    if prob <= p
        then return term0
        else return term1

-- Deterministic quantum mixture: p*(u' v u'^dagger) + (1-p)*(u v u^dagger)
quantumChoiceMix :: Monad m => Matrix -> Matrix -> Double -> Matrix -> m Matrix
quantumChoiceMix u' u p v =
    return $ matAdd (scalarMul p (matMul (matMul u' v) (dagger u')))
                    (scalarMul (1 - p) (matMul (matMul u v) (dagger u)))

-- Example functions classic
addOne :: Int -> Int
addOne x = x + 1

subtractOne :: Int -> Int
subtractOne x = x - 1

-- Helper to round a complex number to n decimal places
roundComplex :: Int -> Complex Double -> Complex Double
roundComplex n (r :+ i) = (roundTo n r) :+ (roundTo n i)
  where
    roundTo d x = fromInteger (round $ x * (10^d)) / (10.0^^d)

-- Helper to round all elements in a matrix
roundMatrix :: Int -> [[Complex Double]] -> [[Complex Double]]
roundMatrix n = map (map (roundComplex n))

-- Example usage
test_prob_comb :: IO ()
test_prob_comb = do
    putStrLn "Classic combinator results (100 trials):"
    classicResults <- sequence [choice addOne subtractOne 0.7 5 | _ <- [1..100]]
    let counts = Map.toList $ Map.fromListWith (+) [(r, 1 :: Int) | r <- classicResults]
    putStrLn "Value | Count"
    mapM_ (\(v, c) -> putStrLn $ show v ++ "     | " ++ show c) counts

    putStrLn "\nQuantum combinator results (100 trials) x_p + h ground:"
    quantumResults <- sequence [quantumChoice x h 0.1 ground_state_density | _ <- [1..100]]
    let quantumCounts = Map.toList $ Map.fromListWith (+) [(show (roundMatrix 3 m), 1 :: Int) | m <- quantumResults]
    putStrLn "Result Matrix (as string) | Count"
    mapM_ (\(matStr, c) -> putStrLn $ matStr ++ " | " ++ show c) quantumCounts
    putStrLn "\nQuantum combinator results (100 trials)  x_p + h excited:"
    quantumResults <- sequence [quantumChoice x h 0.1 excited_state_density | _ <- [1..100]]
    let quantumCounts = Map.toList $ Map.fromListWith (+) [(show (roundMatrix 3 m), 1 :: Int) | m <- quantumResults]
    putStrLn "Result Matrix (as string) | Count"
    mapM_ (\(matStr, c) -> putStrLn $ matStr ++ " | " ++ show c) quantumCounts