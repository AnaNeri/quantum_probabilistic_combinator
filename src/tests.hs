module Main where

import Matrices
import ProbabilisticCombinator
import Data.Complex
import qualified Data.Map as Map

-- matrices 
runMatrixTests :: IO ()
runMatrixTests = do
    putStrLn "=== Matrix Tests ==="
    test_matrices

-- probabilistic combinator 
runProbCombTests :: IO ()
runProbCombTests = do
    putStrLn "\n=== Probabilistic Combinator Tests ==="
    test_prob_comb

-- gate with depolarising noise
gateWithDepolarizing :: [[Complex Double]] -> [[Complex Double]] -> Double -> IO [[Complex Double]]
gateWithDepolarizing gate state prob = do
    -- Apply the gate: U ρ U†
    let state' = matMul gate state
    -- Apply depolarizing noise: X, Y, Z with prob/3 each
    s1 <- quantumChoice x id_m (prob/3) state'
    s2 <- quantumChoice y id_m (prob/3) s1
    s3 <- quantumChoice z id_m (prob/3) s2
    return s3

-- n-qubit identity matrix
identityN :: Int -> [[Complex Double]]
identityN n = foldl1 tensor_prod (replicate n id_m)

-- Helper: apply a single-qubit depolarizing channel to qubit k in an n-qubit system
applyDepolarizingToQubit :: Int -> Double -> [[Complex Double]] -> IO [[Complex Double]]
applyDepolarizingToQubit k p state = do
    -- Build X, Y, Z acting on qubit k (tensor with id_m elsewhere)
    let n = round (logBase 2 (fromIntegral (length state))) -- number of qubits
        opOnQubit op = foldl1 tensor_prod [if i == k then op else id_m | i <- [0..n-1]]
    s1 <- quantumChoice (opOnQubit x) id_m (p/3) state
    s2 <- quantumChoice (opOnQubit y) id_m (p/3) s1
    s3 <- quantumChoice (opOnQubit z) id_m (p/3) s2
    return s3

-- General function: apply gate, then depolarizing noise only to qubits in qubitsList
gateWithDepInQubit :: [[Complex Double]] -> [[Complex Double]] -> Double -> [Int] -> IO [[Complex Double]]
gateWithDepInQubit gate state prob qubitsList = do
    let n = round (logBase 2 (fromIntegral (length state)))
        buildNoiseOp op = foldl1 tensor_prod [if i `elem` qubitsList then op else id_m | i <- [0..n-1]]
        idN = identityN n
    let s1 = matMul (matMul gate state) (dagger gate)
    s2 <- quantumChoice (buildNoiseOp x) idN (prob/3) s1
    s3 <- quantumChoice (buildNoiseOp y) idN (prob/3) s2
    s4 <- quantumChoice (buildNoiseOp z) idN (prob/3) s3
    return s4


-- Extend a single-qubit state to n qubits (all others in |0⟩)
extendToNQubits :: Int -> [[Complex Double]] -> [[Complex Double]]
extendToNQubits n s
    | n <= 1    = s
    | otherwise = foldl (\acc _ -> tensor_prod acc ground_state_density) s [2..n]

extendSWAP :: Int -> Int -> [[Complex Double]]
extendSWAP n nQubits
    | n < 0 || n+1 >= nQubits = error "extendSWAP: invalid qubit indices"
    | otherwise = foldl1 tensor_prod ops
  where
    ops = buildOps 0
    buildOps i
      | i >= nQubits = []
      | i == n       = swap : buildOps (i+2)  -- place swap at position n (acts on n and n+1), skip n+1
      | otherwise    = id_m : buildOps (i+1)  -- identity elsewhere

-- test 1
test1 :: Double -> [[Complex Double]] -> Int -> IO ()
test1 p s n = do
    putStrLn $ "Test: X, with error X_p/3 ∘ Y_p/3 ∘ Z_p/3 with probability " ++ show p ++
               ", initial state " ++ show (roundMatrix 2 s) ++ ", " ++ show n ++ " iterations"
    results <- sequence [ do
        s1 <- gateWithDepolarizing x s p
        return s1
        | _ <- [1..n]]
    let rounded = map (roundMatrix 3) results
        counts = Map.toList $ Map.fromListWith (+) [(show m, 1 :: Int) | m <- rounded]
    putStrLn "Result Matrix | Count"
    mapM_ (\(matStr, c) -> do
        putStrLn $ show matStr ++ "Count: " ++ show c
        ) counts
    putStrLn $ "Test: HZH, with error X_p/3 ∘ Y_p/3 ∘ Z_p/3 with probability " ++ show p ++
               ", initial state " ++ show (roundMatrix 2 s) ++ ", " ++ show n ++ " iterations"
    results <- sequence [ do
        s1 <- gateWithDepolarizing h s p
        s2 <- gateWithDepolarizing z s1 p
        s3 <- gateWithDepolarizing h s2 p
        return s3
        | _ <- [1..n]]
    let rounded = map (roundMatrix 3) results
        counts = Map.toList $ Map.fromListWith (+) [(show m, 1 :: Int) | m <- rounded]
    putStrLn "Result Matrix | Count"
    mapM_ (\(matStr, c) -> do
        putStrLn $ show matStr ++ "Count: " ++ show c
        ) counts

-- test 2
test2 :: Double -> [[Complex Double]] -> Int -> IO ()
test2 p s n = do
    putStrLn $ "Test: X, with error X_p/3 ∘ Y_p/3 ∘ Z_p/3 with probability " ++ show p ++
               ", initial state " ++ show (roundMatrix 2 s) ++ ", " ++ show n ++ " iterations"
    results <- sequence [ do
        s1 <- gateWithDepolarizing x s p
        return s1
        | _ <- [1..n]]
    let rounded = map (roundMatrix 3) results
        counts = Map.toList $ Map.fromListWith (+) [(show m, 1 :: Int) | m <- rounded]
    putStrLn "Result Matrix | Count"
    mapM_ (\(matStr, c) -> do
        putStrLn $ show matStr ++ "Count: " ++ show c
        ) counts
    putStrLn $ "Test: X with bit-flip error correction, with error X_p/3 ∘ Y_p/3 ∘ Z_p/3 with probability " ++ show p ++
               ", initial state " ++ show (roundMatrix 2 s) ++ ", " ++ show n ++ " iterations"
    let n_qubits = 3
        s_n = extendToNQubits n_qubits s
    putStrLn "\n--- Initial 3-qubit state s_n ---"
    printMatrix (roundMatrix 3 s_n)
    results <- sequence [ do
        putStrLn "\n--- Step 1: CX on qubits 0,1 ---"
        let cx1 = tensor_prod cx id_m
        printMatrix cx1
        s1 <- gateWithDepInQubit cx1 s_n p [0,1]
        putStrLn "State after CX1:"
        printMatrix s1

        putStrLn "\n--- Step 2: CX on qubits 0,2 ---"
        let cx2 = matMul (matMul (extendSWAP 1 3) (tensor_prod cx id_m)) (extendSWAP 1 3)
        printMatrix cx2
        s2 <- gateWithDepInQubit cx2 s1 p [0,2]
        putStrLn "State after CX2:"
        printMatrix (roundMatrix 3 s2)

        putStrLn "\n--- Step 3: X⊗X⊗X ---"
        let xxx = tensor_prod (tensor_prod x x) x
        printMatrix xxx
        s3 <- gateWithDepInQubit xxx s2 p [0,1,2]
        putStrLn "State after XXX:"
        printMatrix (roundMatrix 3 s3)

        putStrLn "\n--- Step 4: CX on qubits 0,1 ---"
        let cx1 = tensor_prod cx id_m
        printMatrix cx1
        s4 <- gateWithDepInQubit cx1 s3 p [0,1]
        putStrLn "State after CX1 (again):"
        printMatrix (roundMatrix 3 s4)

        putStrLn "\n--- Step 5: CX on qubits 0,2 ---"
        let cx2 = matMul (matMul (extendSWAP 1 3) (tensor_prod cx id_m)) (extendSWAP 1 3)
        printMatrix cx2
        s5 <- gateWithDepInQubit cx2 s4 p [0,2]
        putStrLn "State after CX2 (again):"
        printMatrix (roundMatrix 3 s5)

        putStrLn "\n--- Step 6: Toffoli (CCX) ---"
        let invccx = matMul (matMul (matMul (extendSWAP 0 3) (extendSWAP 1 3)) ccx) (matMul (extendSWAP 0 3) (extendSWAP 1 3))
        printMatrix invccx
        s6 <- gateWithDepInQubit invccx s5 p [0,1,2]
        putStrLn "State after Toffoli (CCX):"
        printMatrix (roundMatrix 3 s6)

        return s6
        | _ <- [1..n]]
    let rounded = map (roundMatrix 3) results
        counts = Map.toList $ Map.fromListWith (+) [(show m, 1 :: Int) | m <- rounded]
    putStrLn "Result Matrix | Count"
    mapM_ (\(matStr, c) -> do
        putStrLn $ show matStr ++ "Count: " ++ show c
        ) counts

-- test 3

main :: IO()
main = do 
    --runMatrixTests
    --runProbCombTests
    --test1 0.1 ground_state_density 1000
    test2 1.0 excited_state_density 1

