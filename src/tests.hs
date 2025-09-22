module Main where

import Matrices
import ProbabilisticCombinator
import Data.Complex
import qualified Data.Map as Map
import System.IO (withFile, IOMode(WriteMode), hPutStrLn)
import System.Directory (createDirectoryIfMissing)
import Text.Printf (printf)
import System.Environment (getArgs)

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
    let state' = matMul (matMul gate state) (dagger gate)
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

-- test 1 ---------------------------------------------------------------------------
-- different implementation of the same system may have different noise
-- X = HZH theoretically but the noise depends on circuit's depth 
-------------------------------------------------------------------------------------
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

-- test 2 -----------------------------------------------------------------
-- X vs X with bit-flip error correction scheme
-- bit-flip is not enough for error correction of faulty gates
-- may be helpful if it is in noise transmission, 
-- where correction gates have less noise than the working system
--------------------------------------------------------------------------- 

-- Projector for |0⟩ on a single qubit
proj0 :: [[Complex Double]]
proj0 = [[1 :+ 0, 0 :+ 0],
         [0 :+ 0, 0 :+ 0]]

-- Projector for |1⟩ on a single qubit
proj1 :: [[Complex Double]]
proj1 = [[0 :+ 0, 0 :+ 0],
         [0 :+ 0, 1 :+ 0]]

-- Build the projector for measuring qubit q in n-qubit system, for outcome 0 or 1
buildProjector :: Int -> Int -> Int -> [[Complex Double]]
buildProjector nQubits q outcome =
    foldl1 tensor_prod [if i == q then (if outcome == 0 then proj0 else proj1) else id_m | i <- [0..nQubits-1]]

-- Trace of a square matrix
traceM :: [[Complex Double]] -> Complex Double
traceM m = sum [m !! i !! i | i <- [0..length m - 1]]

-- Measure density matrix s in qubit q, return counts for '0' and '1'
measure :: [[Complex Double]] -> Int -> [(Char, Double)]
measure s q =
    let nQubits = round (logBase 2 (fromIntegral (length s)))
        p0 = realPart $ traceM $ matMul (buildProjector nQubits q 0) s
        p1 = realPart $ traceM $ matMul (buildProjector nQubits q 1) s
    in [('0', p0), ('1', p1)]

test2 :: Double -> [[Complex Double]] -> Int -> IO ()
test2 p s n = do
    putStrLn $ "Test: X, with error X_p/3 ∘ Y_p/3 ∘ Z_p/3 with probability " ++ show p ++
               ", initial state " ++ show (roundMatrix 2 s) ++ ", " ++ show n ++ " iterations"
    results <- sequence [ do
        s1 <- gateWithDepolarizing x s p
        return s1
        | _ <- [1..n]]
    -- Present results as measurement statistics on qubit 0
    let measured = map (`measure` 0) results
        total0 = sum [p | [('0',p),('1',_)] <- measured]
        total1 = sum [p | [('0',_),('1',p)] <- measured]
        norm = total0 + total1
    putStrLn "Measurement on qubit 0:"
    putStrLn $ "'0': " ++ show (total0 / norm)
    putStrLn $ "'1': " ++ show (total1 / norm)

    putStrLn $ "Test: X with bit-flip error correction, with error X_p/3 ∘ Y_p/3 ∘ Z_p/3 with probability " ++ show p ++
               ", initial state " ++ show (roundMatrix 2 s) ++ ", " ++ show n ++ " iterations"
    let n_qubits = 3
        s_n = extendToNQubits n_qubits s
    results <- sequence [ do
        let cx1 = tensor_prod cx id_m
        s1 <- gateWithDepInQubit cx1 s_n p [0,1]
        let cx2 = matMul (matMul (extendSWAP 1 3) (tensor_prod cx id_m)) (extendSWAP 1 3)
        s2 <- gateWithDepInQubit cx2 s1 p [0,2]
        let xxx = tensor_prod (tensor_prod x x) x
        s3 <- gateWithDepInQubit xxx s2 p [0,1,2]
        let cx1 = tensor_prod cx id_m
        s4 <- gateWithDepInQubit cx1 s3 p [0,1]
        let cx2 = matMul (matMul (extendSWAP 1 3) (tensor_prod cx id_m)) (extendSWAP 1 3)
        s5 <- gateWithDepInQubit cx2 s4 p [0,2]
        let invccx = matMul (matMul (matMul (extendSWAP 0 3) (extendSWAP 1 3)) ccx) (matMul (extendSWAP 0 3) (extendSWAP 1 3))
        s6 <- gateWithDepInQubit invccx s5 p [0,1,2]
        return s6
        | _ <- [1..n]]
    let measured = map (`measure` 0) results
        total0 = sum [p | [('0',p),('1',_)] <- measured]
        total1 = sum [p | [('0',_),('1',p)] <- measured]
        norm = total0 + total1
    putStrLn "Measurement on qubit 0:"
    putStrLn $ "'0': " ++ show (total0 / norm)
    putStrLn $ "'1': " ++ show (total1 / norm)

runTest2Sweep :: Int -> String -> IO ()
runTest2Sweep nTest idOfTest = do
    let ps = [0.4,0.395..0.001]
        dir = "./data"
        fname = dir ++ "/out_test2_" ++ show nTest ++ "_" ++ idOfTest ++ ".csv"
    createDirectoryIfMissing True dir
    withFile fname WriteMode $ \h -> do
        hPutStrLn h "prob_error,1_no_corr,0_no_corr,1_corr,0_corr"
        mapM_ (\p -> do
            putStrLn $ "Running test2 for p = " ++ show p
            -- X without correction
            resultsNoCorr <- sequence [gateWithDepolarizing x excited_state_density p | _ <- [1..nTest]]
            let measuredNoCorr = map (`measure` 0) resultsNoCorr
                total0_no_corr = sum [v | [('0',v),('1',_)] <- measuredNoCorr]
                total1_no_corr = sum [v | [('0',_),('1',v)] <- measuredNoCorr]
                norm_no_corr = total0_no_corr + total1_no_corr
                p0_no_corr = if norm_no_corr == 0 then 0 else total0_no_corr / norm_no_corr
                p1_no_corr = if norm_no_corr == 0 then 0 else total1_no_corr / norm_no_corr
            -- X with correction
            let n_qubits = 3
                s_n = extendToNQubits n_qubits excited_state_density
            resultsCorr <- sequence [ do
                let cx1 = tensor_prod cx id_m
                s1 <- gateWithDepInQubit cx1 s_n (p*0.25) [0,1]
                let cx2 = matMul (matMul (extendSWAP 1 3) (tensor_prod cx id_m)) (extendSWAP 1 3)
                s2 <- gateWithDepInQubit cx2 s1 (p*0.25) [0,2]
                let xxx = tensor_prod (tensor_prod x x) x
                s3 <- gateWithDepInQubit xxx s2 p [0,1,2]
                let cx1 = tensor_prod cx id_m
                s4 <- gateWithDepInQubit cx1 s3 (p*0.25) [0,1]
                let cx2 = matMul (matMul (extendSWAP 1 3) (tensor_prod cx id_m)) (extendSWAP 1 3)
                s5 <- gateWithDepInQubit cx2 s4 (p*0.25) [0,2]
                let invccx = matMul (matMul (matMul (extendSWAP 0 3) (extendSWAP 1 3)) ccx) (matMul (extendSWAP 0 3) (extendSWAP 1 3))
                s6 <- gateWithDepInQubit invccx s5 (p*0.25) [0,1,2]
                return s6
                | _ <- [1..nTest]]
            let measuredCorr = map (`measure` 0) resultsCorr
                total0_corr = sum [v | [('0',v),('1',_)] <- measuredCorr]
                total1_corr = sum [v | [('0',_),('1',v)] <- measuredCorr]
                norm_corr = total0_corr + total1_corr
                p0_corr = if norm_corr == 0 then 0 else total0_corr / norm_corr
                p1_corr = if norm_corr == 0 then 0 else total1_corr / norm_corr
            hPutStrLn h $ printf "%.5f,%.8f,%.8f,%.8f,%.8f" p p1_no_corr p0_no_corr p1_corr p0_corr
            ) ps
    putStrLn $ "Results saved to " ++ fname

-- test 3 ---------------------------------------------------------------

test3 :: Double -> [[Complex Double]] -> Int -> IO ()
test3 p s n = do 
    putStrLn $ "Test: X, with error X_p/3 ∘ Y_p/3 ∘ Z_p/3 with probability " ++ show p ++
               ", initial state " ++ show (roundMatrix 2 s) ++ ", " ++ show n ++ " iterations"
    let expected = matMul (matMul x s) (dagger x)
    results <- sequence [gateWithDepolarizing x s p | _ <- [1..n]]
    -- Present results as measurement statistics on qubit 0
    let measured = map (`measure` 0) results
        total0 = sum [p | [('0',p),('1',_)] <- measured]
        total1 = sum [p | [('0',_),('1',p)] <- measured]
        norm = total0 + total1
    putStrLn "Measurement on qubit 0:"
    putStrLn $ "'0': " ++ show (total0 / norm)
    putStrLn $ "'1': " ++ show (total1 / norm)
    -- Fidelity
    let fidelities = [fidelity s1 expected | s1 <- results]
        avgFid = sum fidelities / fromIntegral n
    putStrLn $ "Average fidelity with expected state: " ++ show avgFid

-- Helper for fidelity
fidelity :: [[Complex Double]] -> [[Complex Double]] -> Double
fidelity rho sigma = realPart $ traceM $ matMul rho sigma

main :: IO()
main = do 
    args <- getArgs
    case args of
      ["test1", pStr, nStr] -> 
        let p = read pStr
            n = read nStr
        in test1 p excited_state_density n
      ["test2", pStr, nStr] -> 
        let p = read pStr
            n = read nStr
        in test2 p excited_state_density n
      ["test3", pStr, nStr] -> 
        let p = read pStr
            n = read nStr
        in test3 p excited_state_density n
      ["runTest2Sweep", nStr, idStr] ->
        let n = read nStr
        in runTest2Sweep n idStr
      _ -> putStrLn "Usage:\n  test1 <prob> <n>\n  test2 <prob> <n>\n  test3 <prob> <n>\n  runTest2Sweep <nTest> <idOfTest>"

