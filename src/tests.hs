module Main where

import Matrices
import ProbabilisticCombinator
import DensityMatrices
import Utils
import Distance
import Quantamorphism
import Noise
import Data.Complex
import qualified Data.Map as Map
import System.IO (withFile, IOMode(WriteMode), hPutStrLn)
import System.Directory (createDirectoryIfMissing)
import Text.Printf (printf)
import System.Environment (getArgs)
import Data.Time (getCurrentTime, formatTime, defaultTimeLocale)

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
        s1 <- gateWithBFInQubit cx1 s_n p [0,1]
        let cx2 = matMul (matMul (extendSWAP 1 3) (tensor_prod cx id_m)) (extendSWAP 1 3)
        s2 <- gateWithBFInQubit cx2 s1 p [0,2]
        let xxx = tensor_prod (tensor_prod x x) x
        s3 <- gateWithBFInQubit xxx s2 p [0,1,2]
        let cx1 = tensor_prod cx id_m
        s4 <- gateWithBFInQubit cx1 s3 p [0,1]
        let cx2 = matMul (matMul (extendSWAP 1 3) (tensor_prod cx id_m)) (extendSWAP 1 3)
        s5 <- gateWithBFInQubit cx2 s4 p [0,2]
        let invccx = matMul (matMul (matMul (extendSWAP 0 3) (extendSWAP 1 3)) ccx) (matMul (extendSWAP 0 3) (extendSWAP 1 3))
        s6 <- gateWithBFInQubit invccx s5 p [0,1,2]
        return s6
        | _ <- [1..n]]
    let measured = map (`measure` 0) results
        total0 = sum [p | [('0',p),('1',_)] <- measured]
        total1 = sum [p | [('0',_),('1',p)] <- measured]
        norm = total0 + total1
    putStrLn "Measurement on qubit 0:"
    putStrLn $ "'0': " ++ show (total0 / norm)
    putStrLn $ "'1': " ++ show (total1 / norm)

runTest2Sweep :: Int -> String -> Double -> IO ()
runTest2Sweep nTest idOfTest decreaseProb = do
    let ps = [0.4,0.395..0.001]
        dir = "./data"
        fname = dir ++ "/out_test2_" ++ show nTest ++ "_" ++ idOfTest ++ "_" ++ show decreaseProb ++ ".csv"
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
                s1 <- gateWithBFInQubit cx1 s_n (p*(1-decreaseProb)) [0,1]
                let cx2 = matMul (matMul (extendSWAP 1 3) (tensor_prod cx id_m)) (extendSWAP 1 3)
                s2 <- gateWithBFInQubit cx2 s1 (p*(1-decreaseProb)) [0,2]
                let xxx = tensor_prod (tensor_prod x x) x
                s3 <- gateWithBFInQubit xxx s2 p [0,1,2]
                let cx1 = tensor_prod cx id_m
                s4 <- gateWithBFInQubit cx1 s3 (p*(1-decreaseProb)) [0,1]
                let cx2 = matMul (matMul (extendSWAP 1 3) (tensor_prod cx id_m)) (extendSWAP 1 3)
                s5 <- gateWithBFInQubit cx2 s4 (p*(1-decreaseProb)) [0,2]
                let invccx = matMul (matMul (matMul (extendSWAP 0 3) (extendSWAP 1 3)) ccx) (matMul (extendSWAP 0 3) (extendSWAP 1 3))
                s6 <- gateWithBFInQubit invccx s5 (p*(1-decreaseProb)) [0,1,2]
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
-- Fidelity calculation
-----------------------------------------------------------------------

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

-- test 4 ----------------------------------------------------------------
-- test quantamorphism with noise
-- target qubit is qubit 0, the control qubits are qubit 1 and 2
-- qubit 3 is an ancilla qubit.
-- M is Z gate and the number of controls is k = 2
-- noise probability of each gate has phase flip is 0.05%
--------------------------------------------------------------------------

test4a :: IO()
test4a = do
    putStrLn "Test 4: Quantummorphism with noise, without error correction"
    let p = 0.05
    let n_qubits = 6
    let s_n = extendToNQubits n_qubits projPlus
    -- Get current date/time for unique file id
    t <- getCurrentTime
    let timeStr = formatTime defaultTimeLocale "%Y%m%d%H%M%S" t
        fname = "./data/out_test4a_" ++ show p ++ "_" ++ timeStr ++ ".csv"
    createDirectoryIfMissing True "./data"
    -- initial state |+⟩|0⟩|0⟩|0⟩
    res_ideal <- quantamorphism_b0_n2 0 projPlus n_qubits (-1)
    withFile fname WriteMode $ \h -> do
        hPutStrLn h "fidelity"
        fidelity_with_noise <- sequence [
            do
                res_with_noise <- quantamorphism_b0_n2 p projPlus n_qubits (-1)
                let fid_res_with_noise = fidelity_density res_with_noise res_ideal
                hPutStrLn h (show fid_res_with_noise)
                putStrLn $ "Fidelity of run with noise: " ++ show fid_res_with_noise
                return fid_res_with_noise
            | _ <- [1..100]
            ]
        return ()
    putStrLn $ "Fidelity results saved to " ++ fname

test4b :: IO()
test4b = do
    putStrLn "Test 4b: Quantummorphism with noise and error correction"
    let p = 0.05   
    let n_qubits = 6
    let s_n = extendToNQubits n_qubits projPlus
    -- Get current date/time for unique file id
    t <- getCurrentTime
    let timeStr = formatTime defaultTimeLocale "%Y%m%d%H%M%S" t
        fname = "./data/out_test4b_" ++ show p ++ "_" ++ timeStr ++ ".csv"
    createDirectoryIfMissing True "./data"
    putStrLn $ "Start ideal run"
    res_ideal <- quantamorphism_b0_n2 0 projPlus n_qubits 0
    withFile fname WriteMode $ \h -> do
        hPutStrLn h "fidelity_full,fidelity_partial"
        fidelity_with_correction <- sequence [
            do
                putStrLn $ "Start run with correction, noise probability: " ++ show p
                res_with_correction <- quantamorphism_b0_n2_wcorrection p projPlus 0
                let fid_full = fidelity_density res_with_correction res_ideal
                    fid_partial = fidelityQubits res_with_correction res_ideal [0] n_qubits
                hPutStrLn h $ show fid_full ++ "," ++ show fid_partial
                putStrLn $ "Fidelity of run with correction (full): " ++ show fid_full
                putStrLn $ "Fidelity of run with correction (qubits 0-3): " ++ show fid_partial
                return fid_partial
            | _ <- [1..100]
            ]
        return ()
    putStrLn $ "Fidelity results saved to " ++ fname
        
--------------------------------------------------------------------------
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
      ["runTest2Sweep", nStr, idStr, dpStr] ->
        let n = read nStr
            dp = read dpStr :: Double
        in runTest2Sweep n idStr dp
      ["test4a"] ->
        test4a
      ["test4b"] ->
        test4b
      _ -> putStrLn "Usage:\n  test1 <prob> <n>\n  test2 <prob> <n>\n  test3 <prob> <n>\n  runTest2Sweep <nTest> <idOfTest> <probabilityDecrease>\n  test4a \n  test4b"

