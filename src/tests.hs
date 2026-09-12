module Main where

import Matrices
import ProbabilisticCombinator
import DensityMatrices
import Utils
import Distance
import Quantamorphism
import Noise
import MQfor
import Data.Complex
import Control.Monad (foldM)
import qualified Data.Map as Map
import System.IO (withFile, IOMode(WriteMode), hPutStrLn)
import System.Directory (createDirectoryIfMissing)
import Text.Printf (printf)
import System.Environment (getArgs)
import Data.Time (getCurrentTime, formatTime, defaultTimeLocale)
import Data.List (intercalate, transpose)
import Data.Bits (xor)

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

-- test 2a ----------------------------------------------------------------
-- X vs X with bit-flip error correction scheme
-- bit-flip is not enough for error correction of faulty gates
-- may be helpful if it is in noise transmission, 
-- where correction gates have less noise than the working system
--------------------------------------------------------------------------- 

test2a :: Double -> [[Complex Double]] -> Int -> IO ()
test2a p s n = do
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

-- test 2b --------------------------------------------------------------
-- Fidelity calculation
-----------------------------------------------------------------------

test2b :: Double -> [[Complex Double]] -> Int -> IO ()
test2b p s n = do 
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

-- test 5 ----------------------------------------------------------------
-- No correction ancillas are used. Initial state is |1>|1>|1>|0|.
-- The measured target is qubit 0.
--------------------------------------------------------------------------
test5 :: IO ()
test5 = do
    let probabilities = [0.0005, 0.005, 0.05]
        targetProbability p = p * 0.1
        n_qubits = 4
        trials = 100
        initialState = excite_qubits [0,1,2] n_qubits
    t <- getCurrentTime
    let timeStr = formatTime defaultTimeLocale "%Y%m%d%H%M%S" t
        fname = "./data/test5_raw_Combined_1110_" ++ timeStr ++ ".csv"
    createDirectoryIfMissing True "./data"
    putStrLn "Test 5: combined reset/Z noise (p/2 reset, p/2 Z)"
    withFile fname WriteMode $ \h -> do
        hPutStrLn h "noise_model,probability,low_error_qubit,low_error_probability,trial,fidelity_target_qubit_0"
        mapM_ (\p -> mapM_ (\lowErrorQubit -> do
            ideal <- quantamorphism_b0_n2_test5 0 0 Combined 0 initialState
            mapM_ (\trial -> do
                result <- quantamorphism_b0_n2_test5 p (targetProbability p) Combined lowErrorQubit initialState
                let fidelity = fidelityQubits result ideal [0] n_qubits
                hPutStrLn h $ intercalate "," [show Combined, show p, show lowErrorQubit,
                    show (targetProbability p), show trial, show fidelity]
                ) [1..trials]
            ) [0,1,2]
            ) probabilities
    putStrLn $ "Raw Test 5 results saved to " ++ fname
        
-- test 6 ----------------------------------------------------------------
-- mqfor with m=Dist (exact branch enumeration) vs m=IO (Monte Carlo
-- sampling of the same Dist model), comparing depolarizing noise applied
-- to the target qubit vs a control qubit of a Toffoli (CCX) gate.
-- Self-contained 3-qubit system: raw ccx has controls=qubits 0,1 and
-- target=qubit 2 (8x8 density matrix); controls must start at |1> to fire.
-- Both noisy results are compared against the ideal (noiseless, p=0)
-- quantamorphism to see which qubit's noise causes more error.
--------------------------------------------------------------------------
test6 :: IO ()
test6 = do
    let n_qubits = 3
        steps = 2
        trials = 1000
        probabilities = [0.0005, 0.005, 0.05]
        initialState = excite_qubits [0,1] n_qubits
        cases = [("target", [2]), ("control", [0])]
        idealRho = collapse (fmap snd (mqfor (stepWithGate ccx [0] n_qubits 0) (steps, initialState)))
    t <- getCurrentTime
    let timeStr = formatTime defaultTimeLocale "%Y%m%d%H%M%S" t
        fname = "./data/out_test6_" ++ timeStr ++ ".csv"
    createDirectoryIfMissing True "./data"
    putStrLn "Test 6: mqfor with Dist (exact) vs IO (Monte Carlo) depolarizing noise, compared against the ideal (noiseless) quantamorphism"
    withFile fname WriteMode $ \h -> do
        hPutStrLn h "case,noise_qubit,probability,steps,trials,fidelity_exact_vs_ideal,fidelity_mc_vs_ideal,fidelity_exact_vs_mc"
        mapM_ (\(caseName, qubits) ->
            mapM_ (\p -> do
                let exactRho = collapse (fmap snd (mqfor (stepWithGate ccx qubits n_qubits p) (steps, initialState)))
                mcResults <- sequence
                    [ snd <$> mqfor (stepWithGateIO ccx qubits n_qubits p) (steps, initialState)
                    | _ <- [1 .. trials] ]
                let mcRho = averageMatrices mcResults
                    fidExactIdeal = fidelity_density exactRho idealRho
                    fidMcIdeal = fidelity_density mcRho idealRho
                    fidExactMc = fidelity_density exactRho mcRho
                hPutStrLn h $ intercalate "," [caseName, show qubits, show p, show steps, show trials,
                    show fidExactIdeal, show fidMcIdeal, show fidExactMc]
                putStrLn $ caseName ++ " qubits=" ++ show qubits ++ " p=" ++ show p ++
                    " fidelity(exact vs ideal)=" ++ show fidExactIdeal ++
                    " fidelity(mc vs ideal)=" ++ show fidMcIdeal
                ) probabilities
            ) cases
    putStrLn $ "Test 6 results saved to " ++ fname

-- test 7 ----------------------------------------------------------------
-- SPAM (State Preparation And Measurement) noise: bit-flip X applied once
-- during state preparation only, before any circuit gates. No noise is
-- applied on any other step; measurement noise is not modeled since the
-- quantamorphism defers measurement (only state-prep SPAM is tested here).
-- Same target/control comparison and ideal reference as test6.
--------------------------------------------------------------------------
test7 :: IO ()
test7 = do
    let n_qubits = 3
        steps = 2
        trials = 1000
        probabilities = [0.0005, 0.005, 0.05]
        cleanInitialState = excite_qubits [0,1] n_qubits
        cases = [("target", [2]), ("control", [0])]
        applyCleanSteps rho0 = snd (qfor (\rho -> matMul (matMul ccx rho) (dagger ccx)) (steps, rho0))
        idealRho = applyCleanSteps cleanInitialState
    t <- getCurrentTime
    let timeStr = formatTime defaultTimeLocale "%Y%m%d%H%M%S" t
        fname = "./data/out_test7_spam_" ++ timeStr ++ ".csv"
    createDirectoryIfMissing True "./data"
    putStrLn "Test 7: SPAM noise (bit-flip in state preparation only), mqfor Dist (exact) vs IO (Monte Carlo), vs ideal"
    withFile fname WriteMode $ \h -> do
        hPutStrLn h "case,noise_qubit,probability,steps,trials,fidelity_exact_vs_ideal,fidelity_mc_vs_ideal,fidelity_exact_vs_mc"
        mapM_ (\(caseName, qubits) ->
            mapM_ (\p -> do
                let prepD = bitFlipStep qubits n_qubits p cleanInitialState
                    exactRho = collapse (fmap applyCleanSteps prepD)
                mcResults <- sequence
                    [ applyCleanSteps <$> sampleDist prepD
                    | _ <- [1 .. trials] ]
                let mcRho = averageMatrices mcResults
                    fidExactIdeal = fidelity_density exactRho idealRho
                    fidMcIdeal = fidelity_density mcRho idealRho
                    fidExactMc = fidelity_density exactRho mcRho
                hPutStrLn h $ intercalate "," [caseName, show qubits, show p, show steps, show trials,
                    show fidExactIdeal, show fidMcIdeal, show fidExactMc]
                putStrLn $ caseName ++ " qubits=" ++ show qubits ++ " p=" ++ show p ++
                    " fidelity(exact vs ideal)=" ++ show fidExactIdeal ++
                    " fidelity(mc vs ideal)=" ++ show fidMcIdeal
                ) probabilities
            ) cases
    putStrLn $ "Test 7 results saved to " ++ fname

-- test 8 ----------------------------------------------------------------
-- SPAM noise with independent per-qubit probabilities: two qubits at p,
-- one "protected" qubit at p*0.1. Compares which protected qubit gives a
-- more robust circuit: testA protects the target (qubit 2), testB protects
-- control 0 (qubit 0). Bit-flip SPAM only, applied once at state
-- preparation; rest of the circuit (2x CCX) is noiseless, as in test7.
--------------------------------------------------------------------------
test8 :: IO ()
test8 = do
    let n_qubits = 3
        steps = 2
        trials = 1000
        probabilities = [0.0005, 0.005, 0.05]
        cleanInitialState = excite_qubits [0,1] n_qubits
        applyCleanSteps rho0 = snd (qfor (\rho -> matMul (matMul ccx rho) (dagger ccx)) (steps, rho0))
        idealRho = applyCleanSteps cleanInitialState
        cases = [ ("testA_target_protected", \p -> [(0, p), (1, p), (2, p * 0.1)])
                , ("testB_control0_protected", \p -> [(0, p * 0.1), (1, p), (2, p)])
                ]
    t <- getCurrentTime
    let timeStr = formatTime defaultTimeLocale "%Y%m%d%H%M%S" t
        fname = "./data/out_test8_spam_multi_" ++ timeStr ++ ".csv"
    createDirectoryIfMissing True "./data"
    putStrLn "Test 8: SPAM noise, per-qubit probabilities p,p,p*0.1 - target-protected vs control1-protected robustness"
    withFile fname WriteMode $ \h -> do
        hPutStrLn h "case,qubit_probs,probability,steps,trials,fidelity_exact_vs_ideal,fidelity_mc_vs_ideal,fidelity_exact_vs_mc"
        mapM_ (\(caseName, mkProbs) ->
            mapM_ (\p -> do
                let qubitProbs = mkProbs p
                    prepD = bitFlipStepMulti qubitProbs n_qubits cleanInitialState
                    exactRho = collapse (fmap applyCleanSteps prepD)
                mcResults <- sequence
                    [ applyCleanSteps <$> sampleDist prepD
                    | _ <- [1 .. trials] ]
                let mcRho = averageMatrices mcResults
                    fidExactIdeal = fidelity_density exactRho idealRho
                    fidMcIdeal = fidelity_density mcRho idealRho
                    fidExactMc = fidelity_density exactRho mcRho
                hPutStrLn h $ intercalate "," [caseName, show qubitProbs, show p, show steps, show trials,
                    show fidExactIdeal, show fidMcIdeal, show fidExactMc]
                putStrLn $ caseName ++ " probs=" ++ show qubitProbs ++ " p=" ++ show p ++
                    " fidelity(exact vs ideal)=" ++ show fidExactIdeal ++
                    " fidelity(mc vs ideal)=" ++ show fidMcIdeal
                ) probabilities
            ) cases
    putStrLn $ "Test 8 results saved to " ++ fname

-- test 9 ----------------------------------------------------------------
-- Demonstrates error propagation/accumulation through a single CCX.
-- With only 1 gate application (odd firing count), a control-qubit SPAM
-- error prevents the CCX from firing, so BOTH the control (already wrong)
-- AND the target (never flipped as it should have been) end up wrong: a
-- 2-qubit mismatch vs ideal. A target-qubit SPAM error only ever affects
-- the target: a 1-qubit mismatch. This asymmetry is invisible at steps=2
-- (test6/test8) because two firings vs zero firings are both even parity.
--------------------------------------------------------------------------
test9 :: IO ()
test9 = do
    let n_qubits = 3
        steps = 1
        trials = 1000
        probabilities = [0.0005, 0.005, 0.05]
        cleanInitialState = excite_qubits [0,1] n_qubits
        applyCleanSteps rho0 = snd (qfor (\rho -> matMul (matMul ccx rho) (dagger ccx)) (steps, rho0))
        idealRho = applyCleanSteps cleanInitialState
        cases = [("target", [2]), ("control", [0])]
    t <- getCurrentTime
    let timeStr = formatTime defaultTimeLocale "%Y%m%d%H%M%S" t
        fname = "./data/out_test9_spam_steps1_" ++ timeStr ++ ".csv"
    createDirectoryIfMissing True "./data"
    putStrLn "Test 9: SPAM noise, single CCX (steps=1), target-only vs control-only error propagation"
    withFile fname WriteMode $ \h -> do
        hPutStrLn h "case,noise_qubit,probability,steps,trials,fidelity_exact_vs_ideal,fidelity_mc_vs_ideal,fidelity_exact_vs_mc,expected_hamming_distance"
        mapM_ (\(caseName, qubits) ->
            mapM_ (\p -> do
                let prepD = bitFlipStep qubits n_qubits p cleanInitialState
                    branchesAfterCircuit = fmap applyCleanSteps prepD
                    exactRho = collapse branchesAfterCircuit
                    idealIdx = basisIndex idealRho
                    expHamming = expectedHammingDistance n_qubits idealIdx branchesAfterCircuit
                mcResults <- sequence
                    [ applyCleanSteps <$> sampleDist prepD
                    | _ <- [1 .. trials] ]
                let mcRho = averageMatrices mcResults
                    fidExactIdeal = fidelity_density exactRho idealRho
                    fidMcIdeal = fidelity_density mcRho idealRho
                    fidExactMc = fidelity_density exactRho mcRho
                hPutStrLn h $ intercalate "," [caseName, show qubits, show p, show steps, show trials,
                    show fidExactIdeal, show fidMcIdeal, show fidExactMc, show expHamming]
                putStrLn $ caseName ++ " qubits=" ++ show qubits ++ " p=" ++ show p ++
                    " fidelity(exact vs ideal)=" ++ show fidExactIdeal ++
                    " expected wrong qubits=" ++ show expHamming
                ) probabilities
            ) cases
    putStrLn $ "Test 9 results saved to " ++ fname

-- test 10 ---------------------------------------------------------------
-- Same propagation check as test9, but with post-gate depolarizing noise
-- (X,Y,Z, as in test6) instead of pre-gate SPAM. Since this noise is
-- applied AFTER the single CCX fires, it cannot influence that firing, so
-- (unlike test9) target and control noise are expected to show the SAME
-- (symmetric, non-propagating) expected Hamming distance to ideal.
--------------------------------------------------------------------------
test10 :: IO ()
test10 = do
    let n_qubits = 3
        steps = 1
        trials = 1000
        probabilities = [0.0005, 0.005, 0.05]
        initialState = excite_qubits [0,1] n_qubits
        cases = [("target", [2]), ("control", [0])]
        idealRho = collapse (fmap snd (mqfor (stepWithGate ccx [0] n_qubits 0) (steps, initialState)))
        idealIdx = basisIndex idealRho
    t <- getCurrentTime
    let timeStr = formatTime defaultTimeLocale "%Y%m%d%H%M%S" t
        fname = "./data/out_test10_depolarizing_steps1_" ++ timeStr ++ ".csv"
    createDirectoryIfMissing True "./data"
    putStrLn "Test 10: post-gate depolarizing noise, single CCX (steps=1), target-only vs control-only error propagation"
    withFile fname WriteMode $ \h -> do
        hPutStrLn h "case,noise_qubit,probability,steps,trials,fidelity_exact_vs_ideal,fidelity_mc_vs_ideal,fidelity_exact_vs_mc,expected_hamming_distance"
        mapM_ (\(caseName, qubits) ->
            mapM_ (\p -> do
                let branchesAfterCircuit = fmap snd (mqfor (stepWithGate ccx qubits n_qubits p) (steps, initialState))
                    exactRho = collapse branchesAfterCircuit
                    expHamming = expectedHammingDistance n_qubits idealIdx branchesAfterCircuit
                mcResults <- sequence
                    [ snd <$> mqfor (stepWithGateIO ccx qubits n_qubits p) (steps, initialState)
                    | _ <- [1 .. trials] ]
                let mcRho = averageMatrices mcResults
                    fidExactIdeal = fidelity_density exactRho idealRho
                    fidMcIdeal = fidelity_density mcRho idealRho
                    fidExactMc = fidelity_density exactRho mcRho
                hPutStrLn h $ intercalate "," [caseName, show qubits, show p, show steps, show trials,
                    show fidExactIdeal, show fidMcIdeal, show fidExactMc, show expHamming]
                putStrLn $ caseName ++ " qubits=" ++ show qubits ++ " p=" ++ show p ++
                    " fidelity(exact vs ideal)=" ++ show fidExactIdeal ++
                    " expected wrong qubits=" ++ show expHamming
                ) probabilities
            ) cases
    putStrLn $ "Test 10 results saved to " ++ fname

-- test 11 ---------------------------------------------------------------
-- T1 (relaxation) + T2 (dephasing) noise model from the dissertation
-- (Georgopoulos et al., 2021), T2(q) <= T1(q) regime: Kraus operators
-- K_I = sqrt(pI) Id, K_Z = sqrt(pZ) Z, K_reset = sqrt(preset) |0><0|.
-- Applied to ALL THREE qubits (both controls and the target) at once,
-- after a single atomic CCX, since every qubit involved in a multi-qubit
-- gate decoheres for the same gate duration Tg. Sweeps Tg to show how
-- longer atomic-gate execution time increases the error.
--------------------------------------------------------------------------
test11 :: IO ()
test11 = do
    let n_qubits = 3
        steps = 1
        trials = 1000
        t1 = 50000.0 -- ns
        t2 = 70000.0 -- ns
        tgValues = [20.0, 100.0, 300.0] -- ns, atomic-CCX gate duration sweep
        initialState = excite_qubits [0,1] n_qubits
        qubits = [0,1,2]
        idealRho = snd (qfor (\rho -> matMul (matMul ccx rho) (dagger ccx)) (steps, initialState))
    t <- getCurrentTime
    let timeStr = formatTime defaultTimeLocale "%Y%m%d%H%M%S" t
        fname = "./data/out_test11_thermal_" ++ timeStr ++ ".csv"
    createDirectoryIfMissing True "./data"
    putStrLn "Test 11: T1/T2 relaxation+dephasing noise, atomic CCX, applied to control(s) and target simultaneously"
    withFile fname WriteMode $ \h -> do
        hPutStrLn h "t1_ns,t2_ns,tg_ns,steps,trials,fidelity_exact_vs_ideal,fidelity_mc_vs_ideal,fidelity_exact_vs_mc"
        mapM_ (\tg -> do
            let exactRho = collapse (stepWithGateThermal ccx qubits n_qubits t1 t2 tg initialState)
            mcResults <- sequence
                [ stepWithGateThermalIO ccx qubits n_qubits t1 t2 tg initialState
                | _ <- [1 .. trials] ]
            let mcRho = averageMatrices mcResults
                fidExactIdeal = fidelity_density exactRho idealRho
                fidMcIdeal = fidelity_density mcRho idealRho
                fidExactMc = fidelity_density exactRho mcRho
            hPutStrLn h $ intercalate "," [show t1, show t2, show tg, show steps, show trials,
                show fidExactIdeal, show fidMcIdeal, show fidExactMc]
            putStrLn $ "Tg=" ++ show tg ++ "ns fidelity(exact vs ideal)=" ++ show fidExactIdeal
            ) tgValues
    putStrLn $ "Test 11 results saved to " ++ fname

-- test 12 ---------------------------------------------------------------
-- CCX decomposed into the standard 6-CNOT + H/T/T-dagger circuit (Nielsen
-- & Chuang, Fig. 4.9), for hardware with no native Toffoli gate. T1/T2
-- noise (same model as test11) is applied after each of the 6 CNOTs, on
-- the two qubits it touches, using a shorter per-CX gate duration Tg_cx.
-- Single-qubit gates (H,T,T-dagger) are treated as noiseless here, a
-- simplifying assumption since two-qubit gates dominate hardware error
-- budgets. Compared against test11's atomic CCX (one noise application at
-- a longer Tg_ccx) to see whether decomposing into more, shorter gates
-- increases or decreases the total error relative to a native Toffoli.
--------------------------------------------------------------------------
test12 :: IO ()
test12 = do
    let n_qubits = 3
        trials = 1000
        t1 = 50000.0 -- ns
        t2 = 70000.0 -- ns
        tgCx = 40.0   -- ns, single CX gate duration
        tgCcx = 300.0 -- ns, atomic CCX gate duration (comparison point)
        initialState = excite_qubits [0,1] n_qubits
        -- controls = qubits 0,1 ; target = qubit 2
        cx01 = extendCX 0 n_qubits
        cx12 = extendCX 1 n_qubits
        swap12 = extendSWAP 1 n_qubits
        cx02 = matMul (matMul swap12 cx01) swap12
        hOn k = applySingleQubitGate h k n_qubits
        tGate = [[1 :+ 0, 0 :+ 0], [0 :+ 0, cis (pi / 4)]]
        tDaggerGate = [[1 :+ 0, 0 :+ 0], [0 :+ 0, cis (-(pi / 4))]]
        tOn k = applySingleQubitGate tGate k n_qubits
        tdgOn k = applySingleQubitGate tDaggerGate k n_qubits
        -- Nielsen & Chuang Toffoli decomposition, controls=0,1 target=2
        gateSeq =
            [ hOn 2, cx12, tdgOn 2, cx02, tOn 2, cx12, tdgOn 2, cx02
            , tOn 1, tOn 2, hOn 2, cx01, tOn 0, tdgOn 1, cx01
            ]
        cxAtIndex = [(1, [1,2]), (3, [0,2]), (5, [1,2]), (7, [0,2]), (11, [0,1]), (14, [0,1])]
        applyClean rho0 = foldl (\r g -> matMul (matMul g r) (dagger g)) rho0 gateSeq
        idealRho = matMul (matMul ccx initialState) (dagger ccx)
        decomposedIdealRho = applyClean initialState
        sanityFidelity = fidelity_density decomposedIdealRho idealRho
        applyDecomposedWithNoise :: Matrix -> Dist Matrix
        applyDecomposedWithNoise rho0 = foldM step rho0 (zip [0 ..] gateSeq)
          where
            step r (i, g) =
                let r' = matMul (matMul g r) (dagger g)
                in case lookup i cxAtIndex of
                    Just qs -> foldM (\rr q -> thermalRelaxStep q n_qubits t1 t2 tgCx rr) r' qs
                    Nothing -> return r'
    putStrLn $ "Test 12 sanity check (decomposition == raw CCX, noiseless): fidelity=" ++ show sanityFidelity
    t <- getCurrentTime
    let timeStr = formatTime defaultTimeLocale "%Y%m%d%H%M%S" t
        fname = "./data/out_test12_decomposed_vs_atomic_" ++ timeStr ++ ".csv"
    createDirectoryIfMissing True "./data"
    putStrLn "Test 12: CCX decomposed into 6 CNOTs (Tg_cx each) vs atomic CCX (Tg_ccx), T1/T2 noise"
    let decomposedDist = applyDecomposedWithNoise initialState
        exactDecomposedRho = collapse decomposedDist
        exactAtomicRho = collapse (stepWithGateThermal ccx [0,1,2] n_qubits t1 t2 tgCcx initialState)
    mcDecomposed <- sequence [ sampleDist decomposedDist | _ <- [1 .. trials] ]
    mcAtomic <- sequence [ stepWithGateThermalIO ccx [0,1,2] n_qubits t1 t2 tgCcx initialState | _ <- [1 .. trials] ]
    let mcDecomposedRho = averageMatrices mcDecomposed
        mcAtomicRho = averageMatrices mcAtomic
        fidDecomposedExact = fidelity_density exactDecomposedRho idealRho
        fidDecomposedMc = fidelity_density mcDecomposedRho idealRho
        fidAtomicExact = fidelity_density exactAtomicRho idealRho
        fidAtomicMc = fidelity_density mcAtomicRho idealRho
    withFile fname WriteMode $ \h -> do
        hPutStrLn h "circuit,num_two_qubit_gates,tg_ns,t1_ns,t2_ns,trials,fidelity_exact_vs_ideal,fidelity_mc_vs_ideal"
        hPutStrLn h $ intercalate "," ["decomposed_6cx", "6", show tgCx, show t1, show t2, show trials, show fidDecomposedExact, show fidDecomposedMc]
        hPutStrLn h $ intercalate "," ["atomic_ccx", "1", show tgCcx, show t1, show t2, show trials, show fidAtomicExact, show fidAtomicMc]
    putStrLn $ "decomposed (6xCX, Tg=" ++ show tgCx ++ "ns): fidelity(exact vs ideal)=" ++ show fidDecomposedExact
    putStrLn $ "atomic CCX  (1xCCX, Tg=" ++ show tgCcx ++ "ns): fidelity(exact vs ideal)=" ++ show fidAtomicExact
    putStrLn $ "Test 12 results saved to " ++ fname

-- test 13 ---------------------------------------------------------------
-- Isolates the T1/T2 relaxation+dephasing noise (test11's model) to ONLY
-- the target qubit or ONLY one control qubit (not all three at once, as
-- test11 does), on a single CCX (steps=1), to see which qubit's own
-- decoherence hurts the circuit more. Reports both fidelity vs ideal and
-- expected Hamming distance (number of wrong qubits), since fidelity alone
-- can hide propagation asymmetry (as found with test9/test10).
--------------------------------------------------------------------------
test13 :: IO ()
test13 = do
    let n_qubits = 3
        steps = 1
        trials = 1000
        t1 = 50000.0 -- ns
        t2 = 70000.0 -- ns
        tgValues = [20.0, 100.0, 300.0] -- ns
        initialState = excite_qubits [0,1] n_qubits
        cases = [("target", 2), ("control", 0)]
        idealRho = snd (qfor (\rho -> matMul (matMul ccx rho) (dagger ccx)) (steps, initialState))
        idealIdx = basisIndex idealRho
    t <- getCurrentTime
    let timeStr = formatTime defaultTimeLocale "%Y%m%d%H%M%S" t
        fname = "./data/out_test13_thermal_target_vs_control_" ++ timeStr ++ ".csv"
    createDirectoryIfMissing True "./data"
    putStrLn "Test 13: T1/T2 noise isolated to target-only vs control-only, single CCX (steps=1)"
    withFile fname WriteMode $ \h -> do
        hPutStrLn h "case,noise_qubit,tg_ns,t1_ns,t2_ns,trials,fidelity_exact_vs_ideal,fidelity_mc_vs_ideal,expected_hamming_distance"
        mapM_ (\(caseName, q) ->
            mapM_ (\tg -> do
                let branchesAfterGate = stepWithGateThermal ccx [q] n_qubits t1 t2 tg initialState
                    exactRho = collapse branchesAfterGate
                    expHamming = expectedHammingDistance n_qubits idealIdx branchesAfterGate
                mcResults <- sequence
                    [ stepWithGateThermalIO ccx [q] n_qubits t1 t2 tg initialState
                    | _ <- [1 .. trials] ]
                let mcRho = averageMatrices mcResults
                    fidExactIdeal = fidelity_density exactRho idealRho
                    fidMcIdeal = fidelity_density mcRho idealRho
                hPutStrLn h $ intercalate "," [caseName, show q, show tg, show t1, show t2, show trials,
                    show fidExactIdeal, show fidMcIdeal, show expHamming]
                putStrLn $ caseName ++ " qubit=" ++ show q ++ " Tg=" ++ show tg ++ "ns" ++
                    " fidelity(exact vs ideal)=" ++ show fidExactIdeal ++
                    " expected wrong qubits=" ++ show expHamming
                ) tgValues
            ) cases
    putStrLn $ "Test 13 results saved to " ++ fname

-- test 3a --------------------------------------------------------------
-- Builds the matrix for qfor H over labels (n, Bool), n=0..3, then applies
-- it to a random density matrix rho_0. Finally, applies SPAM bit-flip noise
-- to rho_0 with quantumChoice, runs the same circuit 100 times, and reports
-- the average infidelity from the clean final density matrix.
--------------------------------------------------------------------------
test3a :: IO ()
test3a = do
    let maxN = 3
        p = 0.1
        trials = 100
        labels = [(n, b) | n <- [0 .. maxN], b <- [False, True]]
        dim = length labels
        zero = 0 :+ 0
        one = 1 :+ 0

        targetVector False = [one, zero]
        targetVector True = [zero, one]

        applyH target = do
            let col = matMul h [[target !! 0], [target !! 1]]
            return [head (col !! 0), head (col !! 1)]

        vectorToBlock n target = concat
            [ if n' == n then target else [zero, zero]
            | n' <- [0 .. maxN]
            ]

        vectorOf label =
            [ if rowLabel == label then one else zero
            | rowLabel <- labels
            ]
        permutationMatrix move = transpose [vectorOf (move label) | label <- labels]
        id8 = permutationMatrix id
        targetQ0Flip = permutationMatrix (\(n, b) -> (n, not b))
        controlQ1Flip = permutationMatrix (\(n, b) -> (n `xor` 2, b))
        controlQ2Flip = permutationMatrix (\(n, b) -> (n `xor` 1, b))
        basisDensity label =
            let v = vectorOf label
            in [[a * conjugate b | b <- v] | a <- v]

        applySpamModel targetProb control1Prob control2Prob rho = do
            rho1 <- quantumChoice targetQ0Flip id8 targetProb rho
            rho2 <- quantumChoice controlQ1Flip id8 control1Prob rho1
            quantumChoice controlQ2Flip id8 control2Prob rho2

        targetReduced rho =
            [ [sum [rho !! (2 * n + b) !! (2 * n + b') | n <- [0 .. maxN]]
              | b' <- [0, 1] ]
            | b <- [0, 1] ]
        controlReduced rho =
            [ [sum [rho !! (2 * n + b) !! (2 * n' + b) | b <- [0, 1]]
              | n' <- [0 .. maxN] ]
            | n <- [0 .. maxN] ]

        formatComplex (r :+ i)
            | abs i < 1e-9 = printf "%.6f" r
            | otherwise = printf "%.6f%+.6fi" r i
        formatRow row = intercalate "\t" (map formatComplex row)
        labelText (n, b) = "(" ++ show n ++ "," ++ show b ++ ")"
        labeledMatrixText matrix = unlines $
            ("\t" ++ intercalate "\t" (map labelText labels)) :
            [ labelText label ++ "\t" ++ formatRow row
            | (label, row) <- zip labels matrix
            ]
        matrixText title matrix = title ++ "\n" ++ unlines (map formatRow matrix)
        models =
            [ ("target_only", p, 0, 0)
            , ("control_q1_only", 0, p, 0)
            , ("control_q2_only", 0, 0, p)
            ]
        modelHeader = "model\tfull_fidelity\ttarget_fidelity\tcontrol_fidelity\tfull_distance\ttarget_distance\tcontrol_distance"
        modelRow (name, fullFid, targetFid, controlFid, fullDist, targetDist, controlDist) =
            printf "%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f"
                name fullFid targetFid controlFid fullDist targetDist controlDist

    columns <- mapM
        (\(n, b) -> do
            (_, outTarget) <- mqfor applyH (n, targetVector b)
            return (vectorToBlock n outTarget)
        )
        labels

    let hQforMatrix = transpose columns

    let rho0Label = (1, False)
        rho0 = basisDensity rho0Label
        cleanFinal = matMul (matMul hQforMatrix rho0) (dagger hQforMatrix)
        cleanTarget = targetReduced cleanFinal
        cleanControl = controlReduced cleanFinal

    modelResults <- sequence
        [ do
            measurements <- sequence
                [ do
                    rhoWithSpam <- applySpamModel targetProb control1Prob control2Prob rho0
                    let noisyFinal = matMul (matMul hQforMatrix rhoWithSpam) (dagger hQforMatrix)
                        fullFid = fidelity_density cleanFinal noisyFinal
                        targetFid = fidelity_density cleanTarget (targetReduced noisyFinal)
                        controlFid = fidelity_density cleanControl (controlReduced noisyFinal)
                    return (fullFid, targetFid, controlFid)
                | _ <- [1 .. trials]
                ]
            let avg (fullFids, targetFids, controlFids) =
                    ( sum fullFids / fromIntegral trials
                    , sum targetFids / fromIntegral trials
                    , sum controlFids / fromIntegral trials
                    )
                (avgFull, avgTarget, avgControl) = avg (unzip3 measurements)
            return (name, avgFull, avgTarget, avgControl,
                1 - avgFull, 1 - avgTarget, 1 - avgControl)
        | (name, targetProb, control1Prob, control2Prob) <- models
        ]

    let report = unlines
            [ "Test 3a: qfor H matrix with labels"
            , "Target qubit q0 is the Bool component in labels (n, Bool)."
            , "Control qubits q1 and q2 are encoded by n in binary."
            , "rho_0 label=" ++ labelText rho0Label
            , "p_spam=" ++ show p
            , "trials=" ++ show trials
            , ""
            , labeledMatrixText hQforMatrix
            , matrixText "rho_0:" rho0
            , matrixText "Final density matrix H_qfor rho_0 H_qfor^dagger:" cleanFinal
            , "Isolated SPAM model comparison (100 trials):"
            , unlines (modelHeader : map modelRow modelResults)
            ]
        fname = "./data/out_test3a_qfor_h.txt"

    createDirectoryIfMissing True "./data"
    writeFile fname report
    putStr report
    putStrLn $ "Test 3a results saved to " ++ fname

-- test 3b --------------------------------------------------------------
-- Same qfor H circuit and SPAM gates as test3a, but the SPAM error is
-- applied with the deterministic quantumChoiceMix combinator instead of
-- the Monte Carlo quantumChoice: rho' = p*(u' rho u'^dagger) + (1-p)*(u rho u^dagger).
-- No trials are needed since the result is exact; we print the resulting
-- density matrix directly to inspect how the error spreads across entries.
--------------------------------------------------------------------------
test3b :: IO ()
test3b = do
    let maxN = 3
        p = 0.1
        labels = [(n, b) | n <- [0 .. maxN], b <- [False, True]]
        zero = 0 :+ 0
        one = 1 :+ 0

        targetVector False = [one, zero]
        targetVector True = [zero, one]

        applyH target = do
            let col = matMul h [[target !! 0], [target !! 1]]
            return [head (col !! 0), head (col !! 1)]

        vectorToBlock n target = concat
            [ if n' == n then target else [zero, zero]
            | n' <- [0 .. maxN]
            ]

        vectorOf label =
            [ if rowLabel == label then one else zero
            | rowLabel <- labels
            ]
        permutationMatrix move = transpose [vectorOf (move label) | label <- labels]
        id8 = permutationMatrix id
        targetQ0Flip = permutationMatrix (\(n, b) -> (n, not b))
        controlQ1Flip = permutationMatrix (\(n, b) -> (n `xor` 2, b))
        controlQ2Flip = permutationMatrix (\(n, b) -> (n `xor` 1, b))
        basisDensity label =
            let v = vectorOf label
            in [[a * conjugate b | b <- v] | a <- v]

        applySpamModelMix targetProb control1Prob control2Prob rho =
            let rho1 = quantumChoiceMix targetQ0Flip id8 targetProb rho
                rho2 = quantumChoiceMix controlQ1Flip id8 control1Prob rho1
            in quantumChoiceMix controlQ2Flip id8 control2Prob rho2

        formatComplex (r :+ i)
            | abs i < 1e-9 = printf "%.6f" r
            | otherwise = printf "%.6f%+.6fi" r i
        formatRow row = intercalate "\t" (map formatComplex row)
        labelText (n, b) = "(" ++ show n ++ "," ++ show b ++ ")"
        labeledMatrixText title matrix = title ++ "\n" ++ unlines
            (("\t" ++ intercalate "\t" (map labelText labels)) :
            [ labelText label ++ "\t" ++ formatRow row
            | (label, row) <- zip labels matrix
            ])
        matrixText title matrix = title ++ "\n" ++ unlines (map formatRow matrix)
        models =
            [ ("target_only", p, 0, 0)
            , ("control_q1_only", 0, p, 0)
            , ("control_q2_only", 0, 0, p)
            , ("all_spam", p, p, p)
            ]

    columns <- mapM
        (\(n, b) -> do
            (_, outTarget) <- mqfor applyH (n, targetVector b)
            return (vectorToBlock n outTarget)
        )
        labels

    let hQforMatrix = transpose columns

        rho0Label = (1, False)
        rho0 = basisDensity rho0Label
        cleanFinal = matMul (matMul hQforMatrix rho0) (dagger hQforMatrix)

        modelReports =
            [ let rhoWithSpam = applySpamModelMix targetProb control1Prob control2Prob rho0
                  noisyFinal = matMul (matMul hQforMatrix rhoWithSpam) (dagger hQforMatrix)
              in unlines
                    [ "Model: " ++ name
                        ++ " (target_p=" ++ show targetProb
                        ++ ", control_q1_p=" ++ show control1Prob
                        ++ ", control_q2_p=" ++ show control2Prob ++ ")"
                    , labeledMatrixText "rho_0 after deterministic SPAM mixture:" rhoWithSpam
                    , labeledMatrixText "Final density matrix H_qfor rho_spam H_qfor^dagger:" noisyFinal
                    ]
            | (name, targetProb, control1Prob, control2Prob) <- models
            ]

    let report = unlines
            [ "Test 3b: qfor H matrix with deterministic SPAM mixture (quantumChoiceMix)"
            , "rho' = p*(u' rho u'^dagger) + (1-p)*(u rho u^dagger), applied exactly (no sampling)."
            , "Target qubit q0 is the Bool component in labels (n, Bool)."
            , "Control qubits q1 and q2 are encoded by n in binary."
            , "rho_0 label=" ++ labelText rho0Label
            , "p_spam=" ++ show p
            , ""
            , labeledMatrixText "H_qfor matrix:" hQforMatrix
            , matrixText "rho_0:" rho0
            , matrixText "Final density matrix H_qfor rho_0 H_qfor^dagger (clean):" cleanFinal
            , intercalate "\n" modelReports
            ]
        fname = "./data/out_test3b_qfor_h_deterministic.txt"

    createDirectoryIfMissing True "./data"
    writeFile fname report
    putStr report
    putStrLn $ "Test 3b results saved to " ++ fname

--------------------------------------------------------------------------
main :: IO()
main = do
        args <- getArgs
        case args of {
            ["test1", pStr, nStr] -> let p = read pStr; n = read nStr in test1 p excited_state_density n;
            ["test2a", pStr, nStr] -> let p = read pStr; n = read nStr in test2a p excited_state_density n;
            ["test2b", pStr, nStr] -> let p = read pStr; n = read nStr in test2b p excited_state_density n;
            ["runTest2Sweep", nStr, idStr, dpStr] -> let n = read nStr; dp = read dpStr :: Double in runTest2Sweep n idStr dp;
            ["test4a"] -> test4a;
            ["test4b"] -> test4b;
            ["test5"] -> test5;
            ["test6"] -> test6;
            ["test7"] -> test7;
            ["test8"] -> test8;
            ["test9"] -> test9;
            ["test10"] -> test10;
            ["test11"] -> test11;
            ["test12"] -> test12;
            ["test13"] -> test13;
            ["test3a"] -> test3a;
            ["test3b"] -> test3b;
            _ -> putStrLn "Usage:\n  test1 <prob> <n>\n  test2a <prob> <n>\n  test2b <prob> <n>\n  runTest2Sweep <nTest> <idOfTest> <probabilityDecrease>\n  test3a\n  test3b\n  test4a\n  test4b\n  test5 [PhaseFlip|Reset]\n  test6\n  test7\n  test8\n  test9\n  test10\n  test11\n  test12\n  test13"
        }

