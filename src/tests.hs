module Main where

import Matrices
import ProbabilisticCombinator
import DensityMatrices
import Utils
import Distance
import Quantamorphism
import Noise
import MQfor
import qualified HCore as HC
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
import Data.Functor.Identity (runIdentity)

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
    putStrLn $ "Test: X, with error X_p/3 âˆ˜ Y_p/3 âˆ˜ Z_p/3 with probability " ++ show p ++
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
    putStrLn $ "Test: HZH, with error X_p/3 âˆ˜ Y_p/3 âˆ˜ Z_p/3 with probability " ++ show p ++
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
    putStrLn $ "Test: X, with error X_p/3 âˆ˜ Y_p/3 âˆ˜ Z_p/3 with probability " ++ show p ++
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

    putStrLn $ "Test: X with bit-flip error correction, with error X_p/3 âˆ˜ Y_p/3 âˆ˜ Z_p/3 with probability " ++ show p ++
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
    putStrLn $ "Test: X, with error X_p/3 âˆ˜ Y_p/3 âˆ˜ Z_p/3 with probability " ++ show p ++
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

-- test 5a ----------------------------------------------------------------
-- test quantamorphism with noise
-- target qubit is qubit 0, the control qubits are qubit 1 and 2
-- qubit 3 is an ancilla qubit.
-- M is Z gate and the number of controls is k = 2
-- noise probability of each gate has phase flip is 0.05%
--------------------------------------------------------------------------

test5a :: IO()
test5a = do
    putStrLn "Test 5a: Quantummorphism with noise, without error correction"
    let p = 0.05
    let n_qubits = 6
    let s_n = extendToNQubits n_qubits projPlus
    -- Get current date/time for unique file id
    t <- getCurrentTime
    let timeStr = formatTime defaultTimeLocale "%Y%m%d%H%M%S" t
        fname = "./data/out_test5a_" ++ show p ++ "_" ++ timeStr ++ ".csv"
    createDirectoryIfMissing True "./data"
    -- initial state |+âŸ©|0âŸ©|0âŸ©|0âŸ©
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

test5b :: IO()
test5b = do
    putStrLn "Test 5b: Quantummorphism with noise and error correction"
    let p = 0.05   
    let n_qubits = 6
    let s_n = extendToNQubits n_qubits projPlus
    -- Get current date/time for unique file id
    t <- getCurrentTime
    let timeStr = formatTime defaultTimeLocale "%Y%m%d%H%M%S" t
        fname = "./data/out_test5b_" ++ show p ++ "_" ++ timeStr ++ ".csv"
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

-- test 5c ----------------------------------------------------------------
-- No correction ancillas are used. Initial state is |1>|1>|1>|0|.
-- The measured target is qubit 0.
--------------------------------------------------------------------------
test5c :: IO ()
test5c = do
    let probabilities = [0.0005, 0.005, 0.05]
        targetProbability p = p * 0.1
        n_qubits = 4
        trials = 100
        initialState = excite_qubits [0,1,2] n_qubits
    t <- getCurrentTime
    let timeStr = formatTime defaultTimeLocale "%Y%m%d%H%M%S" t
        fname = "./data/test5c_raw_Combined_1110_" ++ timeStr ++ ".csv"
    createDirectoryIfMissing True "./data"
    putStrLn "Test 5c: combined reset/Z noise (p/2 reset, p/2 Z)"
    withFile fname WriteMode $ \h -> do
        hPutStrLn h "noise_model,probability,low_error_qubit,low_error_probability,trial,fidelity_target_qubit_0"
        mapM_ (\p -> mapM_ (\lowErrorQubit -> do
            ideal <- quantamorphism_b0_n2_test5c 0 0 Combined 0 initialState
            mapM_ (\trial -> do
                result <- quantamorphism_b0_n2_test5c p (targetProbability p) Combined lowErrorQubit initialState
                let fidelity = fidelityQubits result ideal [0] n_qubits
                hPutStrLn h $ intercalate "," [show Combined, show p, show lowErrorQubit,
                    show (targetProbability p), show trial, show fidelity]
                ) [1..trials]
            ) [0,1,2]
            ) probabilities
    putStrLn $ "Raw Test 5 results saved to " ++ fname
        
-- test 3a --------------------------------------------------------------
-- Builds the matrix for qfor H over labels (n, Bool), n=0..3, then applies
-- it to a random density matrix rho_0. Finally, applies SPAM bit-flip noise
-- to rho_0 with quantumChoice, runs the same circuit 100 times, and reports
-- the average distance from the clean final density matrix.
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
        modelHeader = "model\tfull_fidelity\ttarget_fidelity\tcontrol_fidelity"
        modelRow (name, fullFid, targetFid, controlFid) =
            printf "%s\t%.6f\t%.6f\t%.6f"
                name fullFid targetFid controlFid

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
            return (name, avgFull, avgTarget, avgControl)
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
            let rho1 = runIdentity (quantumChoiceMix targetQ0Flip id8 targetProb rho)
                rho2 = runIdentity (quantumChoiceMix controlQ1Flip id8 control1Prob rho1)
            in runIdentity (quantumChoiceMix controlQ2Flip id8 control2Prob rho2)

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

-- test 4a ----------------------------------------------------------------
-- Same qfor H circuit as test3a/test3b, but instead of a single SPAM
-- bit-flip applied once before the circuit, a full depolarizing channel
-- (X, Y, Z with equal probability p/3 each) is applied after EVERY H gate
-- firing (the qfor loop fires the gate n times, where n is the control
-- register value). Only the target qubit (q0, the Bool component of the
-- label) is affected by the noise; on a multi-qubit gate only its target
-- would be hit (H here is single-qubit, so it is simply the target).
-- Monte Carlo: mqfor with the IO-sampled per-gate depolarizing step,
-- averaged over trials against the noiseless H_qfor result (average
-- fidelity over trials, the same way test3a averages over quantumChoice
-- trials).
--------------------------------------------------------------------------
test4a :: IO ()
test4a = do
    let maxN = 3
        p = 0.1
        trials = 100
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
        basisDensity label =
            let v = vectorOf label
            in [[a * conjugate b | b <- v] | a <- v]

        -- H acting on the target qubit (index 2, LSB) tensored with Id on the two control qubits
        hTargetGate = HC.fromHMatrix $ HC.tensorListFromLists [id_m, id_m, h]
        idTargetGate = HC.fromHMatrix $ HC.tensorListFromLists [id_m, id_m, id_m]
        xTargetGate = HC.fromHMatrix $ HC.tensorListFromLists [id_m, id_m, x]
        yTargetGate = HC.fromHMatrix $ HC.tensorListFromLists [id_m, id_m, y]
        zTargetGate = HC.fromHMatrix $ HC.tensorListFromLists [id_m, id_m, z]

        stepWithQuantumChoice rho = do
            let afterGate = matMul (matMul hTargetGate rho) (dagger hTargetGate)
            afterX <- quantumChoice xTargetGate idTargetGate (p / 3) afterGate
            afterY <- quantumChoice yTargetGate idTargetGate (p / 3) afterX
            quantumChoice zTargetGate idTargetGate (p / 3) afterY

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
        labeledMatrixText title matrix = title ++ "\n" ++ unlines
            (("\t" ++ intercalate "\t" (map labelText labels)) :
            [ labelText label ++ "\t" ++ formatRow row
            | (label, row) <- zip labels matrix
            ])

        models = [ ("n1_single_gate", 1), ("n2_two_gates", 2), ("n3_three_gates", 3) ]
        modelHeader = "model\tn_gate_firings\tfull_fidelity\ttarget_fidelity\tcontrol_fidelity"
        modelRow (name, n, fullFid, targetFid, controlFid) =
            printf "%s\t%d\t%.6f\t%.6f\t%.6f"
                name (n :: Int) fullFid targetFid controlFid

    columns <- mapM
        (\(n, b) -> do
            (_, outTarget) <- mqfor applyH (n, targetVector b)
            return (vectorToBlock n outTarget)
        )
        labels

    let hQforMatrix = transpose columns

    modelResultsAndMatrices <- sequence
        [ do
            let rho0 = basisDensity (n, False)
                cleanFinal = matMul (matMul hQforMatrix rho0) (dagger hQforMatrix)
            measurements <- sequence
                [ do
                    (_, noisyFinal) <- mqfor stepWithQuantumChoice (n, rho0)
                    let fullFid = fidelity_density cleanFinal noisyFinal
                        targetFid = fidelity_density (targetReduced cleanFinal) (targetReduced noisyFinal)
                        controlFid = fidelity_density (controlReduced cleanFinal) (controlReduced noisyFinal)
                    return (fullFid, targetFid, controlFid)
                | _ <- [1 .. trials]
                ]
            let (fullFids, targetFids, controlFids) = unzip3 measurements
                avgFull = sum fullFids / fromIntegral trials
                avgTarget = sum targetFids / fromIntegral trials
                avgControl = sum controlFids / fromIntegral trials
                matricesText = unlines
                    [ "Model: " ++ name ++ " (n_gate_firings=" ++ show n ++ ")"
                    , labeledMatrixText "rho_0:" rho0
                    , labeledMatrixText "Clean (noiseless) final density matrix:" cleanFinal
                    ]
            return ((name, n, avgFull, avgTarget, avgControl), matricesText)
        | (name, n) <- models
        ]
    let (modelResults, modelMatrices) = unzip modelResultsAndMatrices

    let report = unlines
            [ "Test 4a: qfor H circuit with post-gate depolarizing noise (X,Y,Z equal prob p/3)"
            , "applied after EVERY gate firing (Monte Carlo, mqfor IO sampling)."
            , "Unlike test3a/test3b (SPAM applied once before the circuit), the noise here"
            , "is injected after each of the n H-gate firings driven by the qfor loop, and"
            , "only ever touches the gate's target qubit (H is single-qubit here, so that is"
            , "simply the target qubit q0 each time)."
            , "p=" ++ show p
            , "trials=" ++ show trials
            , ""
            , labeledMatrixText "H_qfor matrix:" hQforMatrix
            , intercalate "\n" modelMatrices
            , unlines (modelHeader : map modelRow modelResults)
            ]
        fname = "./data/out_test4a_qfor_h_depolarizing_mc.txt"

    createDirectoryIfMissing True "./data"
    writeFile fname report
    putStr report
    putStrLn $ "Test 4a results saved to " ++ fname

-- test 4b ----------------------------------------------------------------
-- Same setup as test4a, but the post-gate depolarizing noise is combined
-- exactly with the quantumChoiceMix implementation of the dissertation's
-- nested ((Z_(1/2) Diamond Y)_(2/3) Diamond X)_p Diamond I model.
-- No trials are needed because quantumChoiceMix returns the mixed density matrix.
--------------------------------------------------------------------------
test4b :: IO ()
test4b = do
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
        basisDensity label =
            let v = vectorOf label
            in [[a * conjugate b | b <- v] | a <- v]

        hTargetGate = HC.fromHMatrix $ HC.tensorListFromLists [id_m, id_m, h]
        idTargetGate = HC.fromHMatrix $ HC.tensorListFromLists [id_m, id_m, id_m]
        xTargetGate = HC.fromHMatrix $ HC.tensorListFromLists [id_m, id_m, x]
        yTargetGate = HC.fromHMatrix $ HC.tensorListFromLists [id_m, id_m, y]
        zTargetGate = HC.fromHMatrix $ HC.tensorListFromLists [id_m, id_m, z]

        stepWithQuantumChoiceMix rho =
            let afterGate = matMul (matMul hTargetGate rho) (dagger hTargetGate)
                zyChoice = runIdentity
                    (quantumChoiceMix zTargetGate yTargetGate (1 / 2) afterGate)
                xChoice = runIdentity
                    (quantumChoiceMix xTargetGate idTargetGate 1 afterGate)
                noisyChoice = matAdd
                    (scalarMul (2 / 3) zyChoice)
                    (scalarMul (1 / 3) xChoice)
            in matAdd
                (scalarMul p noisyChoice)
                (scalarMul (1 - p) afterGate)

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

        models = [ ("n1_single_gate", 1), ("n2_two_gates", 2), ("n3_three_gates", 3) ]

    columns <- mapM
        (\(n, b) -> do
            (_, outTarget) <- mqfor applyH (n, targetVector b)
            return (vectorToBlock n outTarget)
        )
        labels

    let hQforMatrix = transpose columns

        modelReports =
            [ let rho0 = basisDensity (n, False)
                  cleanFinal = matMul (matMul hQforMatrix rho0) (dagger hQforMatrix)
                  (_, noisyFinal) = qfor stepWithQuantumChoiceMix (n, rho0)
                  fid = fidelity_density cleanFinal noisyFinal
              in unlines
                    [ "Model: " ++ name ++ " (n_gate_firings=" ++ show n ++ ", p=" ++ show p ++ ")"
                    , labeledMatrixText "Final density matrix after n noisy H firings (exact):" noisyFinal
                    , "fidelity_vs_ideal=" ++ show fid
                    ]
            | (name, n) <- models
            ]

    let report = unlines
            [ "Test 4b: qfor H circuit with exact post-gate depolarizing noise (quantumChoiceMix)"
            , "Depolarizing channel encoded as ((Z_(1/2) Diamond Y)_(2/3) Diamond X)_p Diamond I,"
            , "applied after each H-gate firing driven by the qfor loop."
            , "p=" ++ show p
            , ""
            , intercalate "\n" modelReports
            ]
        fname = "./data/out_test4b_qfor_h_depolarizing_exact.txt"

    createDirectoryIfMissing True "./data"
    writeFile fname report
    putStr report
    putStrLn $ "Test 4b results saved to " ++ fname

--------------------------------------------------------------------------
main :: IO()
main = do
        args <- getArgs
        case args of {
            ["test1", pStr, nStr] -> let p = read pStr; n = read nStr in test1 p excited_state_density n;
            ["test2a", pStr, nStr] -> let p = read pStr; n = read nStr in test2a p excited_state_density n;
            ["test2b", pStr, nStr] -> let p = read pStr; n = read nStr in test2b p excited_state_density n;
            ["runTest2Sweep", nStr, idStr, dpStr] -> let n = read nStr; dp = read dpStr :: Double in runTest2Sweep n idStr dp;
            ["test3a"] -> test3a;
            ["test3b"] -> test3b;
            ["test4a"] -> test4a;
            ["test4b"] -> test4b;
            ["test5a"] -> test5a;
            ["test5b"] -> test5b;
            ["test5c"] -> test5c;
            _ -> putStrLn "Usage:\n  test1 <prob> <n>\n  test2a <prob> <n>\n  test2b <prob> <n>\n  runTest2Sweep <nTest> <idOfTest> <probabilityDecrease>\n  test3a\n  test3b\n  test4a\n  test4b\n  test5a\n  test5b\n  test5c"
        }
