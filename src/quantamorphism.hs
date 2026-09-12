module Quantamorphism (
    quantamorphism_b0_n2,
    quantamorphism_b0_n2_wcorrection,
    quantamorphism_b0_n2_test5c,
    NoiseModel(..)
) where


import Data.Complex
import Matrices
import Noise
import Utils
import System.Mem (performGC)
import Control.Monad (when, foldM)
import Text.Printf (printf)
import GHC.Stats (getRTSStatsEnabled, getRTSStats)
import qualified HCore as HC
import qualified Numeric.LinearAlgebra as LA
import System.Random (randomRIO)
import Control.DeepSeq (deepseq)

import Distance

data NoiseModel = PhaseFlip | Reset | Combined deriving (Eq, Show, Read)

-- Memory check: estimate bytes used by a matrix (approx.)
defaultMemoryLimitBytes :: Int
defaultMemoryLimitBytes = 200000000 -- 200 MB, adjust as needed

estimateMatrixBytes :: [[Complex Double]] -> Int
estimateMatrixBytes m =
    let rows = length m
        cols = if null m then 0 else length (head m)
    in rows * cols * 16 -- approx 16 bytes per Complex Double

checkMatrixMemory :: String -> [[Complex Double]] -> IO ()
checkMatrixMemory name m = do
    let rows = length m
        cols = if null m then 0 else length (head m)
        estBytes = estimateMatrixBytes m
    putStrLn $ "[MEMCHECK] " ++ name ++ ": " ++ show rows ++ "x" ++ show cols ++ " ≈ " ++ show estBytes ++ " bytes"
    when (estBytes > defaultMemoryLimitBytes) $ do
        putStrLn $ "[MEMCHECK] Warning: estimated memory (" ++ show estBytes ++ " bytes) exceeds limit (" ++ show defaultMemoryLimitBytes ++ " bytes). Running GC."
        performGC
        putStrLn "[MEMCHECK] Performed GC."

-- RTS stats helper. Requires running the program with RTS `-T` enabled.
printRTSStats :: String -> IO ()
printRTSStats label = do
        enabled <- getRTSStatsEnabled
        if enabled
            then do
                stats <- getRTSStats
                putStrLn $ "[RTS] " ++ label ++ ": " ++ show stats
            else
                putStrLn "[RTS] RTS stats not enabled; run with +RTS -T to enable."



-- quantamorphism by gates 
quantamorphism_b0_n2 :: Double -> [[Complex Double]] -> Int -> Int -> IO [[Complex Double]]
quantamorphism_b0_n2 p s n_qubits protected = do
    let s_h = HC.toHMatrix $ extendToNQubits n_qubits s
        -- If protected < 0, treat as 'no qubit protected' (i.e., protected = 0)
        protected' = if protected < 0 then 0 else protected
    r_h <- quantamorphism_b0_n2_h p s_h n_qubits protected'
    return (HC.fromHMatrix r_h)

-- HMatrix-native implementation to avoid conversions when callers already have HMatrix
quantamorphism_b0_n2_h :: Double -> HC.HMatrix -> Int -> Int -> IO HC.HMatrix
quantamorphism_b0_n2_h p s_n_h n_qubits protected = do
    quantamorphism_b0_n2_h_with_noise
        (\gate state qubits -> gateWithPFInQubitH gate state p qubits)
        s_n_h n_qubits protected

quantamorphism_b0_n2_h_with_noise :: (HC.HMatrix -> HC.HMatrix -> [Int] -> IO HC.HMatrix) -> HC.HMatrix -> Int -> Int -> IO HC.HMatrix
quantamorphism_b0_n2_h_with_noise applyNoise s_n_h n_qubits protected = do
    let s_n_b = HC.fromHMatrix s_n_h
    let one_h = HC.toHMatrix $ applySingleQubitGate id_m 0 n_qubits
    s1_h <- applyNoise one_h s_n_h [0]
    -- If protected == 0, also apply id (with error) to 4 and 5
    s1a_h <- if protected == 0
        then do
            s' <- applyNoise (HC.toHMatrix $ applySingleQubitGate id_m 4 n_qubits) s1_h [4]
            s'' <- applyNoise (HC.toHMatrix $ applySingleQubitGate id_m 5 n_qubits) s' [5]
            return s''
        else return s1_h

    let q1plusq2a_h = HC.tensorListFromLists [if i == 1 || i == 2 then x else id_m | i <- [0..n_qubits-1]]
    s2a_h <- applyNoise q1plusq2a_h s1a_h [1,2]

    let q1plusq2b_h = HC.toHMatrix $ applyMultiQubitGate ccx 1 n_qubits
    s2b_h <- applyNoise q1plusq2b_h s2a_h [1,2,3]

    let q1plusq2c_h = HC.toHMatrix $ applySingleQubitGate x 3 n_qubits
    s2c_h <- applyNoise q1plusq2c_h s2b_h [3]

    let q1plusq2d_core_h = foldl1 HC.matMulH
            ([ HC.toHMatrix (extendSWAP 2 n_qubits)
             , HC.toHMatrix (extendSWAP 1 n_qubits)
             , HC.toHMatrix (extendSWAP 0 n_qubits)
             ]
             ++ (if protected == 0 then [HC.toHMatrix (applyMultiQubitGate cz 3 n_qubits),
                                         HC.toHMatrix (extendSWAP 4 n_qubits),
                                         HC.toHMatrix (applyMultiQubitGate cz 4 n_qubits),
                                         HC.toHMatrix (extendSWAP 4 n_qubits)] else [])
             ++ [ HC.toHMatrix (applyMultiQubitGate cz 0 n_qubits)
                , HC.toHMatrix (extendSWAP 0 n_qubits)
                , HC.toHMatrix (extendSWAP 1 n_qubits)
                , HC.toHMatrix (extendSWAP 2 n_qubits)
                ]
            )

    let q1plusq2d = HC.fromHMatrix q1plusq2d_core_h
    let q1plusq2d_targets = if protected == 0 then [0,3,4,5] else [0,3]
    s2d_h <- applyNoise (HC.toHMatrix q1plusq2d) s2c_h q1plusq2d_targets

    s2e_h <- applyNoise q1plusq2c_h s2d_h [3]

    s2f_h <- applyNoise q1plusq2b_h s2e_h [1,2,3]

    s2g_h <- applyNoise q1plusq2a_h s2f_h [1,2]

    let q2_core_h = foldl1 HC.matMulH
            ([ HC.toHMatrix (extendSWAP 1 n_qubits)
            , HC.toHMatrix (extendSWAP 0 n_qubits)
            ]
            ++ (if protected == 0 then [ HC.toHMatrix (extendSWAP 3 n_qubits),
                                         HC.toHMatrix (applyMultiQubitGate cz 3 n_qubits),
                                         HC.toHMatrix (extendSWAP 3 n_qubits),

                                         HC.toHMatrix (extendSWAP 4 n_qubits),
                                         HC.toHMatrix (extendSWAP 3 n_qubits),
                                         HC.toHMatrix (applyMultiQubitGate cz 3 n_qubits),
                                         HC.toHMatrix (extendSWAP 3 n_qubits),
                                         HC.toHMatrix (extendSWAP 4 n_qubits)] else [])
            ++ [ HC.toHMatrix (applyMultiQubitGate cz 0 n_qubits)
               , HC.toHMatrix (extendSWAP 0 n_qubits)
               , HC.toHMatrix (extendSWAP 1 n_qubits)
            ]
            )
    let q2 = HC.fromHMatrix q2_core_h
    let q2_targets = if protected == 0 then [2,0,4,5] else [2,0]
    s3_h <- applyNoise (HC.toHMatrix q2) s2g_h q2_targets

    s4a_h <- applyNoise q1plusq2b_h s3_h [1,2,3]

    s4b_h <- applyNoise (HC.toHMatrix q1plusq2d) s4a_h q1plusq2d_targets

    s4c_h <- applyNoise q1plusq2b_h s4b_h [1,2,3]

    return s4c_h

-- Test 5 noise: each affected qubit has an independent fault probability.
applyTest5Noise :: Double -> Double -> NoiseModel -> Int -> HC.HMatrix -> HC.HMatrix -> [Int] -> IO HC.HMatrix
applyTest5Noise otherProb lowErrorProb noiseModel lowErrorQubit gate state qubits = do
    let noisyState = applyUnitaryH gate state
        n = round (logBase 2 (fromIntegral (LA.rows noisyState)))
        probability q = if q == lowErrorQubit then lowErrorProb else otherProb
        noiseOperator op q = HC.tensorListH
            [if i == q then op else HC.toHMatrix id_m | i <- [0..n-1]]
        reset0 = HC.toHMatrix [[1 :+ 0, 0 :+ 0], [0 :+ 0, 0 :+ 0]]
        reset1 = HC.toHMatrix [[0 :+ 0, 1 :+ 0], [0 :+ 0, 0 :+ 0]]
        applyKraus k rho = HC.matMulH (HC.matMulH k rho) (LA.tr' k)
        applyReset q rho = applyKraus (noiseOperator reset0 q) rho
                        + applyKraus (noiseOperator reset1 q) rho
        applyFault q rho = case noiseModel of
            PhaseFlip -> applyUnitaryH (noiseOperator (HC.toHMatrix z) q) rho
            Reset -> applyReset q rho
            Combined -> applyReset q rho
        applyOne rho q = do
            fault <- randomRIO (0, 1)
            let p = probability q
            if fault <= p / 2
                then case noiseModel of
                    Combined -> return (applyReset q rho)
                    _ -> return (applyFault q rho)
                else if noiseModel == Combined && fault <= p
                    then return (applyUnitaryH (noiseOperator (HC.toHMatrix z) q) rho)
                    else if noiseModel == Combined
                        then return rho
                        else if fault <= p
                            then return (applyFault q rho)
                            else return rho
    foldM applyOne noisyState qubits

quantamorphism_b0_n2_test5c :: Double -> Double -> NoiseModel -> Int -> [[Complex Double]] -> IO [[Complex Double]]
quantamorphism_b0_n2_test5c otherProb lowErrorProb noiseModel lowErrorQubit state_i = do
    let applyNoise gate state qubits =
            applyTest5Noise otherProb lowErrorProb noiseModel lowErrorQubit gate state qubits
        s_n_h = HC.toHMatrix state_i
        n_qubits = round (logBase 2 (fromIntegral (LA.rows s_n_h)))
    result <- quantamorphism_b0_n2_h_with_noise applyNoise s_n_h n_qubits 1
    return (HC.fromHMatrix result)

-- quantamorphism by gates with correction       
quantamorphism_b0_n2_wcorrection :: Double -> [[Complex Double]] -> Int -> IO [[Complex Double]]
quantamorphism_b0_n2_wcorrection prob state_i qubit_wcorrection = do
    -- quantamorphism with correction; debug prints with fidelity
    let n_c_qubits = 6
    let s_n_b = extendToNQubits n_c_qubits state_i
    let s_n_b_h = HC.toHMatrix s_n_b

    let gates_hadamard = [ if i == qubit_wcorrection then h
             else if i == 4 || i == 5 then h
             else id_m
           | i <- [0..n_c_qubits-1] ]

    let swapForward swaps = foldl
            (\acc i ->
                if i < (n_c_qubits - 1)
                    then let m_h = HC.toHMatrix (extendSWAP i n_c_qubits)
                         in HC.matMulH m_h acc
                    else acc
            )
            (HC.toHMatrix $ identityN n_c_qubits)
            swaps

    let swapBackward swaps = foldl
            (\acc i ->
                if i < (n_c_qubits - 1)
                    then let m_h = HC.toHMatrix (extendSWAP i n_c_qubits)
                         in HC.matMulH m_h acc
                    else acc
            )
            (HC.toHMatrix $ identityN n_c_qubits)
            (reverse swaps)

    -- Generalized encode1: CX from qubit_wcorrection to 4
    let swaps1 = [qubit_wcorrection+1 .. 3]  -- 3 is one before 4
    let encode1_h = HC.matMulH (swapBackward swaps1) (HC.matMulH (HC.toHMatrix (extendCX qubit_wcorrection n_c_qubits)) (swapForward swaps1))
    k2a_h <- gateWithPFInQubitH encode1_h s_n_b_h prob [qubit_wcorrection, 4]

    let swaps2 = [qubit_wcorrection+1 .. 4]  -- 4 is one before 5
    let encode2_h = HC.matMulH (swapBackward swaps2) (HC.matMulH (HC.toHMatrix (extendCX qubit_wcorrection n_c_qubits)) (swapForward swaps2))
    k2b_h <- gateWithPFInQubitH encode2_h k2a_h prob [qubit_wcorrection, 5]

    let encode3_h = HC.tensorListFromLists gates_hadamard
    k2c_h <- gateWithPFInQubitH encode3_h k2b_h prob [qubit_wcorrection, 4,5]

    -- init state preparation
    -- prepare operator: x on qubits 1 and 2, id_m elsewhere
    let prepare = tensor_prod id_m (tensor_prod x (tensor_prod x (tensor_prod id_m (tensor_prod id_m id_m))))
    let prepare_h = HC.toHMatrix prepare
    k1_h <- gateWithPFInQubitH prepare_h k2c_h prob [qubit_wcorrection]

    -- Force evaluation of k2c and run GC before the expensive inner call
    k2c_h `seq` performGC

    k3_h <- quantamorphism_b0_n2_h prob k2c_h n_c_qubits qubit_wcorrection
    -- Force evaluation of k3 and run GC to release temporaries
    k3_h `seq` performGC

    -- Decoding
    let decode1_h = encode3_h
    k4a_h <- gateWithPFInQubitH decode1_h k3_h prob [qubit_wcorrection, 4,5]

    let decode2_h = encode2_h
    k4b_h <- gateWithPFInQubitH decode2_h k4a_h prob [qubit_wcorrection, 5]

    let decode3_h = encode1_h
    k4c_h <- gateWithPFInQubitH decode3_h k4b_h prob [qubit_wcorrection, 4]

    k4c_h `seq` performGC
    -- correction step:
    let swapForCorrection = [qubit_wcorrection..5]
    let correction_h = HC.matMulH (swapBackward swapForCorrection) (HC.matMulH (HC.toHMatrix (extendCCX 3 n_c_qubits)) (swapForward swapForCorrection))
    k5 <- gateWithPFInQubitH correction_h k4c_h prob [qubit_wcorrection, 4, 5]

    return (HC.fromHMatrix k5)

-- Local HMatrix-aware helpers to avoid repeated list<->HMatrix conversions
applyUnitaryH :: HC.HMatrix -> HC.HMatrix -> HC.HMatrix
applyUnitaryH u rho = HC.matMulH (HC.matMulH u rho) (LA.tr' u)

quantumChoiceH :: HC.HMatrix -> HC.HMatrix -> Double -> HC.HMatrix -> IO HC.HMatrix
quantumChoiceH u' u p v = do
    -- Handle edge probabilities deterministically to avoid randomness when p is 0 or 1
    if p <= 0 then return $ applyUnitaryH u v
    else if p >= 1 then return $ applyUnitaryH u' v
    else do
        prob <- randomRIO (0,1)
        if prob <= p
            then return $ applyUnitaryH u' v
            else return $ applyUnitaryH u v

gateWithPFInQubitH :: HC.HMatrix -> HC.HMatrix -> Double -> [Int] -> IO HC.HMatrix
gateWithPFInQubitH gateH stateH prob qubitsList = do
    let n = round (logBase 2 (fromIntegral (LA.rows stateH)))
        buildNoiseOp op = HC.tensorListFromLists [if i `elem` qubitsList then op else id_m | i <- [0..n-1]]
        idN = HC.tensorListFromLists (replicate n id_m)
    let s1 = applyUnitaryH gateH stateH
    -- If probability is zero, skip noise deterministically
    if prob <= 0 then return s1
    else quantumChoiceH (buildNoiseOp z) idN prob s1
