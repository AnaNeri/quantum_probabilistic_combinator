module Noise (
    gateWithDepolarizing,
    gateWithBFInQubit,
    gateWithPFInQubit,
    gateWithDepolarizingExact,
    gateWithBFInQubitExact
) where

import Data.Complex
import Matrices
import Utils
import qualified HCore as HC
import ProbabilisticCombinator
------------------------------------------------------


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

-- Helper: apply a single-qubit bit-flip channel to qubit k in an n-qubit system
applyBitFlipToQubit :: Int -> Double -> [[Complex Double]] -> IO [[Complex Double]]
applyBitFlipToQubit k p state = do
    -- Build X acting on qubit k (tensor with id_m elsewhere)
    let n = round (logBase 2 (fromIntegral (length state))) -- number of qubits
        opOnQubit op = HC.fromHMatrix $ HC.tensorListFromLists [if i == k then op else id_m | i <- [0..n-1]]
    s1 <- quantumChoice (opOnQubit x) id_m (p) state
    return s1
    
-- General function: apply gate, then bit-flip noise only to qubits in qubits List
gateWithBFInQubit :: [[Complex Double]] -> [[Complex Double]] -> Double -> [Int] -> IO [[Complex Double]]
gateWithBFInQubit gate state prob qubitsList = do
    let n = round (logBase 2 (fromIntegral (length state)))
        buildNoiseOp op = HC.fromHMatrix $ HC.tensorListFromLists [if i `elem` qubitsList then op else id_m | i <- [0..n-1]]
        idN = identityN n
    let s1 = matMul (matMul gate state) (dagger gate)
    s2 <- quantumChoice (buildNoiseOp x) idN (prob) s1
    return s2

-- General function: apply gate, then phase-flip noise only to qubits in qubits List
gateWithPFInQubit :: [[Complex Double]] -> [[Complex Double]] -> Double -> [Int] -> IO [[Complex Double]]
gateWithPFInQubit gate state prob qubitsList = do
    let n = round (logBase 2 (fromIntegral (length state)))
        buildNoiseOp op = HC.fromHMatrix $ HC.tensorListFromLists [if i `elem` qubitsList then op else id_m | i <- [0..n-1]]
        idN = identityN n
    let s1 = matMul (matMul gate state) (dagger gate)
    s2 <- quantumChoice (buildNoiseOp z) idN (prob) s1
    return s2

-- Exact (non-sampling) counterpart of gateWithDepolarizing: computes the
-- exact channel output as a weighted sum of density matrices.
gateWithDepolarizingExact :: [[Complex Double]] -> [[Complex Double]] -> Double -> [[Complex Double]]
gateWithDepolarizingExact gate state prob =
    let state' = matMul (matMul gate state) (dagger gate)
        s1 = quantumChoiceExact x id_m (prob/3) state'
        s2 = quantumChoiceExact y id_m (prob/3) s1
        s3 = quantumChoiceExact z id_m (prob/3) s2
    in s3

-- Exact (non-sampling) counterpart of gateWithBFInQubit.
gateWithBFInQubitExact :: [[Complex Double]] -> [[Complex Double]] -> Double -> [Int] -> [[Complex Double]]
gateWithBFInQubitExact gate state prob qubitsList =
    let n = round (logBase 2 (fromIntegral (length state)))
        buildNoiseOp op = HC.fromHMatrix $ HC.tensorListFromLists [if i `elem` qubitsList then op else id_m | i <- [0..n-1]]
        idN = identityN n
        s1 = matMul (matMul gate state) (dagger gate)
    in quantumChoiceExact (buildNoiseOp x) idN prob s1