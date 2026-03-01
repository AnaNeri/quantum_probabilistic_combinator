module Utils (
    tensor_prod,
    identityN,
    extendToNQubits,
    excite_qubits,
    extendSWAP,
    extendCX,
    buildProjector, 
    applyMultiQubitGate,
    applySingleQubitGate,
    extendCCX
) where

import Data.Complex
import Matrices 
import DensityMatrices 
import qualified HCore as HC

------------------------------------------------------

-- n-qubit identity matrix
identityN :: Int -> [[Complex Double]]
identityN n = HC.fromHMatrix $ HC.tensorListFromLists (replicate n id_m)

-- Extend a single-qubit state to n qubits (all others in |0⟩)
extendToNQubits :: Int -> [[Complex Double]] -> [[Complex Double]]
extendToNQubits n s
    | n <= 1    = s
    | otherwise = HC.fromHMatrix $ HC.tensorListFromLists (s : replicate (n-1) ground_state_density)

-- Helper to apply a single-qubit gate at position k (0-based, leftmost is 0)
applySingleQubitGate :: [[Complex Double]] -> Int -> Int -> [[Complex Double]]
applySingleQubitGate gate k n_qubits = HC.fromHMatrix $ HC.tensorListFromLists [if i == k then gate else id_m | i <- [0..n_qubits-1]]

-- Helper to apply a multi-qubit gate at position k (replaces k..k+gate_qubits-1)
applyMultiQubitGate :: [[Complex Double]] -> Int -> Int -> [[Complex Double]]
applyMultiQubitGate gate k n_qubits =
    let gate_qubits = round (logBase 2 (fromIntegral (length gate)))
        before = replicate k id_m
        after = replicate (n_qubits - k - gate_qubits) id_m
    in HC.fromHMatrix $ HC.tensorListFromLists (before ++ [gate] ++ after)


-- Excite qubits [k] in an n-qubit system
excite_qubits :: [Int] -> Int -> [[Complex Double]]
excite_qubits qubitsList nQubits =
    HC.fromHMatrix $ HC.tensorListFromLists [if i `elem` qubitsList then excited_state_density else ground_state_density | i <- [0..nQubits-1]]

-- Extend SWAP gate to nQubits, swapping qubits n and n+1
extendSWAP :: Int -> Int -> [[Complex Double]]
extendSWAP n nQubits
        | n < 0 || n+1 >= nQubits = error "extendSWAP: invalid qubit indices"
        | otherwise = HC.fromHMatrix $ HC.tensorListFromLists ops
  where
    ops = buildOps 0
    buildOps i
      | i >= nQubits = []
      | i == n       = swap : buildOps (i+2)  -- place swap at position n (acts on n and n+1), skip n+1
      | otherwise    = id_m : buildOps (i+1)  -- identity elsewhere

-- Extend cx gate to nQubits, with the contol in position c and target in c+1
extendCX :: Int -> Int -> [[Complex Double]]
extendCX c nQubits
        | c < 0 || c+1 >= nQubits = error "extendCX: invalid qubit indices"
        | otherwise = HC.fromHMatrix $ HC.tensorListFromLists ops
  where
    ops = buildOps 0
    buildOps i
      | i >= nQubits = []
      | i == c       = cx : buildOps (i + 2)
      | otherwise    = id_m : buildOps (i + 1)

-- extend CCX gate to nQubits, with controls in positions c1 and c1+1 and target in c1+1
extendCCX :: Int -> Int -> [[Complex Double]]
extendCCX c1 nQubbits
        | c1 < 0 || c1+2 >= nQubbits = error "extendCCX: invalid qubit indices"
        | otherwise = HC.fromHMatrix $ HC.tensorListFromLists ops
  where
    ops = buildOps 0
    buildOps i
      | i >= nQubbits = []
      | i == c1       = ccx : buildOps (i + 3)
      | otherwise     = id_m : buildOps (i + 1)

-- Build the projector for measuring qubit q in n-qubit system, for outcome 0 or 1
buildProjector :: Int -> Int -> Int -> [[Complex Double]]
buildProjector nQubits q outcome =
    HC.fromHMatrix $ HC.tensorListFromLists [if i == q then (if outcome == 0 then proj0 else proj1) else id_m | i <- [0..nQubits-1]]