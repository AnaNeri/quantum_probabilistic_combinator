module MQfor (
    Dist(..),
    qfor,
    mqfor,
    collapse,
    averageMatrices,
    sampleDist,
    depolarizingStep,
    stepWithGate,
    stepWithGateIO,
    bitFlipStep,
    bitFlipStepMulti,
    stepWithGateMulti,
    stepWithGateMultiIO,
    thermalRelaxStep,
    resetQubit,
    stepWithGateThermal,
    stepWithGateThermalIO,
    basisIndex,
    hammingDistance,
    expectedHammingDistance
) where

import Data.Complex
import Data.Bits (testBit, setBit)
import Control.Monad (foldM)
import Matrices (Matrix, matMul, dagger, scalarMul, matAdd, x, y, z, id_m)
import qualified HCore as HC
import System.Random (randomRIO)

-- Finite weighted-branch (Giry-style) probability monad
newtype Dist a = Dist { runDist :: [(Double, a)] }

instance Functor Dist where
    fmap f (Dist xs) = Dist [(p, f a) | (p, a) <- xs]

instance Applicative Dist where
    pure x' = Dist [(1.0, x')]
    Dist fs <*> Dist xs = Dist [(pf * px, f x') | (pf, f) <- fs, (px, x') <- xs]

instance Monad Dist where
    return = pure
    Dist xs >>= k = Dist [(p * q, b) | (p, a) <- xs, (q, b) <- runDist (k a)]

-- Plain quantamorphism loop: apply f n times, threading the state
qfor :: (b -> b) -> (Int, b) -> (Int, b)
qfor f (0, b) = (0, b)
qfor f (n, b) = let (m, b') = qfor f (n - 1, f b) in (m + 1, b')

-- Monadic quantamorphism loop, polymorphic in the monad (m = Dist or m = IO)
mqfor :: Monad m => (b -> m b) -> (Int, b) -> m (Int, b)
mqfor f (0, b) = return (0, b)
mqfor f (n, b) = do
    b' <- f b
    (m, b'') <- mqfor f (n - 1, b')
    return (m + 1, b'')

-- Weighted sum of branch density matrices: recovers the exact mixed state
collapse :: Dist Matrix -> Matrix
collapse (Dist xs) = foldl1 matAdd [scalarMul p m | (p, m) <- xs]

-- Average a list of density matrices (Monte Carlo estimate of the mixed state)
averageMatrices :: [Matrix] -> Matrix
averageMatrices ms = scalarMul (1 / fromIntegral (length ms)) (foldl1 matAdd ms)

-- Weighted random draw: generalizes quantumChoice's binary choice to n branches
sampleDist :: Dist a -> IO a
sampleDist (Dist xs) = do
    r <- randomRIO (0, 1) :: IO Double
    return (pick r 0 xs)
  where
    pick _ _ [(_, a)] = a
    pick r acc ((p, a) : rest)
        | r <= acc + p = a
        | otherwise = pick r (acc + p) rest
    pick _ _ [] = error "sampleDist: empty distribution"

-- One depolarizing-channel application as a 4-branch Dist (Kraus terms I,X,Y,Z)
-- on the given qubits within an n-qubit system.
depolarizingStep :: [Int] -> Int -> Double -> Matrix -> Dist Matrix
depolarizingStep qubitsList n p rho = Dist
    [ (1 - p, rho)
    , (p / 3, conjugate' (buildOp x))
    , (p / 3, conjugate' (buildOp y))
    , (p / 3, conjugate' (buildOp z))
    ]
  where
    buildOp op = HC.fromHMatrix $ HC.tensorListFromLists [if i `elem` qubitsList then op else id_m | i <- [0 .. n - 1]]
    conjugate' op = matMul (matMul op rho) (dagger op)

-- Apply a fixed gate, then depolarizing noise on qubitsList, as a Dist step
stepWithGate :: Matrix -> [Int] -> Int -> Double -> Matrix -> Dist Matrix
stepWithGate gate qubitsList n p rho =
    depolarizingStep qubitsList n p (matMul (matMul gate rho) (dagger gate))

-- Same step, sampled once in IO (Monte Carlo trajectory of the same Dist model)
stepWithGateIO :: Matrix -> [Int] -> Int -> Double -> Matrix -> IO Matrix
stepWithGateIO gate qubitsList n p rho = sampleDist (stepWithGate gate qubitsList n p rho)

-- SPAM-style state-preparation bit-flip: X w.p. p, I w.p. 1-p, on qubitsList only
bitFlipStep :: [Int] -> Int -> Double -> Matrix -> Dist Matrix
bitFlipStep qubitsList n p rho = Dist
    [ (1 - p, rho)
    , (p, conjugate' (buildOp x))
    ]
  where
    buildOp op = HC.fromHMatrix $ HC.tensorListFromLists [if i `elem` qubitsList then op else id_m | i <- [0 .. n - 1]]
    conjugate' op = matMul (matMul op rho) (dagger op)

-- Apply a gate, then independent depolarizing noise per qubit (each with its
-- own probability). The Dist monad's bind composes the per-qubit branches
-- (via foldM), so probabilities multiply and all combinations are enumerated.
stepWithGateMulti :: Matrix -> [(Int, Double)] -> Int -> Matrix -> Dist Matrix
stepWithGateMulti gate qubitProbs n rho =
    foldM (\r (q, p) -> depolarizingStep [q] n p r) (matMul (matMul gate rho) (dagger gate)) qubitProbs

-- Same step, sampled once in IO
stepWithGateMultiIO :: Matrix -> [(Int, Double)] -> Int -> Matrix -> IO Matrix
stepWithGateMultiIO gate qubitProbs n rho = sampleDist (stepWithGateMulti gate qubitProbs n rho)

-- SPAM-style state-preparation bit-flip with independent per-qubit probabilities
bitFlipStepMulti :: [(Int, Double)] -> Int -> Matrix -> Dist Matrix
bitFlipStepMulti qubitProbs n rho =
    foldM (\r (q, p) -> bitFlipStep [q] n p r) rho qubitProbs

-- Thermal relaxation (T1) + dephasing (T2) channel for a single qubit q,
-- valid in the T2(q) <= T1(q) regime (Georgopoulos et al., 2021):
--   p_T1 = exp(-Tg/T1), p_T2 = exp(-Tg/T2)
--   p_reset = 1 - p_T1                          (relax to |0>)
--   p_Z     = (1 - p_reset) * (1 - p_T2/sqrt(p_T1))   (dephasing)
--   p_I     = 1 - p_Z - p_reset                 (no error)
-- Kraus operators: K_I = sqrt(p_I) Id, K_Z = sqrt(p_Z) Z; the reset branch
-- uses the two-operator reset-to-|0> channel {sqrt(p_reset)|0><0|,
-- sqrt(p_reset)|0><1|} (via resetQubit) so it is trace-preserving for a
-- qubit found in |1> too, not just |0>.
thermalRelaxStep :: Int -> Int -> Double -> Double -> Double -> Matrix -> Dist Matrix
thermalRelaxStep q n t1 t2 tg rho = Dist
    [ (pI, conjugate' (buildOp id_m))
    , (pZ, conjugate' (buildOp z))
    , (pReset, resetQubit q n rho)
    ]
  where
    pT1 = exp (-tg / t1)
    pT2 = exp (-tg / t2)
    pReset = 1 - pT1
    pZ = (1 - pReset) * (1 - pT2 / sqrt pT1)
    pI = 1 - pZ - pReset
    buildOp op = HC.fromHMatrix $ HC.tensorListFromLists [if i == q then op else id_m | i <- [0 .. n - 1]]
    conjugate' op = matMul (matMul op rho) (dagger op)

-- Deterministic "reset qubit q to |0>" channel: K0=|0><0|_q, K1=|0><1|_q
-- (tensored with Id elsewhere). Trace-preserving regardless of qubit q's
-- state (unlike a single |0><0| conjugation, which kills the |1> component).
resetQubit :: Int -> Int -> Matrix -> Matrix
resetQubit q n rho =
    [ [ if testBit i bitPos || testBit j bitPos
          then 0 :+ 0
          else (rho !! i !! j) + (rho !! setBit i bitPos !! setBit j bitPos)
      | j <- [0 .. dim - 1] ]
    | i <- [0 .. dim - 1] ]
  where
    dim = 2 ^ n
    bitPos = n - 1 - q

-- Apply a gate, then independent thermal relaxation/dephasing noise on each
-- qubit in qubitsList (e.g. both control(s) and target of a multi-qubit
-- gate, since all qubits involved decohere during the same gate duration Tg).
stepWithGateThermal :: Matrix -> [Int] -> Int -> Double -> Double -> Double -> Matrix -> Dist Matrix
stepWithGateThermal gate qubitsList n t1 t2 tg rho =
    foldM (\r q -> thermalRelaxStep q n t1 t2 tg r) (matMul (matMul gate rho) (dagger gate)) qubitsList

-- Same step, sampled once in IO
stepWithGateThermalIO :: Matrix -> [Int] -> Int -> Double -> Double -> Double -> Matrix -> IO Matrix
stepWithGateThermalIO gate qubitsList n t1 t2 tg rho = sampleDist (stepWithGateThermal gate qubitsList n t1 t2 tg rho)

-- Index (0..2^n-1) of the basis state carrying the diagonal weight of a
-- classical/diagonal density matrix (valid for bit-flip-only circuits,
-- where every branch stays a computational basis state, no superposition).
basisIndex :: Matrix -> Int
basisIndex rho = snd $ maximum [(realPart (rho !! i !! i), i) | i <- [0 .. length rho - 1]]

-- Hamming distance (number of differing qubits) between two basis indices over n qubits
hammingDistance :: Int -> Int -> Int -> Int
hammingDistance n a b = length [() | k <- [0 .. n - 1], testBit a k /= testBit b k]

-- Expected number of wrong qubits (Hamming distance to a fixed ideal basis
-- index), averaged over a Dist's branches. This reveals error propagation
-- that full-state fidelity cannot: fidelity treats any wrong branch as
-- equally orthogonal, regardless of how many qubits it got wrong.
expectedHammingDistance :: Int -> Int -> Dist Matrix -> Double
expectedHammingDistance n idealIdx (Dist xs) =
    sum [p * fromIntegral (hammingDistance n idealIdx (basisIndex m)) | (p, m) <- xs]
