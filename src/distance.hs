module Distance (
    fidelity,
    measure,
    fidelity_density,
    fidelityQubits
) where

import qualified Numeric.LinearAlgebra as LA
import Data.Complex

import Data.Complex
import Matrices 
import Matrices (sqrtm)
import Utils 
------------------------------------------------------

fidelity :: [[Complex Double]] -> [[Complex Double]] -> Double
fidelity rho sigma = realPart $ traceM $ matMul rho sigma

fidelity_density :: [[Complex Double]] -> [[Complex Double]] -> Double
fidelity_density rho sigma = realPart $ traceM $ matMul (sqrtm rho) (matMul sigma (sqrtm rho))

measure :: [[Complex Double]] -> Int -> [(Char, Double)]
measure s q =
    let nQubits = round (logBase 2 (fromIntegral (length s)))
        p0 = realPart $ traceM $ matMul (buildProjector nQubits q 0) s
        p1 = realPart $ traceM $ matMul (buildProjector nQubits q 1) s
    in [('0', p0), ('1', p1)]

-- | Partial trace over all qubits except those in keepQubits.
--   Returns the reduced density matrix as a list of lists (matrix).
partialTraceQubits :: [[Complex Double]] -> [Int] -> Int -> [[Complex Double]]
partialTraceQubits rho keepQubits nQubits =
    let keepSet = [ q | q <- [0..nQubits-1], q `elem` keepQubits ]
        tracedSet = [ q | q <- [0..nQubits-1], q `notElem` keepQubits ]
        nKeep = length keepSet
        nTraced = length tracedSet
        dimKeep = 2^nKeep
        -- Convert index to bit list of given length (MSB-first)
        toBitsLen len x = reverse $ take len $ map (\i -> (x `div` 2^i) `mod` 2) [0..]
        -- Convert bit list (MSB-first) to index
        fromBits bits = sum [ b * 2^i | (i,b) <- zip [0..] (reverse bits) ]
        -- For each (i',j') in reduced space, sum over all traced-out indices
        reducedElem i' j' =
            let bitsKeepI = toBitsLen nKeep i'
                bitsKeepJ = toBitsLen nKeep j'
                tracedCombinations = sequence [ [0,1] | _ <- [1..nTraced] ]
                idxInKeep q = case lookup q (zip keepSet [0..]) of Just ix -> ix; _ -> error "idxInKeep"
                idxInTraced q = case lookup q (zip tracedSet [0..]) of Just ix -> ix; _ -> error "idxInTraced"
                buildFull bitsK bitsT =
                    let fullBits = [ if q `elem` keepSet
                                     then bitsK !! idxInKeep q
                                     else bitsT !! idxInTraced q
                                   | q <- [0..nQubits-1] ]
                    in fromBits fullBits
            in sum [ rho !! (buildFull bitsKeepI bitsT) !! (buildFull bitsKeepJ bitsT) | bitsT <- tracedCombinations ]
    in [ [ reducedElem i j | j <- [0..dimKeep-1] ] | i <- [0..dimKeep-1] ]

-- | Fidelity for a subset of qubits: computes the fidelity between two density matrices, after tracing out all other qubits.
fidelityQubits :: [[Complex Double]] -> [[Complex Double]] -> [Int] -> Int -> Double
fidelityQubits rho1 rho2 keepQubits nQubits =
    let r1 = partialTraceQubits rho1 keepQubits nQubits
        r2 = partialTraceQubits rho2 keepQubits nQubits
    in fidelity_density r1 r2