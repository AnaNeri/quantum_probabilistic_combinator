module DensityMatrices (
    proj0,
    proj1,
    projPlus,
    projMinus,
    ground_state_density,
    excited_state_density,
    -- hmatrix variants
    proj0H,
    proj1H
) where

import Data.Complex
import qualified HCore as HC
------------------------------------------------------

-- Projector for |0⟩ on a single qubit
-- hmatrix-backed projectors (preferred)
proj0H :: HC.HMatrix
proj0H = HC.toHMatrix [[1 :+ 0, 0 :+ 0], [0 :+ 0, 0 :+ 0]]

-- List-of-rows compatibility exports (existing API)
proj0 :: [[Complex Double]]
proj0 = HC.fromHMatrix proj0H

-- Projector for |1⟩ on a single qubit
proj1H :: HC.HMatrix
proj1H = HC.toHMatrix [[0 :+ 0, 0 :+ 0], [0 :+ 0, 1 :+ 0]]

proj1 :: [[Complex Double]]
proj1 = HC.fromHMatrix proj1H

-- Projector for |+>
projPlus :: [[Complex Double]]
projPlus = HC.fromHMatrix $ HC.toHMatrix [[0.5 :+ 0, 0.5 :+ 0],[0.5 :+ 0, 0.5 :+ 0]]

projMinus :: [[Complex Double]]
projMinus = HC.fromHMatrix $ HC.toHMatrix [[0.5 :+ 0, (-0.5) :+ 0],[(-0.5) :+ 0, 0.5 :+ 0]]

ground_state_density :: [[Complex Double]]
ground_state_density = HC.fromHMatrix $ HC.toHMatrix [[1 :+ 0, 0 :+ 0],[0 :+ 0, 0 :+ 0]]

excited_state_density :: [[Complex Double]]
excited_state_density = HC.fromHMatrix $ HC.toHMatrix [[0 :+ 0, 0 :+ 0],[0 :+ 0, 1 :+ 0]]