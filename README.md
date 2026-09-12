# Quantum Probabilistic Combinator

Tests of a quantum probabilistic combinator in Haskell: quantum gates and density-matrix
simulation are used to study how noise/error models and error-correction strategies affect
circuit fidelity. The core test suite is stable; further changes are mostly additions
explored out of curiosity or needed to support related publications.

## Requirements

* GHC and `cabal` (Haskell2010, tested with GHC 9.6.7)
* BLAS/LAPACK development libraries (required by the `hmatrix` dependency)
* Python 3.8+ with `pandas` and `matplotlib` for the analysis/plotting scripts in `data/`
    (dependencies declared in [pyproject.toml](pyproject.toml); e.g. `uv run --with pandas --with matplotlib --no-project data/analyze_test5c.py`)

## Build

```sh
cabal build
```

This builds the `QuantumProbabilisticCombinators` library and the `tests` executable
(see [quantum_probabilistic_combinators.cabal](quantum_probabilistic_combinators.cabal)).

## Modules

* [matrices.hs](src/matrices.hs) - quantum gates
    * Pauli matrices: x, y, z
    * hadamard gate
    * identity gate
    * cx gate
    * cu3 gate
* [probabilisticcombinator.hs](src/probabilisticcombinator.hs) - probabilistic/quantum choice combinators
* [noise.hs](src/noise.hs) - depolarizing / bit-flip / phase-flip noise channels
* [quantamorphism.hs](src/quantamorphism.hs) - quantamorphism (recursive circuit) construction
* [tests.hs](src/tests.hs) - test suite entry point, described below

## Running the tests

Run via `cabal run tests -- <args>`, e.g.:

```sh
cabal run tests -- test1 0.05 20
```

* test 1 — `test1 <prob> <n>`
    - check noise accumulation in circuit X vs circuit HZH
    - shows that more depth usually results in more noise
* test 2a — `test2a <prob> <n>` (also `runTest2Sweep <nTest> <idOfTest> <probabilityDecrease>`)
    - check noise in circuit X with and without error correction bit flip error (implementation uses DOI: [10.13140/RG.2.2.18542.77129](https://www.researchgate.net/publication/334634646_The_first_three-qubit_and_six-qubit_full_quantum_multiple_error-correcting_codes_with_low_quantum_costs?channel=doi&linkId=5d36ffe2a6fdcc370a57ac6a&showFulltext=true))
    - this type of codes is used in message transmission and does not seem to work in programs.
    - allow to see that error correction reduces error when correction gates are better than the working gates.
* test 2b — `test2b <prob> <n>`
    - check fidelity of the noisy X gate against the ideal expected state
* test 3a — `test3a`
    - qfor H matrix with labels: Monte Carlo SPAM bit-flip model over the encoded basis states
    - compares target-only vs control-only noise effects on the final density matrix using repeated sampling
    - writes the report to `data/out_test3a_qfor_h.txt`
* test 3b — `test3b`
    - same qfor H setup as test 3a, but uses the deterministic SPAM mixture `quantumChoiceMix`
    - applies the exact mixture instead of Monte Carlo sampling to inspect how the error spreads through the density matrix
    - writes the report to `data/out_test3b_qfor_h_deterministic.txt`
* test 4a — `test4a`
    - qfor H circuit with post-gate depolarizing noise sampled by Monte Carlo
    - writes the report to `data/out_test4a_qfor_h_depolarizing_mc.txt`
* test 4b — `test4b`
    - qfor H circuit with exact post-gate depolarizing noise using `Dist` and `collapse`
    - writes the report to `data/out_test4b_qfor_h_depolarizing_exact.txt`
* test 5a — `test5a`
    - quantamorphism with and without error correction to each qubit
    - should allows to see noise accumulation and that some error correction strategies are better than others
    - it is a **practical application of the combinator**
* test 5b — `test5b`
    - quantamorphism with and without error correction to each qubit
    - should allows to see noise accumulation and that some error correction strategies are better than others
    - it is a **practical application of the combinator**
* test 5c — `test5c`
    - compares cases where qubit 0, 1, or 2 is assigned 10% of the other gates' error probability; the quantamorphism target remains qubit 0
    - combines independent reset-to-`|0>` and phase-flip (`Z`) errors: probability `p/2` each
    - uses four qubits, starts from `|1>|1>|1>|0>`, and does not use correction ancillas
    - records fidelity of target qubit 0 for probabilities `0.0005`, `0.005`, and `0.05`

Run with no arguments (or an unrecognized test name) to print the current usage string,
since new tests/flags get added over time.

## Data & analysis

Each test writes its raw output as CSV under [data/](data/). Python scripts in the same
folder turn that raw output into the summary CSVs and comparison plots referenced above,
e.g. for test 5c:

```sh
cabal run tests -- test5c
python data/analyze_test5c.py
```

This reads `data/test5c_raw_Combined_1110_*.csv` and produces `data/test5c_clean.csv` and
`data/test5c_comparison.png`. Other tests follow the same raw-CSV-in, script-out pattern:
[analyze_test2.py](data/analyze_test2.py), [analyze_test2_grid.py](data/analyze_test2_grid.py),
and [analyze_test5_averages.py](data/analyze_test5_averages.py).

## License

MIT, see [LICENSE](LICENSE).
