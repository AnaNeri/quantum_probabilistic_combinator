# Quantum Probabilistic Combinator
Tests of quantum probabilistic combinator in haskell.

* matrices.hs - has quantum gates 
    * Pauli matrices: x, y, z
    * hadamard gate 
    * identity gate
    * cx gate
    * cu3 gate
* probabilistic_combinator.hs
* tests.hs
    * test 1    
        - check noise accumulation in circuit X vs circuit HZH 
        - shows that more depth usually results in more noise         
    * test 2 
        - check noise in circuit X with and without error correction bit flip error (implementation uses DOI: [10.13140/RG.2.2.18542.77129](https://www.researchgate.net/publication/334634646_The_first_three-qubit_and_six-qubit_full_quantum_multiple_error-correcting_codes_with_low_quantum_costs?channel=doi&linkId=5d36ffe2a6fdcc370a57ac6a&showFulltext=true))
        - Allow to see that error correction reduces error
    * test 3 
        - quantamorphism with and without error correction to each qubit
        - should allows to see noise accumulation and that some error correction strategies are better than others
        - it is a **practical application of the combinator**
