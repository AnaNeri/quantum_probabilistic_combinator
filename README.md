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
        - should show that more depth usually results in more noise         
    * test 2 
        - check noise in circuit X with and without error correction bit flip error
        - should allow to see that error correction reduces error
    * test 3 
        - quantamorphism with and without error correction to each qubit
        - should allows to see noise accumulation and that some error correction strategies are better than others
        - it is a **practical application of the combinator**