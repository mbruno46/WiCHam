# WiCHam
========

A Mathematica library to compute the running of the Wilson coefficients for the Weak Hamiltonian

## Setup 
========

 * Download the git repository and place the directory WiCHam in your home

 * Add the following line to your Mathematica notebook file

```
<<"/path/to/WiCHam/WiCHam.m"
```

## Documentation
================

The functions contained in this library are documented in the file **docs/manual.pdf**

## Usage
========

To compute the Wilson Coefficients at some scale mu to LO type in your notebook

```
a1 = value_of_alpha_strong_at_weak_scale;
loop_alpha = 1;
z = ComputeZ[mu,a1,loop_alpha];
ReduceOrder[z,"LO"]

y = ComputeY[z,mu,a1,loop_alpha];
ReduceOrder[y,"LO"]

```

## Tests
========

In the file **test.nb** several results presented in the Reference below are correctly 
reproduced and checked. In particular LO and NLO evaluations of the Wilson Coefficients
in Tables X, XII, XIII, XVIII, XIX and XX.

## References
=============

[Weak decays beyond leading logarithms](http://inspirehep.net/record/403867)
