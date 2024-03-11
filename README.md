# WiCHam

A Mathematica library to compute the running of the Wilson coefficients for the Weak Hamiltonian

- [**Documentation**](./docs/manual.pdf)
- **Examples:** [tests](./tests)
- [**Source code**](.//WiCHam)
- **Bug reports:** https://github.com/mbruno46/WiCHam/issues

## Citation

This library implements the results derived in [Weak decays beyond leading logarithms](http://inspirehep.net/record/403867).

If you use this library in your publication cite this article.

## Authors

Mattia Bruno, Copyright (c) 2016-2024

## Setup 

 * Download the git repository and place the directory WiCHam in your home

 * Add the following line to your Mathematica notebook file

```
Get["/path/to/WiCHam/WiCHam.m"]
```

## Example

To compute the Wilson Coefficients at some scale mu to LO type in your notebook

```
a1 = value_of_alpha_strong_at_weak_scale;
loop_alpha = 1;
z = ComputeZ[mu,a1,loop_alpha];
ReduceOrder[z,"LO"]

y = ComputeY[z,mu,a1,loop_alpha];
ReduceOrder[y,"LO"]
```

The notebook file [work.nb](./work.nb) is ready for a calculation, just use it!

## Tests

Several results presented in the main Reference above are correctly 
reproduced and checked. In particular LO and NLO evaluations of the Wilson Coefficients
in Tables X, XVIII and XX.

