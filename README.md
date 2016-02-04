# WiCHam

A Mathematica library to compute the running of the Wilson coefficients for the
Weak Hamiltonian

## Setup 

 * Download the git repository and place the directory WiCHam in your home

 * Add the following line to your Mathematica notebook file

```
<<"~/path/to/WiCHam/WiCHam.m"
```

## Usage

To compute the Wilson Coefficients at some scale mu to LO type in your notebook

```
z = ComputeZ[mu,initial_value_of_alpha_s,loop_alpha];
ReduceOrder[z,LO]

y = ComputeY[z,mu,initial_value_of_alpha_s,loop_alpha];
ReduceOrder[y,LO]

```

## References

[Weak decays beyond leading logarithms](http://inspirehep.net/record/403867)
