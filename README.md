# quadprog: Quadratic Programming Solver (Python)

[![.github/workflows/build-and-test.yaml](https://github.com/quadprog/quadprog/actions/workflows/build-and-test.yaml/badge.svg?branch=master)](https://github.com/quadprog/quadprog/actions/workflows/build-and-test.yaml)

```
Solve a strictly convex quadratic program

Minimize     1/2 x^T G x - a^T x
Subject to   C.T x >= b

This routine uses the the Goldfarb/Idnani dual algorithm [1].

References
---------
... [1] D. Goldfarb and A. Idnani (1983). A numerically stable dual
    method for solving strictly convex quadratic programs.
    Mathematical Programming, 27, 1-33.
```

### Installation
`pip install quadprog`

### Dependencies
- Runtime
   - `numpy`
- Installation
   - `numpy`
   - C compiler if installing from sdist
   - `Cython` if building from source

### Developing

See [docs/DEVELOP.md](docs/DEVELOP.md).