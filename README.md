# quadprog
[![Build Status](https://travis-ci.org/rmcgibbo/quadprog.svg?branch=master)](https://travis-ci.org/rmcgibbo/quadprog)

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
- Build time
   - `numpy`, `cython`, C++ compiler.
  
