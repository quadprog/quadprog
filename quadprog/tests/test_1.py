import numpy as np
from quadprog import solve_qp

def test_1():
    G = np.eye(3, 3)
    a = np.zeros(3)
    print(solve_qp(G, a))



test_1()
