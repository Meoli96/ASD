import numpy as np

def skew(vec: np.ndarray) -> np.ndarray:
    return np.array([[0, -vec[2], vec[1]],
                     [vec[2], 0, -vec[0]],
                     [-vec[1], vec[0], 0]])

