import numpy as np

def skew(vec: np.ndarray) -> np.ndarray:
    # Skew-symmetric matrix
    # Checks on vec.shape
    if vec.shape == (3,1) or vec.shape == (1,3):
        vec = vec.flatten()
        
    return np.array([[0, -vec[2], vec[1]],
                     [vec[2], 0, -vec[0]],
                     [-vec[1], vec[0], 0]])
