import numpy as np

def skew(vec: np.ndarray) -> np.ndarray:
    # Skew-symmetric matrix
    # Checks on vec.shape
    if vec.shape == (3,1) or vec.shape == (1,3):
        vec = vec.flatten()
        
    return np.array([[0, -vec[2], vec[1]],
                     [vec[2], 0, -vec[0]],
                     [-vec[1], vec[0], 0]])

def Rx(theta: float) -> np.ndarray:
    # Rotation matrix around the x-axis
    return np.array([[1, 0, 0],
                     [0, np.cos(theta), -np.sin(theta)],
                     [0, np.sin(theta), np.cos(theta)]])
def Ry(theta: float) -> np.ndarray:
    # Rotation matrix around the y-axis
    return np.array([[np.cos(theta), 0, np.sin(theta)],
                     [0, 1, 0],
                     [-np.sin(theta), 0, np.cos(theta)]])

def Rz(theta: float) -> np.ndarray:
    # Rotation matrix around the z-axis
    return np.array([[np.cos(theta), -np.sin(theta), 0],
                     [np.sin(theta), np.cos(theta), 0],
                     [0, 0, 1]])

