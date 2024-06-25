import numpy as np

class Quat:
    # Class to handle quaternions
    def __init__(self, q: np.ndarray):
        assert q.shape == (4,), "Quaternion must be a 4x1 vector"
        self.q = q
    def __add__(self, other):
        return Quat(self.q + other.q)
    def __sub__(self, other):
        return Quat(self.q - other.q)
    def __mul__(self, other):
        return Quat(np.array([self.q[0]*other.q[0] - np.dot(self.q[1:], other.q[1:]),
                              self.q[0]*other.q[1:] + other.q[0]*self.q[1:] + np.cross(self.q[1:], other.q[1:])]))
    def __str__(self):
        return str(self.q)
    def conj(self):
        return Quat(np.array([self.q[0], -self.q[1], -self.q[2], -self.q[3]]))
    def norm(self):
        return np.linalg.norm(self.q)
    

def skew(vec: np.ndarray) -> np.ndarray:
    # Convert vec to one dimenion
    vec = vec.reshape(3)
    # return the skew-symmetric matrix of vec
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

def Rx_dot (theta: float) -> np.ndarray:
    # Derivative of the rotation matrix around the x-axis
    return np.array([[0, 0, 0],
                     [0, -np.sin(theta), -np.cos(theta)],
                     [0, np.cos(theta), -np.sin(theta)]])

def Ry_dot (theta: float) -> np.ndarray:
    # Derivative of the rotation matrix around the y-axis
    return np.array([[-np.sin(theta), 0, np.cos(theta)],
                     [0, 0, 0],
                     [-np.cos(theta), 0, -np.sin(theta)]])

def Rz_dot (theta: float) -> np.ndarray:    
    # Derivative of the rotation matrix around the z-axis
    return np.array([[-np.sin(theta), -np.cos(theta), 0],
                     [np.cos(theta), -np.sin(theta), 0],
                     [0, 0, 0]])

def rot2quat(rotmat: np.ndarray) -> Quat:
    # Converts a rotation matrix to a quaternion
    assert rotmat.shape == (3, 3), "Rotation matrix must be 3x3"
    q = np.zeros(4)
    if rotmat[2,2] < 0:
        if rotmat[0,0] > rotmat[1,1]:
            t = 1 + rotmat[0,0] - rotmat[1,1] - rotmat[2,2]
            q = np.array([t, rotmat[0,1] + rotmat[1,0], rotmat[2,0] + rotmat[0,2], rotmat[1,2] - rotmat[2,1]])
        else:
            t = 1 - rotmat[0,0] + rotmat[1,1] - rotmat[2,2]
            q = np.array([rotmat[0,1] + rotmat[1,0], t, rotmat[1,2] + rotmat[2,1], rotmat[2,0] - rotmat[0,2]])
    else:
        if rotmat[0,0] < -rotmat[1,1]:
            t = 1 - rotmat[0,0] - rotmat[1,1] + rotmat[2,2]
            q = np.array([rotmat[2,0] + rotmat[0,2], rotmat[1,2] + rotmat[2,1], t, rotmat[0,1] - rotmat[1,0]])
        else:
            t = 1 + rotmat[0,0] + rotmat[1,1] + rotmat[2,2]
            q = np.array([rotmat[1,2] - rotmat[2,1], rotmat[2,0] - rotmat[0,2], rotmat[0,1] - rotmat[1,0], t])
    q = q / (2*np.sqrt(t))
    return Quat(q)

