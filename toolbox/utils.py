import numpy as np

# class Quat:
#     # Class to handle quaternions
#     def __init__(self, q: np.ndarray = np.array([1, 0, 0, 0])):
#         # Sanity check on quaternion
#         assert q.shape == (4,), "Quaternion must be a 4-element array"
#         if not (np.isclose(np.linalg.norm(q), 1)):
#             # Normalize the quaternion
#             q = q/np.linalg.norm(q)
#         self.q = q
        

#     def __add__(self, other):
#         return Quat(self.q + other.q)
#     def __sub__(self, other):
#         return Quat(self.q - other.q)
#     def __eq__(self, other):
#         return np.allclose(self.q, other.q) or np.allclose(self.q, -other.q)
#     # Da fixare moltiplicazione
#     def __mul__(self, other):
#         return Quat(np.array([self.q[0]*other.q[0] + np.dot(self.q[1:], other.q[1:]),
#                               self.q[0]*other.q[1:] + other.q[0]*self.q[1:] + np.cross(self.q[1:], other.q[1:])]))
#     # [ operator
#     def __getitem__(self, key):
#         return self.q[key]
#     def __str__(self):
#         return str(self.q)
    
#     def conj(self):
#         return Quat(np.array([self.q[0], -self.q[1], -self.q[2], -self.q[3]]))
    
#     def dot(self, other):
#         return np.dot(self.q[1:], other.q[1:])
#     def cross(self, other):
#         return Quat(np.array([]))
#     def toR(self):
#         return np.array([[1-2*(self.q[2]**2+self.q[3]**2), 2*(self.q[1]*self.q[2]-self.q[0]*self.q[3]), 2*(self.q[1]*self.q[3]+self.q[0]*self.q[2])],
#                          [2*(self.q[1]*self.q[2]+self.q[0]*self.q[3]), 1-2*(self.q[1]**2+self.q[3]**2), 2*(self.q[2]*self.q[3]-self.q[0]*self.q[1])],
#                          [2*(self.q[1]*self.q[3]-self.q[0]*self.q[2]), 2*(self.q[2]*self.q[3]+self.q[0]*self.q[1]), 1-2*(self.q[1]**2+self.q[2]**2)]])
#     def fromR(self, R: np.ndarray):
#         # Convert a rotation matrix to a quaternion
#         assert R.shape == (3, 3), "Rotation matrix must be 3x3"
#         q = np.zeros(4)
#         T = np.trace(R)
#         if T > 0:
#             s = 0.5/np.sqrt(T)
#             w = 0.25/s
#             x = (R[2,1] - R[1,2])*s
#             y = (R[0,2] - R[2,0])*s
#             z = (R[1,0] - R[0,1])*s
#         else:
#             if R[0,0] > R[1,1] and R[0,0] > R[2,2]:
#                 s = 2*np.sqrt(1 + R[0,0] - R[1,1] - R[2,2])
#                 w = (R[2,1] - R[1,2])/s
#                 x = 0.25*s
#                 y = (R[0,1] + R[1,0])/s
#                 z = (R[0,2] + R[2,0])/s
#             elif R[1,1] > R[2,2]:
#                 s = 2*np.sqrt(1 + R[1,1] - R[0,0] - R[2,2])
#                 w = (R[0,2] - R[2,0])/s
#                 x = (R[0,1] + R[1,0])/s
#                 y = 0.25*s
#                 z = (R[1,2] + R[2,1])/s
#             else:
#                 s = 2*np.sqrt(1 + R[2,2] - R[0,0] - R[1,1])
#                 w = (R[1,0] - R[0,1])/s
#                 x = (R[0,2] + R[2,0])/s
#                 y = (R[1,2] + R[2,1])/s
#                 z = 0.25*s
#         self.q = np.array([w, x, y, z])
        
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

# def rot2quat(rotmat: np.ndarray) -> Quat:
#     # Converts a rotation matrix to a quaternion
#     assert rotmat.shape == (3, 3), "Rotation matrix must be 3x3"
#     q = np.zeros(4)
#     if rotmat[2,2] < 0:
#         if rotmat[0,0] > rotmat[1,1]:
#             t = 1 + rotmat[0,0] - rotmat[1,1] - rotmat[2,2]
#             q = np.array([t, rotmat[0,1] + rotmat[1,0], rotmat[2,0] + rotmat[0,2], rotmat[1,2] - rotmat[2,1]])
#         else:
#             t = 1 - rotmat[0,0] + rotmat[1,1] - rotmat[2,2]
#             q = np.array([rotmat[0,1] + rotmat[1,0], t, rotmat[1,2] + rotmat[2,1], rotmat[2,0] - rotmat[0,2]])
#     else:
#         if rotmat[0,0] < -rotmat[1,1]:
#             t = 1 - rotmat[0,0] - rotmat[1,1] + rotmat[2,2]
#             q = np.array([rotmat[2,0] + rotmat[0,2], rotmat[1,2] + rotmat[2,1], t, rotmat[0,1] - rotmat[1,0]])
#         else:
#             t = 1 + rotmat[0,0] + rotmat[1,1] + rotmat[2,2]
#             q = np.array([rotmat[1,2] - rotmat[2,1], rotmat[2,0] - rotmat[0,2], rotmat[0,1] - rotmat[1,0], t])
#     q = q / (2*np.sqrt(t))
#     return Quat(q)

