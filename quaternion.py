# Class to handle quaternions
import numpy as np

class Quaternion:
    def __init__(self, q: np.ndarray = np.array([1, 0, 0, 0])):
        # Sanity check on quaternion
        assert q.shape == (4,) or q.shape == (1,4) or q.shape == (4,1), "Quaternion must be a 4-element array"
        if not (np.isclose(np.linalg.norm(q), 1)):
            # Normalize the quaternion
            q = q/np.linalg.norm(q)
        self._q0 = q[0]
        self._q = q[1:4] # Reshape to a 3-element column vector

    def q0(self):
        return self._q0

    def q(self):
        # Return the quaternion as a 3-element column vector
        return self._q.reshape((3, 1))

    def __getitem__(self, key):
        if isinstance(key, slice):
            indices = range(*key.indices(4))
            if 0 in indices:
                return [self._q0] + [self._q[i-1] for i in indices if i != 0]
            return [self._q[i-1] for i in indices]
        elif key == 0:
            return self._q0
        return self._q[key-1]

    def __add__(self, other):
        return Quaternion(np.array([self._q0 + other._q0, self._q + other._q]))

    def __sub__(self, other):
        return Quaternion(np.array([self._q0 - other._q0, self._q - other._q]))

    def __eq__(self, other):
        return (np.allclose(self._q0, other._q0) and np.allclose(self._q, other._q)) \
            or (np.allclose(self._q0, -other._q0) and np.allclose(self._q, -other._q))

    def __matmul__(self, other):
        return Quaternion(np.hstack([self._q0*other._q0 - np.dot(self._q, other._q),
                                    self._q0*other._q + other._q0*self._q + np.cross(self._q, other._q)]))

    def conj(self):
        # Conjugate of a quaternion
        return Quaternion(np.hstack([self._q0, -self._q]))

    def toAxisAngle(self):
        # Convert quaternion to axis-angle representation
        angle = 2*np.arccos(self._q0)
        if np.isclose(angle, 0):
            return np.array([0, 0, 0]), 0
        axis = self._q/np.sin(angle/2)
        return axis, angle

    def fromAxisAngle(self, axis, angle):
        # Convert axis-angle representation to quaternion
        self._q0 = np.cos(angle/2)
        self._q = axis*np.sin(angle/2)
        return self

    def toR(self):
        R = np.zeros((3, 3))
        R[0, 0] = 1 - 2*(self._q[1]**2 + self._q[2]**2)
        R[0, 1] = 2*(self._q[0]*self._q[1] - self._q[2]*self._q0)
        R[0, 2] = 2*(self._q[0]*self._q[2] + self._q[1]*self._q0)
        R[1, 0] = 2*(self._q[0]*self._q[1] + self._q[2]*self._q0)
        R[1, 1] = 1 - 2*(self._q[0]**2 + self._q[2]**2)
        R[1, 2] = 2*(self._q[1]*self._q[2] - self._q[0]*self._q0)
        R[2, 0] = 2*(self._q[0]*self._q[2] - self._q[1]*self._q0)
        R[2, 1] = 2*(self._q[1]*self._q[2] + self._q[0]*self._q0)
        R[2, 2] = 1 - 2*(self._q[0]**2 + self._q[1]**2)
        return R

    def fromR(self, R: np.ndarray):
        # Convert a rotation matrix to a quaternion
        assert R.shape == (3, 3), "Rotation matrix must be 3x3"
        q = np.zeros(4)
        T = np.trace(R)
        if T > 0:
            s = 0.5/np.sqrt(T)
            w = 0.25/s
            x = (R[2,1] - R[1,2])*s
            y = (R[0,2] - R[2,0])*s
            z = (R[1,0] - R[0,1])*s
        else:
            if R[0,0] > R[1,1] and R[0,0] > R[2,2]:
                s = 2*np.sqrt(1 + R[0,0] - R[1,1] - R[2,2])
                w = (R[2,1] - R[1,2])/s
                x = 0.25*s
                y = (R[0,1] + R[1,0])/s
                z = (R[0,2] + R[2,0])/s
            elif R[1,1] > R[2,2]:
                s = 2*np.sqrt(1 + R[1,1] - R[0,0] - R[2,2])
                w = (R[0,2] - R[2,0])/s
                x = (R[0,1] + R[1,0])/s
                y = 0.25*s
                z = (R[1,2] + R[2,1])/s
            else:
                s = 2*np.sqrt(1 + R[2,2] - R[0,0] - R[1,1])
                w = (R[1,0] - R[0,1])/s
                x = (R[0,2] + R[2,0])/s
                y = (R[1,2] + R[2,1])/s
                z = 0.25*s
        self._q0 = w
        self._q = np.array([x, y, z])
        return self

    def toArray(self):
        return np.array([self._q0] + list(self._q))


if __name__ == '__main__':
    # Test the quaternion class
    q = Quaternion(np.array([0.1, 0.5, 0.5, 0.5]))
    print(q[0], q[1], q[2], q[3])
    print(q[0:3])
    print(q.toArray())
