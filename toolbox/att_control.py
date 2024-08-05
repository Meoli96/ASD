import numpy as np
class CelestialBody:
    def __init__(self, mu, radius):
        self.mu = mu # km^3/s^2
        self.radius = radius # km

Earth = CelestialBody(398600, 6378)
Mars = CelestialBody(42828, 3396)
Sun = CelestialBody(132712440041.94, 696340)

class KeplOrbit:
    # Set gravitational parameter
    def __init__(self, a, e, i, RAAN, omega, truan0, central_body: CelestialBody = Earth):
        
        # Sanity checks on orbital elements
        assert a > 0, "Semi-major axis must be positive"
        assert e >= 0, "Eccentricity must be non-negative"
        assert i >= 0 and i <= np.pi, "Inclination must be between 0 and pi"
        assert RAAN >= 0 and RAAN <= 2*np.pi, "RAAN must be between 0 and 2*pi"
        assert omega >= 0 and omega <= 2*np.pi, "Argument of perigee must be between 0 and 2*pi"
        assert truan0 >= 0 and truan0 <= 2*np.pi, "True anomaly must be between 0 and 2*pi"

        
        self.central_body = central_body
        
        self.a = a
        self.e = e
        self.i = i
        self.RAAN = RAAN
        self.omega = omega
        self.truan0 = truan0

    def step(self, t):
        # Calculate the true anomaly at time t

        # If circular orbit, it is a special case
        if self.e == 0:
            return np.sqrt(self.central_body.mu/self.a**3)*t + self.truan0
        elif self.e < 1:
            # Valid for eccentric orbits
            # Calculate the eccentric anomaly
            M = np.sqrt(self.central_body.mu/self.a**3)*t
            E = M
            for _ in range(10):
                # Newton-Raphson method
                E = M + self.e*np.sin(E)
            # Calculate the true anomaly
            nu = 2*np.arctan(np.sqrt((1+self.e)/(1-self.e))*np.tan(E/2))
            return nu
            



class Craft:
    def __init__(self, Jb):
        self.inertia = Jb

    def control(self, desired_frame, desired_angular_velocity, external_torque):
        # Desired frame can be a rotation matrix or a quaternion, constant
        # or a function of some parameters
        # Desired angular velocity is a 3x1 vector, constant or a function
        # of some parameters

        target_frame = desired_frame()


        pass 