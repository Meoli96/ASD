import numpy as np
class Craft:
    def __init__(self, mass, Inertia, com = np.zeros(3,1), devices = []) -> None:
        # System matrix
        self.mass = mass
        self.com = com # Center of mass displacement from origin
        
        self.I = Inertia # Inertia matrix
        self.Sp = mass*self.com # Linear momentum

        self.ndim = 6 # vP and omegaP
        if devices:
            self.devices = devices
            for device in devices:
                self.mass += device.mass
                self.I += device.Inertia
        self.devices = []

        # Should also store analysis results (state, control, energy, energy rate, etc.)


    def add_device(self, device):   
        if type(device) == Wheel or type(device) == Damper:
            self.devices.append(device)
            self.ndim += 1
        else:
            # type not recognized
            raise ValueError('Device type not recognized')

    
    def gen_sysmat(self):
        '''Generate system matrix for the craft. Adding devices will change the system matrix
           so make sure to call this function after adding devices and before any simulation run
        '''

    
    def update_sysmat(self, t, y):
        pass
    def simulate(self, tspan, y0, options):
        pass
    def plot(self):
        pass

class Device:
    def __init__(self, mass, position, axis) -> None:
        self.const = 0 # Bool to check if the device adds time-dependency to system matrix
                       # If true, the device first and second order momenta are constant
        self.mass = mass
        self.position = position
        self.axis = axis
        # from axis compute matrix of frame transformation
        self.T = np.outer(axis, np.eye(3)) ### ??? Check this

class Wheel(Device):
    def __init__(self, mass, position, axis, radius, con_fun = lambda: 0) -> None:
        super().__init__(mass, position, axis)
        self.const = 1 # Wheel does not add time-dependency to system matrix
        self.radius = radius
        self.control = con_fun
        Is = 0.5*mass*radius**2
        It = 0.25*mass*radius**2 
        self.inertia = np.array([[Is, 0, 0], [0, It, 0], [0, 0, It]]) # Calculate inertia in wheel frame
    
    def set_control(self, con_fun):
        self.control = con_fun

class Damper(Device):
    def __init__(self, mass, position, inertia, axis, K, Cd) -> None:
        super().__init__(mass, position, axis)
        self.K = K
        self.Cd = Cd
        

