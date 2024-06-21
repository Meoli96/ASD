import numpy as np
from utils import skew

class Device:
    def __init__(self, mass, position: np.ndarray, axis: np.ndarray) -> None:
        self.const = 0 # Bool to check if the device adds time-dependency to system matrix
                       # If true, the device first and second order momenta are constant
        self.mass = mass
        self.b = position
        self.b_tilde = skew(self.b)
        self.axis = axis
        # from axis compute matrix of frame transformation
        self.id = 0 # Device ID --  represents the index of the device in the system matrix
                    #               must be > 6 (6 is the number of states of the craft)
        
        
    def set_id(self, id):
        self.id = id


class Wheel(Device):
    def __init__(self, mass, position: np.ndarray, axis: np.ndarray, radius, con_fun = lambda: 0) -> None:
        super().__init__(mass, position, axis)
        self.const = 1 # Wheel does not add time-dependency to system matrix
        self.locked = False
        self.radius = radius
        self.control = con_fun
        self.Is = 0.5*mass*radius**2
        self.It = 0.25*mass*radius**2 
        
        self.Sp = mass * self.rd
        self.Iw = np.array([[self.Is, 0, 0], [0, self.It, 0], [0, 0, self.It]])  + self.mass * (self.b@self.b.T - np.outer(self.b, self.b))
    
    def set_control(self, con_fun):
        self.control = con_fun
    
    def set_lock(self, lock_bool: bool):
        self.locked = lock_bool



class Damper(Device):
    def __init__(self, mass, position: np.ndarray, axis: np.ndarray, K, Cd) -> None:
        # Here position stands for the position of the damper at rest with respect to the craft frame
        super().__init__(mass, position, axis)
        self.K = K  
        self.Cd = Cd




     







class Craft:
    def __init__(self, mass, Inertia, com = np.zeros(3,1), devices = []) -> None:
        ### Spacecraft parameters
        self.Mp = mass
        self.com = com
        self.Ip = Inertia # Inertia matrix
        self.Sp = mass*self.com # Linear momentum

        ## System constant parameters
        self.Ms_const = self.Mp # Initial mass 
        self.Ss_const = self.Sp # Initial linear momentum
        self.Is_const = self.Ip # Initial inertia matrix

        ## System matrices (to initialize in the initialize function)
        self.Mt0 = None # Mass matrix
        self.Mt_dot = None # Time derivative of the mass matrix



        self.ndim = 6 # State space dimension
        self.ndim_aug = 6 # Augmented state space dimension for damper displacement
        self.dampers = []
        self.wheels = []
        if devices:
            self.add_devices(devices)  
        
        
        

        # Should also store analysis results (state, control, energy, energy rate, etc.)


    def add_devices(self, devices: list):
        '''Add a device to the craft. The device must be a Wheel or a Damper
        '''

        # Sort devices by type -- Damper first, Wheel second
        devices = sorted(devices, key=lambda device: isinstance(device, Damper), reverse=True)
        for device in devices:
            if type(device) == Damper or type(device) == Wheel:
                # Add device mass to the craft
                self.Ms_const += device.mass

                if type(device) == Damper:

                    device.set_id(self.ndim)
                    self.dampers.append(device)
                    self.ndim += 1
                    self.ndim_aug += 2 # Add two states for each damper



                elif type(device) == Wheel:

                    device.set_id(self.ndim)
                    self.wheels.append(device)
                    self.ndim += 1
                    self.ndim_aug += 1 # Add one state for each wheel

                    # Add constant Sw and Iw to the craft
                    self.Ss_const += device.Sp
                    self.Is_const += device.Iw


             
    
    def initialize(self, chi0: np.ndarray):
        '''Generate system matrix for the craft. Adding devices will change the system matrix
           so make sure to call this function after adding devices and before any simulation run
        '''
        self.Mt = np.zeros((self.ndim, self.ndim))
        self.Mt_dot = np.zeros((self.ndim, self.ndim))

        # Extract 


        # Mass matrix
        M = self.Ms_const
        # Static moment
        S = self.Ss_const

        # Inertia matrix
        I = self.Is_const # Inertia matrix -- constant term
        
        # Add dampers static momenta and inertia
        if self.dampers:
            for damper in self.dampers:
                rd = damper.b + damper.axis * chi0[damper.id]
                Sd_d = damper.mass * damper.rd
                Id_d = damper.mass * (damper.rd @ damper.rd.T*np.eye((3,3)) - np.outer(damper.rd, damper.rd))

                S += Sd_d
                I += Id_d
        
        self.Mt[0:3, 0:3] = M*np.eye(3) # Mass matrix
        self.Mt[3:6, 3:6] = I # Inertia matrix
        
        Ss_tilde = skew(S)
        self.Mt[0:3, 3:6] = - Ss_tilde
        self.Mt[3:6, 0:3] = Ss_tilde

        # Sort dampers and wheels for id
        devices = sorted(self.dampers + self.wheels, key=lambda device: device.id)
        for device in devices:
            if type(device) == Damper:
                self.Mt[device.id, device.id] = device.mass
                # vP
                self.Mt[0:3 , device.id:device.id+1] = device.mass*device.axis
                self.Mt[device.id : device.id+1, 0:3] = device.mass*device.axis.T
                
                # omega
                self.Mt[3:6, device.id : device.id+1] = device.mass*(device.b_tilde @ device.axis)
                self.Mt[device.id : device.id+1, 3:6] = -device.mass*(device.axis.T @ device.b_tilde)
            
            elif type(device) == Wheel:
                self.Mt[device.id, device.id] = device.Is
                # omega
                self.Mt[3:6 , device.id:device.id+1] = device.Is*device.axis
                self.Mt[device.id:device.id+1, 3:6] = device.Is*device.axis.T
         
        



    def update_Mt_Mt_dot(self, Mt, chi:np.ndarray, chi_dot:np.ndarray):
        '''Update the mass matrix and its derivative at time t
        '''
        S_temp = self.Ss_const
        I_temp = self.Is_const
        if self.dampers:
            Sdot_temp = np.zeros((3,1))
            Idot_temp = np.zeros((3,3))
        # Assuming chi[0] is of self.dampers[0]
        
        for i in range(len(self.dampers)):
            dmp = self.dampers[i]
            rd = dmp.b + dmp.axis * chi[i]
            rd_dot = dmp.axis * chi_dot[i]
            # Compute the static and dynamic momenta of i-th damper
            Sp_i = dmp.mass * rd
            Id_i = dmp.mass * (rd @ rd.T*np.eye((3,3)) - np.outer(rd, rd))
            # Add to the total
            S_temp += Sp_i
            I_temp += Id_i 
            # Compute the time derivative of the static and dynamic momenta of i-th damper
            Sp_i_dot = dmp.mass * rd_dot
            Id_i_dot = dmp.mass * ( 2*rd @ rd_dot.T - (np.outer(rd_dot, rd) + np.outer(rd, rd_dot)) )
            # Add to the total
            Sdot_temp += Sp_i_dot
            Idot_temp += Id_i_dot

        # Update the mass matrix
        # Static moment
        S_tilde = skew(S_temp)
        self.Mt[0:3, 3:6] = - S_tilde
        self.Mt[3:6, 0:3] = S_tilde

        # Inertia matrix
        self.Mt[3:6, 3:6] = I_temp

        # Time derivative of the mass matrix
        # Static moment
        Sdot_tilde = skew(Sdot_temp)
        self.Mt_dot[0:3, 3:6] = - Sdot_tilde
        self.Mt_dot[3:6, 0:3] = Sdot_tilde

        # Inertia matrix
        self.Mt_dot[3:6, 3:6] = Idot_temp

        return self.Mt, self.Mt_dot
        




    def f_ode(self, t, nu: np.ndarray):
        # unpack nu
        vP = nu[0:3]
        omega = nu[3:6]

        # unpack dampers and wheels
        chi_dot = nu[6:6+len(self.dampers)]
        w_s = nu[6+len(self.dampers): 6+len(self.dampers)+len(self.wheels)]
        
        # Augmented state (chi) always come for last
        chi = nu[-len(self.dampers):]

        Mt, Mtdot = self.update_Mt_Mt_dot(Mt, chi, chi_dot)
    def simulate(self, tspan, y0, options):
        pass
    def plot(self):
        pass


