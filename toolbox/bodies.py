import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from .utils import skew

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
    def __init__(self, mass, position: np.ndarray, axis: np.ndarray, radius, con_fun = lambda x: 0) -> None:
        super().__init__(mass, position, axis)
        self.locked = False
        self.radius = radius
        self.control = con_fun
        self.Is = 0.5*mass*radius**2
        self.It = 0.25*mass*radius**2 
        
        self.Sp = mass * self.b
        self.Iw = np.array([[self.Is, 0, 0], [0, self.It, 0], [0, 0, self.It]])  + self.mass * (self.b.T@self.b*np.eye(3,3) - np.outer(self.b, self.b))
    
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
    def __init__(self, mass, Inertia, CoM = np.zeros([3,1]), devices = []) -> None:
        ### Spacecraft parameters
        self.Mp = mass
        self.CoM = CoM
        self.Ip = Inertia # Inertia matrix
        self.Sp = mass*self.CoM # Linear momentum

        ## System constant parameters
        self.Ms_const = self.Mp # Initial mass 
        self.Ss_const = self.Sp # Initial linear momentum
        self.Is_const = self.Ip # Initial inertia matrix

        ## System matrices and initial state (to initialize in the initialize function)
        self.Mt = None # Mass matrix
        self.Mt_dot = None # Time derivative of the mass matrix



        self.ndim = 6 # State space dimension
        self.ndim_aug = 6 # Augmented state space dimension for damper displacement
        self.dampers = []
        self.wheels = []
        self.results = None # Store simulation results, np.array
        self.T  = [] # Energy
        if devices:
            self.add_devices(devices)  
        
        
        

        # Should also store analysis results (state, control, energy, energy rate, etc.)

    def lock_wheels(self, lock_bool: bool = True, wheel_id: list = []):
        # Lock or unlock wheels, locking all wheels if wheel_id is empty
        if not wheel_id:
            for wheel in self.wheels:
                wheel.set_lock(lock_bool)
        else:
            for id_w in wheel_id:
                self.wheels[id_w].set_lock(lock_bool)

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
            else:
                raise ValueError("Device must be a Wheel or a Damper")
        


             
    
    def initialize(self, chi0: np.ndarray):
        '''Generate system matrix for the craft. Adding devices will change the system matrix
           so make sure to call this function after adding devices and before any simulation run
        '''
        
        # Assert initial state dimension is equal to the augmented state dimension
        assert len(self.dampers) == chi0.shape[0], "Initial state dimension must be equal to the number of dampers"
        # Initialize system matrix and its derivativee  
        self.Mt = np.zeros((self.ndim, self.ndim))
        self.Mt_dot = np.zeros((self.ndim, self.ndim))

        

        M, S, I = 0.0, 0.0, 0.0
        # Mass matrix
        M += self.Ms_const
        # Static moment
        S += self.Ss_const

        # Inertia matrix
        I += self.Is_const # Inertia matrix -- constant term
        
        # Add dampers static momenta and inertia
        if self.dampers:
            for i in range(len(self.dampers)):
                damper = self.dampers[i]
                rd = damper.b + damper.axis * chi0[i]
                Sd_d = damper.mass * rd
                Id_d = damper.mass * (rd.T @ rd*np.eye(3,3) - np.outer(rd, rd))

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
        self.T = [] # Energy
        



    def update_Mt_Mt_dot(self, chi:np.ndarray, chi_dot:np.ndarray):
        '''Update the mass matrix and its derivative at time t
        '''
        S_temp = 0.0
        I_temp = 0.0
        S_temp += self.Ss_const
        I_temp += self.Is_const
        
        Sdot_temp = np.zeros((3,1))
        Idot_temp = np.zeros((3,3))
        # Assuming chi[0] is of self.dampers[0]
        
        for i in range(len(self.dampers)):
            dmp = self.dampers[i]
            rd = dmp.b + dmp.axis * chi[i]
            rd_dot = dmp.axis * chi_dot[i]
            # Compute the static and dynamic momenta of i-th damper
            Sp_i = dmp.mass * rd
            Id_i = dmp.mass * (rd.T @ rd*np.eye(3,3) - np.outer(rd, rd))

            # Add to the total
            S_temp += Sp_i
            I_temp += Id_i 
            
            # Compute the time derivative of the static and dynamic momenta of i-th damper
            Sp_i_dot = dmp.mass * rd_dot
            Id_i_dot = dmp.mass * ( 2*rd.T @ rd_dot*np.eye(3,3) - (np.outer(rd_dot, rd) + np.outer(rd, rd_dot)) )
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

        return 
        




    def f_ode(self, t, nu: np.ndarray, Fext: np.ndarray, Mext: np.ndarray):
        # unpack and unflatten nu
        nu = nu.reshape((self.ndim_aug, 1))
        vP = nu[0:3]
        omega = nu[3:6]

        # unpack dampers and wheels
        chi_dot = nu[6:6+len(self.dampers)]
        w_s = nu[6+len(self.dampers): 6+len(self.dampers)+len(self.wheels)]
        
        # Augmented state (chi) always come for last, we can also strip it from nu
        chi = nu[-len(self.dampers):] 
        _nu = nu[0:self.ndim] # exclude chi


        # Update mass matrix and its derivative
        self.update_Mt_Mt_dot(chi, chi_dot)
        ## Build nu_dot = inv(Mt-B)@((A@Mt - Mt_dot)@om - c)
        A = np.zeros((self.ndim, self.ndim))
        B = np.zeros((self.ndim, self.ndim))
        C = np.zeros((self.ndim, 1))
        # A
        A[0:3, 0:3] = -skew(omega)
        A[3:6, 0:3] = -skew(vP)
        A[3:6, 3:6] = -skew(omega)

        # B

        for wheel in self.wheels:
            if wheel.locked:
                B[wheel.id, 3:6] = wheel.Is * wheel.axis.T
                C[wheel.id, 0] = 0
            else:
                C[wheel.id, 0] = wheel.control(t, _nu)
        
        # C 
        C[0:3, 0:1] = Fext
        C[3:6, 0:1] = Mext
        for i in range(len(self.dampers)):
            damper = self.dampers[i]
            rd = damper.b + damper.axis * chi[i]
            C[damper.id, 0] = damper.mass*omega.T@skew(damper.axis)@(vP - skew(rd)@omega) - damper.K*chi[i] - damper.Cd*chi_dot[i]

        # Solve for nu_dot	
        nu_dot = np.linalg.inv(self.Mt - B) @ ((A @ self.Mt - self.Mt_dot) @ _nu + C)
        # Augment nu_dot
        nu_dot = np.concatenate((nu_dot, chi_dot), axis=0)

        # Compute mechanical energy
        V = 0.0
        for i in range(len(self.dampers)):
            damper = self.dampers[i]
            V += 0.5*damper.K*chi[i]**2

        self.T.append((0.5 * _nu.T @ self.Mt @ _nu)[0][0] + V)
        
        return nu_dot.flatten()




    def simulate(self, nu0: np.ndarray, tspan, Fext = np.zeros((3,1)), Mext = np.zeros((3,1)),  options = None):
        # This function starts the simulation of the craft for the given time span
        # Results should be accumulated in some tyoe of data structure for analysis
        
        # split the time span
        t0, tf = tspan
        # Initialize the results array
        self.results = np.zeros((tf-t0, self.ndim_aug))
        

        # Now call solve_ivp with initial conditions nu0
        res = solve_ivp(self.f_ode, (t0, tf), nu0, t_eval=np.linspace(t0, tf, tf-t0 ), args=(Fext, Mext), method='RK45')
        self.results = res.y.T
        
        


    def plot(self):
        # Plot the results of the simulation
        if self.results is None:
            raise ValueError("No simulation results to plot")
        else:
            ### unpack the results
            # body angular velocities
            om  = self.results[:, 3:6]
            # Damper displacemnts and velocities
            chi = self.results[:, -len(self.dampers):]
            chi_dot = self.results[:, 6:6+len(self.dampers)]
            # Wheel angular velocities
            w_s = self.results[:, 6+len(self.dampers):6+len(self.dampers)+len(self.wheels)]
            # Energy
            T = self.T


            ### Plot angular velocities
            # Convert om to degrees
            om = np.degrees(om)
            plt.figure()
            plt.grid()
            plt.title("Angular velocities")
            plt.plot(om)
            plt.legend(["$\omega_x$", "$\omega_y$", "$\omega_z$"])         
            plt.xlabel("Time")
            plt.ylabel("Angular velocity(deg/s)")
            plt.show()

            ### Plot damper velocities and displacements
            # Displacements
            plt.figure()
            plt.grid()
            plt.title("Damper displacements")
            plt.plot(chi)
            # Display legend based on number of dampers
            leg = [f"$\chi_{i}$" for i in range(len(self.dampers))]
            plt.legend(leg)
            plt.xlabel("Time")
            plt.ylabel("Displacement")
            plt.show()

            # Velocities
            plt.figure()
            plt.grid()
            plt.title("Damper velocities")
            plt.plot(chi_dot)
            # Display legend based on number of dampers
            leg = [f"$\dot\chi_{i}$" for i in range(len(self.dampers))]
            plt.legend(leg)
            plt.xlabel("Time")
            plt.ylabel("Velocity")
            plt.show()

            ### Plot wheel angular velocities
            plt.figure()
            plt.grid()
            plt.title("Wheel angular velocities")
            plt.plot(w_s)
            # Display legend based on number of wheels
            leg = [f"$\omega_{i}$" for i in range(len(self.wheels))]
            plt.legend(leg)
            plt.xlabel("Time")
            plt.ylabel("Angular velocity")
            plt.show()

            ### Plot energy
            plt.figure()
            plt.grid()
            plt.title("Energy")
            plt.plot(T)
            plt.xlabel("Time")
            plt.ylabel("Energy")
            plt.show()



