import numpy as np

# Define body frame unit vectors
b1 = np.array([[1], [0], [0]])
b2 = np.array([[0], [1], [0]])
b3 = np.array([[0], [0], [1]])

b = 0.5 # absolute distance from center of mass to rotor and dampers

Fext = np.array([[0], [0], [0]]) # external force
Mext = np.array([[0], [0], [0]]) # external moment

# Define the spacecraft parameters
Msp = 200 # kg
Ssp = np.array([[0], [0], [0]]) # kg*m
Jsp = np.array([[350, 0, 0], [0, 200, 0], [0, 0, 500]]) # kg*m^2



# Define wheel parameters
Mw = 20 # kg
Rw = 1 # m
b_w = b*b1 # m, position of wheel from center of mass
a_w = b1 # wheel axis direction

Is = 0.5*Mw*Rw**2 # kg*m^2, wheel inertia
It = 0.25*Mw*Rw**2 # kg*m^2, wheel inertia

Sw = Mw*b_w # kg*m, wheel mass moment 
Jww_b = np.array([[Is, 0, 0], [0, It, 0], [0, 0, It]]) + Mw*(np.dot(b_w.T, b_w)*np.eye(3) - np.outer(b_w, b_w)) # kg*m^2, wheel inertia in body frame

wheel_locked = False
ga_control = lambda t, x: np.array([[0], [0], [0]]) # control function

om_s01 = 2*np.pi # rad/s
om_s02 = 3*np.pi # rad/s

Kw = 0.1 # Nm/rad/s

# Define damper parameters
Md = 5 # kg
Kd = 3 # Kg/s^2
Cd = 5 # Kg/s

# Damper positions and directions
b_d1 = b*b3
a_d1 = b3
S_d1 = Md*b_d1

b_d2 = -b*b2
a_d2 = -b2
S_d2 = Md*b_d2

b_d3 = -b*b3
a_d3 = -b3
S_d3 = Md*b_d3

b_d4 = b*b2
a_d4 = b2
S_d4 = Md*b_d4

b_ds = [b_d1, b_d2, b_d3, b_d4]
a_ds = [a_d1, a_d2, a_d3, a_d4]


def skew(vec: np.ndarray) -> np.ndarray:
    # Skew-symmetric matrix
    # Checks on vec.shape
    if vec.shape == (3,1) or vec.shape == (1,3):
        vec = vec.flatten()
        
    return np.array([[0, -vec[2], vec[1]],
                     [vec[2], 0, -vec[0]],
                     [-vec[1], vec[0], 0]])
### System matrix definition
def initialize_Mt():
    # This function resets the global Mt matrix to the initial state
    # This is necessary to avoid the accumulation of values in the matrix due to it being a global variable
    
    global Mt, Mt_dot
    Mt = np.zeros((11, 11))
    Mt_dot = np.zeros((11, 11))

    M = Msp + Mw + 4*Md # Total mass of the system
    Sp0 = skew(Ssp + Sw) # Total static moment of the system
    # Add damper static moments
    for i in range(4):
        Sp0 += skew(Md*b_ds[i])

    Jp0 = Jsp + Jww_b 

    for i in range(4):
        Jp0 += Md*(np.dot(b_ds[i].T, b_ds[i])*np.eye(3) - np.outer(b_ds[i], b_ds[i]))
    ## Mass equations

    Mt[0:3, 0:3] = M * np.eye(3) # Total mass of the system
    Mt[0:3, 3:6] = -Sp0
    # omega
    Mt[3:6, 0:3] = Sp0
    Mt[3:6, 3:6] = Jp0
    ## Damper equations
    # vP
    for i in range(4):
        Mt[6+i, 6+i] = Md

        # vP
        Mt[0:3, 6+i:6+i+1] = Md*a_ds[i]
        Mt[ 6+i:6+i+1, 0:3] = Md*np.transpose(a_ds[i])
        # omega
        Mt[3:6,  6+i:6+i+1] = Md*skew(b_ds[i])@a_ds[i]
        Mt[6+i:6+i+1, 3:6] = -Md*np.transpose(a_ds[i])@skew(b_ds[i])

    ## Wheel equations

    Mt[10, 10] = Is
    Mt[3:6, 10:11] = Is*a_w
    Mt[10:11, 3:6] = Is*np.transpose(a_w)
    
    return Mt

def update_Mt_Mt_dot(Mt_, Mt_dot_, chi: np.ndarray, chi_dot: np.ndarray):
    # This function updates the Mt and Mt_dot matrices with the current state of the system
    # Constants
    Sp_temp = skew(Ssp + Sw) # Total static moment of the system
    I_temp = Jsp + Jww_b

    Spdot_temp = np.zeros((3,3))
    Idot_temp = np.zeros((3,3))

    # Add damper static moments
    for i in range(4):
        rd = b_ds[i] + chi[i]*a_ds[i]
        rd_dot = chi_dot[i]*a_ds[i]

        Sp_temp += skew(Md*rd)
        I_temp += Md*(np.dot(rd.T, rd)*np.eye(3) - np.outer(rd, rd))

        Spdot_temp += skew(Md*rd_dot)
        Idot_temp += Md*(2*np.dot(rd_dot.T, rd) - (np.outer(rd_dot, rd) + np.outer(rd, rd_dot)))


    # Update Mt
    Mt_[0:3, 3:6] = -Sp_temp
    Mt_[3:6, 0:3] = Sp_temp
    Mt_[3:6, 3:6] = I_temp

    # Update Mt_dot
    Mt_dot_[0:3, 3:6] = -Spdot_temp
    Mt_dot_[3:6, 0:3] = Spdot_temp
    Mt_dot_[3:6, 3:6] = Idot_temp
    
    return Mt_, Mt_dot_

def simulate_nu(t, nu: np.ndarray):
    global Mt, Mt_dot, T_arr
    
    nu = nu.reshape((15,1))
    # Unpack state variables
    vP = nu[0:3, 0:1]
    om = nu[3:6, 0:1]
    chi_dot = nu[6:10,0]
    om_s = nu[10][0]
    chi = nu[11:15,0]

    _nu = nu[0:11] # exclude chi

    # Update mass matrix
    Mt, Mt_dot = update_Mt_Mt_dot(Mt, Mt_dot, chi, chi_dot)


    ## Build nu_dot = inv(Mt-B)@((A@Mt - Mt_dot)@om - c)
    # A
    A = np.zeros((11,11))
    A[0:3, 0:3] = -skew(om)
    A[3:6, 0:3] = -skew(vP)
    A[3:6, 3:6] = -skew(om)

    # B
    B = np.zeros((11,11))
    
    # if a wheel is locked, update B
    if wheel_locked:
        B[10:11, 3:6] = Is*a_w.T 
    
    # c 
    c = np.zeros((11,1))
    c[0:3,0:1] = Fext
    c[3:6, 0:1] = Mext
    for i in range(4):
        rd_i = b_ds[i] + chi[i]*a_ds[i] 
        c[6+i] = Md*om.T@skew(a_ds[i])@(vP - skew(rd_i)@om) - Kd*chi[i] - Cd*chi_dot[i]
    


    if wheel_locked:
        # if wheel is locked, B absorbs the control input
        c[10] = 0
    else:
        c[10] = ga_control(t, _nu)
        
    
    nu_dot = np.zeros((15,1))
    # Compute nu_dot (without chi augmentation)
    nu_dot[0:11, 0:1] = np.linalg.inv(Mt-B)@( ( A@Mt - Mt_dot )@_nu + c )
    # Augment nu_dot with chi_dot
    nu_dot[11:15,0] = chi_dot
    
    # Compute energy
    
    T_temp = 0.5*_nu.T@Mt@_nu 
    T_arr.append(T_temp[0][0])

    return nu_dot.flatten()