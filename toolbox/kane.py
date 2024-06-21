import numpy as np

# Done to use type hints
class Joint:
    pass
class Body:
    pass

class Body:
    def __init__(self, id, mass_center, frame: np.ndarray, mass, inertia: np.ndarray):
        # Ancillary data
        self.id = id

        # Physical data
        self.mass = mass
        self.mass_center = mass_center # In what frame? N?
        self.Ib = inertia # Inertia matrix expressed in body frame
        # Frame data
        self.frame = frame

        
        # Parent
        self.parent = None

        # Children (joints attached to this body)
        self.joints = []

    def add_child(self, joint: Joint, child_body: Body, rb1: np.ndarray, rb2: np.ndarray):
        """
        Adds a child joint and body to the current body.
        Parameters:
            - joint (Joint): The joint connecting the current body to the child body.
            - child_body (Body): The child body being added.
            - rb1 (np.ndarray): The position vector of the joint connection point on the current body.
            - rb2 (np.ndarray): The position vector of the joint connection point on the child body.
        """
        
        joint.set_joint(self, child_body, rb1, rb2)
        child_body.parent = self
        self.joints.append(joint)




class Joint:
    def __init__(self, name, 
                 dof, sequence: str, 
                 axis: np.ndarray, theta_lim: list = None):
        self.name = name

        # Assert dof between 1 and 3
        assert dof >= 1 and dof <= 3, "Degree of freedom must be between 1 and 3 for a rotational joint"
        self.dof = dof
        
        assert len(sequence) == dof, "Sequence must have the same length as the degree of freedom"
        # Be sure that each sequence char is between 1 and 3
        for s in sequence:
            assert s in ['1', '2', '3'], "Sequence must be a string of 1, 2, and 3"
        # Now we can assign the sequence    
        self.sequence = sequence # Sequence of the joint transformation matrix
        
        assert axis.shape == (3,1), "Axis must be a 3x1 vector"
        self.axis = axis

        

        self.body1 = None
        self.rb1 = None # Position of joint in body frame 1 
        self.body2 = None
        self.rb2 = None # Position of joint in body frame 2

        # theta_lim is a list of tuples, each tuple containing the lower and upper limits of the joint
        # if lim0 == lim1, then the joint is free to rotate (not locked, it would not be a joint otherwise)
        assert len(theta_lim) == dof, "Joint limits must be a list of dof elements"
        for lim in theta_lim:
            assert len(lim) == 2, "Each joint limit must be a tuple of two elements"
            assert lim[0] <= lim[1], "Lower limit must be less than or equal to upper limit"
        assert theta_lim[0] < theta_lim[1], "Joint limits must be in ascending order"
        self.theta_lim = theta_lim

        self.T = np.zeros((dof, 3))   # Joint transformation matrix

      

    def set_joint(self, body1: Body, body2: Body, 
                  rb1: np.ndarray, rb2: np.ndarray):
        """
        Sets the joint connection between two bodies.
        
        Parameters:
            - body1 (Body): The inner body.
            - body2 (Body): The outer body.
            - rb1 (np.ndarray): The position vector of the joint connection point on body1 wrt body1 frame.
            - rb2 (np.ndarray): The position vector of the joint connection point on body2 wrt body2 frame.
        """
        self.body1 = body1 
        self.body2 = body2
        self.rb1 = rb1  # inner body, wrt body1
        self.rb2 = rb2  # outer body, wrt body2



            
def simulate(t, u: np.ndarray, body: Body):
        """
        Simulates the body at time t with input u.
        Parameters:
            - t (float): The time at which to simulate the body.
            - u (np.ndarray): The input to the body.
            - body (Body): The root body to simulate.
        
        

        """


        # Evaluate joint partials and their time derivatives
        Ti = evaluate_joint_partials(body) # List of joint transformation matrices, one for each joint
        T_dot = evaluate_joint_partials_dot(body) # List of joint transformation matrices, one for each joint

        # For each Bi, map {x, u} to {om1, vp1, x, RN2i*r_pi}

        # Find all relative rotation matrices Rj2i

        # Find all relative positions r_pj2i	

        # Evaluate OM and V and find omi, and vpi for all bodies

        # Evaluate remainder accelerations

        # Evaluate inertia forces and torques for each Bi

        # Evaluate external forces and torques for each Bi (including interaction forces between bodies)

        # Evaluate Kane equations [OM^t*[J]*OM + V^t*[M]*V ] u_dot = OM^t*{T- J*alpha_r - [skew(om)*J*om]} + V^t*{F - M*a_r}	
        # solving for u_dot

        # Check locking conditions for all joints (eliminate the corresponding u_dot component) and solve for u_dot again

        # Use u to integrate also x to the next step, along with rp1


