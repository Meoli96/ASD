import numpy as np

# Done to use type hints
class Joint:
    pass
class Body:
    pass

class Body:
    
    def __init__(self, bid, mass_center, frame: np.ndarray, mass, inertia: np.ndarray):
        # Ancillary data
        self.id = bid

        # Physical data
        self.mass = mass
        self.mass_center = mass_center # In what frame? N?
        self.Ib = inertia # Inertia matrix expressed in body frame
        # Frame data
        self.frame = frame

        
        # Joint parent 
        self.jparent = None

        # Children (joints attached to this body)
        self.jchilds = []

        def dof_n(self):
            # Returns the number of degrees of freedom of the body 
            # and counts the number of bodies in the system
            
            # Base case (no children)
            dof = 0
            n = 1

            # Recursive case
            for joint in self.jchilds:
                # Make recursive call   
                dof_c, n_c = joint.b_child.dof()
                dof += joint.dof + dof_c
                n += n_c
            return dof, n




class Joint:
    def __init__(self, 
                 dof, sequence: str, 
                 axis: np.ndarray, theta_lim: list[tuple] = None):
        self.b_parent = None
        self.b_child = None

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

        

        self.rb1 = None # Position of joint in body frame 1 
        self.rb2 = None # Position of joint in body frame 2

        # theta_lim is a list of tuples, each tuple containing the lower and upper limits of the joint
        # if lim0 == lim1, then the joint is free to rotate (not locked, it would not be a joint otherwise)
        assert len(theta_lim) == dof, "Joint limits must be a list of dof elements"
        for lim in theta_lim:
            assert len(lim) == 2, "Each joint limit must be a tuple of two elements"
            assert lim[0] <= lim[1], "Lower limit must be less than or equal to upper limit"
        self.theta_lim = theta_lim

        self.PHI = np.zeros((3, dof))   # Joint transformation matrix

      

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



def connect(joint: Joint, b_parent: Body, b_child: Body):
    """
    Connects a joint to two bodies.
    Parameters:
        - joint (Joint): The joint to connect.
        - b_parent (Body): The parent body.
        - b_child (Body): The child body.
    """
    joint.b_parent = b_parent
    joint.b_child = b_child

    b_parent.joints.append(joint)
    
    if b_child.jparent is not None:
        b_child.parent = (b_parent,) + b_parent.jparent
    else:
        b_child.parent = (b_parent,)
    
            


def evaluate_OM_V(root: Body, u: np.ndarray):
    # Compute overall dof and number of bodies of the system
    dof, n = root.dof() + 6 # 6 for the base body
    # Initialize the overall OM and V matrices


def simulate(t, u: np.ndarray, root: Body):
        """
        Simulates the body at time t with input u.
        Parameters:
            - t (float): The time at which to simulate the body.
            - u (np.ndarray): The input to the body.
            - body (Body): The root body to simulate.
        
        

        """


        # Evaluate joint partials and their time derivatives
        OM, V = evaluate_OM_V(root, u)

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


def esploraTree(body: Body):
