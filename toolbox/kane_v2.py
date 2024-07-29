import numpy as np
class Body:
    def __init__(self, bid):
        self.id = bid

        self.jparents = None
        self.jchilds = []

    def dof(self):
        # self.dof() returns the degrees of freedom of the system 
        # and the number of joints of the system

        
        dof = 0
        
        for joint in self.jchilds:
            dofc = joint.child.dof()
            dof += joint.dof + dofc
           
        if self.isroot():
            # Add 6 degrees of freedom for the root body
            dof += 6

        return dof
    
    def gen_omega(self, q: np.ndarray):
        # Generate the angular velocity of the tree
        # Need to associate an index at each body of the tree
        # The index is the position of the body in the state vector

        
        if self.isroot():
            # Root body
            pass
        # retrieve childs from joints
       

        # Associate an index to each body and joint

        for joint in self.jchilds:
            child_body = joint.child

    def retrieve_joints(self):
        # Retrieve all the joints in the tree, row by row

        joints = [i for i in self.jchilds]
        for joint in self.jchilds:
            joints += joint.child.retrieve_joints()
        return joints
    
    def retrieve_bodies(self):
        # Retrieve the list of body ids in BFS order
        bodylist = [self]
        for joint in self.jchilds:
            bodylist += joint.child.retrieve_bodies()
        return bodylist
    
    def isroot(self):
        return self.jparents is None



class Joint:
    def __init__(self, 
                 dof: int, sequence: str = "", 
                 axis: np.ndarray = np.zeros([3,1]),  # defaults to empty column vector
                 theta_lim: list[tuple] = None):
        
        # Body objects link is attached to
        self.parent = None
        self.rb_parent = None # Distance from attach point in parent body frame
        
        self.child = None
        self.rb_child = None # Distance from attach point in child body frame

        
        # dof validity check
        assert dof >= 1 and dof <= 3, "Degree of freedom must be between 1 and 3 for a rotational joint"
        self.dof = dof
        

        # theta_lim validation 
        if theta_lim != None:
            assert len(theta_lim) == dof, "Joint limits must be a list of dof elements"
            for q_lim in theta_lim:
                assert len(q_lim) == 2, "Each joint limit must be a tuple of two elements"
                assert q_lim[0] <= q_lim[1], "Lower limit must be less than or equal to upper limit"
            self.theta_lim = theta_lim

    def phi(self, q: np.ndarray):
        if self.dof == 1:
            # One dimensional 
            
            pass
        elif self.dof == 2:
            pass
        else:
            # dof  == 3
            pass

    
    def connect(self,
                parent: Body, child: Body,
                rb_parent:np.ndarray = np.zeros((3,1)), 
                rb_child: np.ndarray= np.zeros((3,1))):
        
        ## Set joint attributes
        # Parent
        self.parent = parent
        self.rb_parent = rb_parent
        # Child
        self.child = child
        self.rb_child = rb_child

        ## Add joint to parent childs
        parent.jchilds.append(self)
        ## Add joint to child parent joints
        if parent.jparents is None:
            child.jparents = (self,)
        else:
            child.jparents = parent.jparents + (self,)
        return self

def connect(joint: Joint, parent: Body, rb_parent:np.ndarray ,child: Body, rb_child: np.ndarray):
    ## Set joint attributes
    # Parent
    joint.parent = parent
    joint.rb_parent = rb_parent
    # Child
    joint.child = child
    joint.rb_child = rb_child
    
    ## Add joint to parent childs
    parent.jchilds.append(joint)
    ## Add joint to child parent joints
    if parent.jparents is None:
        child.jparents = (joint,)
    else:
        child.jparents = parent.jparents + (joint,)

def isroot(body: Body):
    # Helper function to check if a body is the root body
    return body.jparents is None

def simulate_tree(root_body: Body):
    # Extract the degrees of freedom of the system and the number of joints
    assert isinstance(root_body, Body) and isroot(root_body), "The root body must be the root of the tree"
    
    dof, nj = root_body.dof()

    # Create the state vector
    q = np.zeros((dof, 1)) # [om_b, q1, ... , qn, v_b]

    # We need to associate at each joint a number of state variables based on the joint dof
    




def simulate(joints: list[Joint], q: np.ndarray, qdot: np.ndarray, qddot: np.ndarray):
    # Simulate the system
    for joint in joints:
        # Compute the joint transformation matrix
        joint.phi(q)
        # Compute the joint velocity matrix
        joint.phi_dot(q, qdot)
        

def evaluate_OM_V(root: Body, u: np.ndarray):
    