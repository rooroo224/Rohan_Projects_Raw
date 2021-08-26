# Author: Rohan Krishna Balaji
# Course : Simulation Science
# Date  : 25.06.2021
# Project : ICTM Analysis, Master's Thesis at Fraunhofer IPT
# Email :  rohan.balaji@rwth-aachen.de
# Verion : 1.01
# Credits : Inspired from the forward transformation script provided by Viktor Rudel  
# Script Description: To perform the Forward and Inverse Transformation on the given data

import numpy as np
from scipy.optimize import least_squares
from scipy.optimize import root

class Transformation:
    
    # Initial data, offsets from various pivot points and tool
    
    # Opr : Program origin (G54)
    # Om : Machine origin 
    # Ot : Turret origin 
    # Ob : Rotation center of B  
    # IT_x,IT_y,-T_z - initial tool positions

    OprOm_x = 0
    OprOm_y = 0
    OprOm_z = -510.4830#+167.85

    OtOa_x = 0.006
    OtOa_y =-0.0083  
    OtOa_z =-600.2366

    OtOc_x = 0.000
    OtOc_y = -0.0039
    OtOc_z = -79.7629

    l = 0
    tl = 226.4591

    IT_x=0
    IT_y=0
    IT_z=1

    T_x = 0
    T_y = 0
    T_z = tl-l
    T_CentertoTip = l
    
    def __init__(self, size, angle):
        self.size = size
        self.ones  = np.ones((size,))
        self.zeros = np.zeros((size,))
        self.angle = angle
    
    # All the rotation and translation matices are defined
    
    # Translation in X axis  :    -OtOa_x - T_x
    # Rotation on X axis     :    MachineA data in radians
    # Translation in Z axis  :     OtOc_x
    # Rotation on Z axis     :     MachineC data in radians
    # Translation            :     -OprOm_x + OtOa_x + OtOc_x

    def mat(self,a,c):
        self.init_C = np.array([[np.cos(np.radians(-self.angle))*self.ones, -np.sin(np.radians(-self.angle))*self.ones, self.zeros, self.zeros],
                               [np.sin(np.radians(-self.angle))*self.ones,  np.cos(np.radians(-self.angle))*self.ones, self.zeros, self.zeros],
                               [self.zeros,                        self.zeros,                       self.ones,  self.zeros],
                               [self.zeros,                        self.zeros,                       self.zeros, self.ones]])

        
        self.matrix_translationA = np.array([[self.ones,  self.zeros, self.zeros,  np.ones((self.size,))*(-self.OtOa_x-self.T_x)],
                                       [self.zeros, self.ones,  self.zeros,  np.ones((self.size,))*(-self.OtOa_y-self.T_y)],
                                       [self.zeros, self.zeros, self.ones,   np.ones((self.size,))*(-self.OtOa_z-self.T_z)],
                                       [self.zeros, self.zeros, self.zeros,  self.ones  ]])
        
        self.matrix_rotationA = np.array([[self.ones,   self.zeros,       self.zeros,      self.zeros],
                                          [self.zeros,   np.cos(a),      -np.sin(a),       self.zeros],
                                          [self.zeros,   np.sin(a),       np.cos(a),       self.zeros],
                                          [self.zeros,  self.zeros,        self.zeros     ,self.ones ]],)
         
        self.matrix_translationC = np.array([[self.ones,  self.zeros, self.zeros,  -np.ones((self.size,))*self.OtOc_x],
                                             [self.zeros, self.ones,  self.zeros,  -np.ones((self.size,))*self.OtOc_y],
                                             [self.zeros, self.zeros, self.ones,   -np.ones((self.size,))*self.OtOc_z],
                                             [self.zeros, self.zeros, self.zeros,   self.ones  ]])
        
        self.matrix_rotationC = np.array([[np.cos(c), -np.sin(c), self.zeros, self.zeros],
                                          [np.sin(c),   np.cos(c), self.zeros, self.zeros],
                                          [self.zeros,  self.zeros, self.ones,  self.zeros],
                                          [self.zeros,  self.zeros, self.zeros, self.ones]])
        
        self.matrix_back_translation = np.array([[self.ones,  self.zeros, self.zeros,   np.ones((self.size,))*(-self.OprOm_x + self.OtOa_x + self.OtOc_x) ],
                                                 [self.zeros, self.ones,  self.zeros,   np.ones((self.size,))*(-self.OprOm_y + self.OtOa_y + self.OtOc_y) ],
                                                 [self.zeros, self.zeros, self.ones,    np.ones((self.size,))*(-self.OprOm_z + self.OtOa_z + self.OtOc_z) ],
                                                 [self.zeros, self.zeros, self.zeros,   self.ones                       ]])
        
        return self.init_C, self.matrix_translationA, self.matrix_rotationA,self.matrix_translationC, self.matrix_rotationC, self.matrix_back_translation 
    
    # Forward transformation is performed
    def forward(self,x,y,z,a,c):
        
        if(type(x) != np.ndarray):
            x = x.to_numpy()
            y = y.to_numpy()
            z = z.to_numpy()
            a = a.to_numpy()
            c = c.to_numpy()

        
        self.x = x
        self.y = y
        self.z = z
        self.a = a
        self.c = c
        
        # Tool Position (X,Y,Z )<--  Forward transformation * (x ; y; z ; 1)
        # Tool Orientation (I, J, K) <-- Forward transformation * (0 ;0;  1 ; 0)

        
        # initial tool position 0,0,1
        self.initial_tool_position = np.array([[self.ones*[self.IT_x]],
                                               [self.ones*[self.IT_y]],
                                               [self.ones*[self.IT_z]],
                                               [self.ones*[0]]])
        # initial tool position x,y,x
        self.machine_points_xyz   = np.array([[x],
                                              [y],
                                              [z],
                                              [self.ones]])
        
        self.init_C,self.matrix_translationA, self.matrix_rotationA,self.matrix_translationC, self.matrix_rotationC,self.matrix_back_translation = self.mat(self.a,self.c)

        #Forward transformation <-- Back Translation * Rotation on C axis * Translation in C axis * Rotation on A axis * Translation in A

        self.forward_transformation = np.transpose(self.init_C, (2,0,1)) @ \
                                      np.transpose(self.matrix_back_translation,(2,0,1))@ \
                                      np.transpose(self.matrix_rotationC, (2,0,1)) @ \
                                      np.transpose(self.matrix_translationC, (2,0,1))@ \
                                      np.transpose(self.matrix_rotationA, (2,0,1))@ \
                                      np.transpose(self.matrix_translationA, (2,0,1))
        
        self.tool_position_workpiece_CS = self.forward_transformation @ np.transpose(self.machine_points_xyz,(2,0,1))
        self.tool_position_workpiece_CS = np.transpose(self.tool_position_workpiece_CS,(1,2,0))
        
        self.tool_orientation_workpiece_CS = self.forward_transformation @ np.transpose(self.initial_tool_position,(2,0,1))
        self.tool_orientation_workpiece_CS = np.transpose(self.tool_orientation_workpiece_CS,(1,2,0))
        
        return self.tool_position_workpiece_CS, self.tool_orientation_workpiece_CS
     

    # Method to perform inverse transformation
    def backward(self,X,Y,Z,I,J,K):
        if(type(X) != np.ndarray):
            X = X.to_numpy()
            Y = Y.to_numpy()
            Z = Z.to_numpy()
            I = I.to_numpy()
            J = J.to_numpy()
            K = K.to_numpy()
        
        self.I = I
        self.J = J
        self.K = K
        
        self.X = X
        self.Y = Y
        self.Z = Z
        
        self.tool_tip_points_XYZ =  np.array([[self.X],
                                              [self.Y],
                                              [self.Z],
                                              [self.ones]])
        self.inva = np.zeros(self.size,)
        self.invc = np.zeros(self.size,)
        #####################################################################################################################
        
        # Option 1 : easy but may cause some isues since we have a overdetermined system
        
        self.inva = -np.arccos(self.K)
        self.invc = np.arctan2(-self.I,self.J)
        
        # Option 2 : expensive to calculate, since we mimize a,c angles the solve least squares with unconstrained optimization
        #            but results are definitely more robost. None the less the results are still not accurate 
        
        """
        def inverse_optimiztion(x ,I,J,K):
            return np.array([(I - (np.sin(x[0])*np.sin(x[1]))) , (J + (np.sin(x[0])*np.cos(x[1])) ) , (K - np.cos(x[0]))])
        
        initial_angle = np.array([0,0])
        
        for i in np.arange(self.size):
            res_lsq = least_squares(inverse_optimiztion,initial_angle,loss='soft_l1', args=(self.I[i], self.J[i], self.K[i]))
            self.inva[i] = res_lsq.x[0]
            self.invc[i] = res_lsq.x[1]
            print(i,self.inva[i],self.invc[i])
        """
        
        # Option 3 : Root finding to Solve for least squares with Levenberg-Marquardt, most accurate results are obtained, note that sin/cos(360+x)=sin/cos(x)
        """
        def inv_angles(x, I, J ,K):
            f = [ (I - (np.sin(x[0])*np.sin(x[1]))),
                  (J + (np.sin(x[0])*np.cos(x[1]))),
                  (K - np.cos(x[0]))
                ]
            return f

        for i in np.arange(len(self.I)):
            roots = root(inv_angles,[-1,6], method='lm', args=(self.I[i],self.J[i],self.K[i]), options={'xtol': 1.49012e-15, 'ftol': 1.49012e-15, 'factor':1}) 
            self.inva[i] = roots['x'][0]
            self.invc[i] = roots['x'][1] + np.radians(self.angle)
           #print(i, np.rad2deg(self.inva[i]),np.rad2deg(self.a[i]),np.rad2deg(self.invc[i]),np.rad2deg(self.c[i]))
           #print(i, np.rad2deg(self.inva[i]),np.rad2deg(self.invc[i]))
        """
        ######################################################################################################################
  
        self.init_C, self.matrix_translationA, self.matrix_rotationA,self.matrix_translationC, self.matrix_rotationC,self.matrix_back_translation = self.mat(self.inva,self.invc)
        
        # Backward Transformation <-- matrix_pseudo_inverse(Forward Transformation(using new a,c))
        # Machine linear coordinates(x,y,x) <-- Backward Transformation * tool tip points(X;Y;Z;1)

        
        self.backward_transformation = np.transpose(self.init_C, (2,0,1)) @ \
                                       np.transpose(self.matrix_back_translation,(2,0,1))@ \
                                       np.transpose(self.matrix_rotationC, (2,0,1)) @ \
                                       np.transpose(self.matrix_translationC, (2,0,1))@ \
                                       np.transpose(self.matrix_rotationA, (2,0,1))@ \
                                       np.transpose(self.matrix_translationA, (2,0,1))
                                       
    
        self.backward_transformation = np.linalg.pinv(self.backward_transformation)
        
        self.inv_machine_points_xyz = self.backward_transformation@np.transpose(self.tool_tip_points_XYZ,(2,0,1))
        self.inv_machine_points_xyz = np.transpose(self.inv_machine_points_xyz, (1,2,0))
        
        self.inva = np.reshape(self.inva, (1,1,self.size))
        self.invc = np.reshape(self.invc, (1,1,self.size))
        
        self.machine_direction_ac = np.concatenate((self.inva,self.invc),axis=0)
        
        return self.inv_machine_points_xyz, self.machine_direction_ac
