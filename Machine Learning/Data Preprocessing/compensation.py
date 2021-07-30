# Author: Rohan Krishna Balaji
# Course : Simulation Science
# Date  : 25.06.2021
# Project : ICTM Analysis, Master's Thesis at Fraunhofer IPT
# Email :  rohan.balaji@rwth-aachen.de
# Verion : 1.01

# Script Description: To apply compenstion on the supplied machine points (acqired)

import numpy as np

class Compensation:

    # there are 5 points in each direction (x,y) plane but in z direction it is calculated as remaining points
    x = 5
    y = 5
    
    def __init__(self,compensation_values, x_range, y_range, z_range):
        self.compensation_values = compensation_values
        self.x_range = x_range
        self.y_range = y_range
        self.z_range = z_range
    
    def calculate(self,xm,ym,zm):
        # Machine points data from paraquet files
        self.xm = xm  
        self.ym = ym
        self.zm = zm
      
        self.size = self.xm.shape[0]
        
        # creating empty arrays to store the compensation error values bases on machine points
        self.deltaX1 = np.zeros(self.size)
        self.deltaY1 = np.zeros(self.size)
        self.deltaZ1 = np.zeros(self.size)
        self.deltaI1 = np.zeros(self.size)
        self.deltaJ1 = np.zeros(self.size)
        self.deltaK1 = np.zeros(self.size) 
        
        # checking in which range - (x1, x2) the point falls when projected on a single axis seperately
        for ii in np.arange(self.size):    # Calculation of compensation error for each machine point
            
            a1 = self.xm[ii]
            a2 = self.ym[ii]
            a3 = self.zm[ii]
            
            for i in np.arange(len(self.x_range)-1):     # loops from 0 to 4 
                if (a1 >= self.x_range[i] and a1 <= self.x_range[i+1]):  
                    # if the given machine data point is in between two consicutive points of intervel, we can obtain
                    # x - axis position of first and second points of a cube 
                    x1 = i
                    x2 = i+1
                    break
            
            for i in np.arange(len(self.y_range)-1):
                if (a2 >= self.y_range[i] and a2 <= self.y_range[i+1]):
                    # if the given machine data point is in between two consicutive points of intervel, we can obtain
                    # y - axis position of first and second points of a cube 
                    y1 = i
                    y2 = i+1
                    break
                    
            for i in np.arange(len(self.z_range)-1):
                if (a3 >= self.z_range[i] and a3 <= self.z_range[i+1]):
                    # if the given machine data point is in between two consicutive points of intervel, we can obtain
                    # z - axis position of first and second points of a cube 
                    z1 = i
                    z2 = i+1
                    break
            # At this point we have the x,y,z points assiged to a starting points x1,y1,z1 of cube
            # In other words, we identify the position of cube that fits our machine data point
            
            # but these points are individual projection on x,y,z axis respectively. But to get a cube from this we need to 
            # use relation between them to obtain all the neighbouring positions around the cube
                    
            # calculation of neighbourhood indices     
                
            i1 = x1 + self.x*y1 + (self.x*self.y)*z1 
            # if we visualize the cube the first point would be the x projection and the multiplication of y with number of point in x
            # and number of suxh x,y planes in z axis
            
            # rest of the corners are arund the first point
            i2 = i1 + 1
            i4 = i1 + self.x
            i3 = i4 + 1
            i5 = i1 + (self.x*self.y)  # one x,y plane above first corner
            i6 = i5 + 1
            i8 = i5 + self.x
            i7 = i8 + 1
            
            # neighbourhood compensation data from the calculated indices of cube corners
            # Assuming given 
            p1 = self.compensation_values[i1]
            p2 = self.compensation_values[i2]
            p3 = self.compensation_values[i3]
            p4 = self.compensation_values[i4]
            p5 = self.compensation_values[i5]
            p6 = self.compensation_values[i6]
            p7 = self.compensation_values[i7]
            p8 = self.compensation_values[i8]
            
            # from here using interpolation formula literally
        
            ratiox =  (a1 - self.x_range[x1])/(self.x_range[x2] - self.x_range[x1])
            ratioy =  (a2 - self.y_range[y1])/(self.y_range[y2] - self.y_range[y1])
            ratioz =  (a3 - self.z_range[z1])/(self.z_range[z2] - self.z_range[z1])
            
            # interpolation for [ΔX1, ΔY1 ,ΔZ1] due to pml  
            
            self.deltaX1[ii] = p1[0]*(1-ratiox)*(1-ratioy)*(1-ratioz) + p2[0]*ratiox*(1-ratioy)*(1-ratioz) \
                              + p3[0]*ratiox*ratioy*(1-ratioz) + p4[0]*(1-ratiox)*ratioy*(1-ratioz) \
                              + p5[0]*(1-ratiox)*(1-ratioy)*ratioz + p6[0]*ratiox*(1-ratioy)*ratioz \
                              + p7[0]*ratiox*ratioy*ratioz + p8[0]*(1-ratiox)*ratioy*ratioz
            
            self.deltaY1[ii] = p1[1]*(1-ratiox)*(1-ratioy)*(1-ratioz) + p2[1]*ratiox*(1-ratioy)*(1-ratioz) \
                              + p3[1]*ratiox*ratioy*(1-ratioz) + p4[1]*(1-ratiox)*ratioy*(1-ratioz) \
                              + p5[1]*(1-ratiox)*(1-ratioy)*ratioz + p6[1]*ratiox*(1-ratioy)*ratioz \
                              + p7[1]*ratiox*ratioy*ratioz + p8[1]*(1-ratiox)*ratioy*ratioz

            self.deltaZ1[ii] = p1[2]*(1-ratiox)*(1-ratioy)*(1-ratioz) + p2[2]*ratiox*(1-ratioy)*(1-ratioz) \
                              + p3[2]*ratiox*ratioy*(1-ratioz) + p4[2]*(1-ratiox)*ratioy*(1-ratioz) \
                              + p5[2]*(1-ratiox)*(1-ratioy)*ratioz + p6[2]*ratiox*(1-ratioy)*ratioz \
                              + p7[2]*ratiox*ratioy*ratioz + p8[2]*(1-ratiox)*ratioy*ratioz
            
            # interplation for  [ΔI1 ,ΔJ1, ΔK1] based on pml
            self.deltaI1[ii] = p1[3]*(1-ratiox)*(1-ratioy)*(1-ratioz) + p2[3]*ratiox*(1-ratioy)*(1-ratioz) \
                           + p3[3]*ratiox*ratioy*(1-ratioz) + p4[3]*(1-ratiox)*ratioy*(1-ratioz) \
                           + p5[3]*(1-ratiox)*(1-ratioy)*ratioz + p6[3]*ratiox*(1-ratioy)*ratioz \
                           + p7[3]*ratiox*ratioy*ratioz + p8[3]*(1-ratiox)*ratioy*ratioz

            self.deltaJ1[ii] = p1[4]*(1-ratiox)*(1-ratioy)*(1-ratioz) + p2[4]*ratiox*(1-ratioy)*(1-ratioz) \
                           + p3[4]*ratiox*ratioy*(1-ratioz) + p4[4]*(1-ratiox)*ratioy*(1-ratioz) \
                           + p5[4]*(1-ratiox)*(1-ratioy)*ratioz + p6[4]*ratiox*(1-ratioy)*ratioz \
                           + p7[4]*ratiox*ratioy*ratioz + p8[4]*(1-ratiox)*ratioy*ratioz

            self.deltaK1[ii] = p1[5]*(1-ratiox)*(1-ratioy)*(1-ratioz) + p2[5]*ratiox*(1-ratioy)*(1-ratioz) \
                          + p3[5]*ratiox*ratioy*(1-ratioz) + p4[5]*(1-ratiox)*ratioy*(1-ratioz) \
                          + p5[5]*(1-ratiox)*(1-ratioy)*ratioz + p6[5]*ratiox*(1-ratioy)*ratioz \
                          + p7[5]*ratiox*ratioy*ratioz + p8[5]*(1-ratiox)*ratioy*ratioz
            
        return self.deltaX1,self.deltaY1,self.deltaZ1,self.deltaI1,self.deltaJ1,self.deltaK1