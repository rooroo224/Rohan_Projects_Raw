# Author: Rohan Krishna Balaji
# Course : Simulation Science
# Date  : 25.06.2021
# Project : ICTM Analysis, Master's Thesis at Fraunhofer IPT
# Email : rohan.balaji@rwth-aachen.de
# Verion : 1.01

# Script Description: To generate the datasets to perform transformation, apply comensation, cluser, filer and create robost dataset

import transformation
import compensation
import data_imports
import pandas as pd
import os
import glob
import math
import re
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from mpl_toolkits import mplot3d
import numpy as np
from scipy.signal import find_peaks
np.set_printoptions(threshold=np.inf)

# This function can be called to generate a dataset given the blade number and blade angle
def generator(dir_paraquet, dir_planning, dir_final_save, block, angle, tol):
    
    df_m,df_p,compensation_values_df = data_imports.data_out(block,angle,dir_paraquet,dir_planning)
    
    # Just copying machine data points linear(x,y,x) and rotart (a,c) into arrays.
    # These machine points (x,y,x,a,c) are transformed from machine coordinate system to workpiece coordinate system by forward tramsformation

    x = df_m['MachineX'].copy(deep=True)
    y = df_m['MachineY'].copy(deep=True)
    z = df_m['MachineZ'].copy(deep=True)
    a = np.radians(df_m['MachineA']).copy(deep=True)
    c = np.radians(df_m['MachineC']).copy(deep=True)

    size1 = x.shape[0]

    # converting pandas series to numpy array
    x = x.to_numpy()
    y = y.to_numpy()
    z = z.to_numpy()
    a = a.to_numpy()
    c = c.to_numpy()
    
    compensation_values = compensation_values_df.to_numpy()
    
    # Within each cube we have ranges defined in x,y,z for the machine position

    x_range = np.arange(-200,201,100)
    y_range = np.arange(-300,301,150)
    z_range = np.arange(-500,1,50)

    obj3 = compensation.Compensation(compensation_values,x_range,y_range,z_range)

    # Caclculation of compensation error values based on machine positions obtained through inverse transformation
    deltaX1, deltaY1, deltaZ1,deltaI1,deltaJ1,deltaK1 = obj3.calculate(x,y,z)  
    size3 = x.shape[0]
    
    conc3 = np.concatenate((deltaX1.reshape(size3,1),deltaY1.reshape(size3,1),deltaZ1.reshape(size3,1),deltaI1.reshape(size3,1),deltaJ1.reshape(size3,1),deltaK1.reshape(size3,1)),axis=1)
    
    
    # prininting the compensation error values
    df_obj3 = pd.DataFrame(conc3, columns=['deltaX1','deltaY1','deltaZ1','deltaI1','deltaJ1','deltaK1'])
    
    x_compensated = x + deltaX1*10**-3    # since given compensation is to be converted from microns to mm (10**-6 x 10**3 = 10**-3)
    y_compensated = y + deltaY1*10**-3
    z_compensated = z + deltaZ1*10**-3

    conc4 = np.concatenate((x_compensated.reshape(size3,1),y_compensated.reshape(size3,1),z_compensated.reshape(size3,1)),axis=1)
    df_obj4 = pd.DataFrame(conc4, columns=['x_compensated','y_compensated','z_compensated'])

    df_m['compensation_x'] = deltaX1*10**-3
    df_m['compensation_y'] = deltaY1*10**-3
    df_m['compensation_z'] = deltaZ1*10**-3
    
    
    obj = transformation.Transformation(size1,angle)
    # Forward Transformation fuction:
    # Input : Machine points in machine coordinate system
    # Output: returns too tip points and orientation in workpiece coordinate system
    tool_position_workpiece_CS, tool_orientation_workpiece_CS = obj.forward(x_compensated,y_compensated,z_compensated,a,c)

    X = tool_position_workpiece_CS[0,0,:]
    Y = tool_position_workpiece_CS[1,0,:]
    Z = tool_position_workpiece_CS[2,0,:]

    I = tool_orientation_workpiece_CS[0,0,:]
    J = tool_orientation_workpiece_CS[1,0,:]
    K = tool_orientation_workpiece_CS[2,0,:]

    # Verification the correctness of code, i.e on applying reverse transformation on the forward transformation we should get same values or equivalent values i.e (cos(a or c),sin(a or c))
    #*********************************************** Nothing Significant (used for testing) *****************************************
    machine_points_xyz, machine_direction_ac = obj.backward(X,Y,Z,I,J,K)

    x_out = machine_points_xyz[0,0,:]
    y_out = machine_points_xyz[1,0,:]
    z_out = machine_points_xyz[2,0,:]

    a_out = machine_direction_ac[0,0,:]
    c_out = machine_direction_ac[1,0,:]

    conc1 = np.concatenate((x.reshape(size1,1),x_out.reshape(size1,1),y.reshape(size1,1),y_out.reshape(size1,1),z.reshape(size1,1),z_out.reshape(size1,1),a.reshape(size1,1),a_out.reshape(size1,1),c.reshape(size1,1),c_out.reshape(size1,1)),axis=1)

    df_obj1 = pd.DataFrame(conc1, columns=['x','x_out','y','y_out','z','z_out','a','a_out','c','c_out'])
    
    conc11 = np.concatenate((x.reshape(size1,1),y.reshape(size1,1),z.reshape(size1,1),a.reshape(size1,1),c.reshape(size1,1),X.reshape(size1,1),Y.reshape(size1,1),Z.reshape(size1,1),I.reshape(size1,1),J.reshape(size1,1),K.reshape(size1,1)),axis=1)
    dfout11 = pd.DataFrame(conc11, columns=['x','y','z','a','c','X','Y','Z','I','J','K'])
    #***************************************************************************************************************************
    
    tool_tip_X = df_p['Tool Tip Point X'].to_numpy()
    tool_tip_Y = df_p['Tool Tip Point Y'].to_numpy()
    tool_tip_Z = df_p['Tool Tip Point Z'].to_numpy()
    X_inv=tool_tip_X
    Y_inv=tool_tip_Y
    Z_inv=tool_tip_Z
    
    # Here the distances are calculated, for each aquired data point, forward transformation was performed above, and now for each of those poits the disace for all planning points are calculated.
    lst1 = []
    lst2 = []

    dist  = np.zeros(len(tool_tip_X))
    count = 0 
    k = 5

    for i in np.arange(len(X)):   # 47917

        dist = (((tool_tip_X-X[i])**2+(tool_tip_Y-Y[i])**2+(tool_tip_Z-Z[i])**2)**(1/2))
        
        # partition the data into 5 chuncks and find the closest distance points in those chucks, if we use full set at once, we may match up with far way points, which were causing lots of outliers
        pos = np.argpartition(dist, k)
        pos = pos[:k]
        min_val = dist[pos]
        
        # the position closest the iterator is considered, since the far away values are thus avoided from matching 
        pos = pos[np.argmin(abs(pos - i))]
        min_val = dist[pos]
        
        # store both position and the value in two separate lists
        if(abs(min_val<=float(tol))): 
            lst1.append(pos)
            lst2.append(min_val)
            count = count+1
        else:
            lst1.append(np.nan)
            lst2.append(np.nan)
            
    #*******************************************Nothing Significant (used for testing)*******************************************
    (unique,count) = np.unique(lst1, return_counts=True)
    df_freq = pd.DataFrame({'Unique Cluster Points':np.array(unique),'Count':np.array(count)})
    df_freq = df_freq.dropna()
    
    df_cluster = pd.DataFrame({'cluster index in planning data':np.array(lst1),'distance error':np.array(lst2)})
             
    #******************************************************************************************************************************
    # Clustering: for a given planning point, there are none or multiple matching acquried points, those are found in this section of code 
    lst3 = []
    tcp_val = []
    
    # For each planning datapoint, see whih index match with the stored index for minimum distance.
    # Those can be considered as the cluster points and the corresponding distance
    
    for i in np.arange(len(tool_tip_Y)):     # each planning point

        matching = np.where(np.array(lst1)==i)    # see what all machine points matches
        tcp_val.append([lst2[index] for index in matching[0]])          # get the same minmum tcp error distance 
        lst3.append(matching[0])
        
    # for the clusterd points the acquired points are averaged
    lst4 = []
    for i in np.arange(len(lst3)):
        lst4.append(df_m.iloc[list(lst3[i])].mean(axis=0))

    mean_m = pd.concat(lst4,axis=1).T
    
    # distances are averaged
    tcp_avg = [(lambda x: sum(x)/len(x))(item) if len(item)!=0 else np.nan for item in tcp_val]
    
    # Now the final dataframe with plannind data and the corresponding averaged acquired data is obtained 
   
    final_df = pd.concat([df_p,mean_m, pd.DataFrame({'tcp_error':tcp_avg})], axis=1)
    final_df = final_df.drop(['Level','Step'],axis=1)
    final_df = final_df.dropna()
    
    # Inspite of the average, there are some high impuse peaks observed in the data (Run the visualiation.ipnby), these must be removed to obtain a flawless, reliable dataset
    
    def remove_peaks(final_df, compensation_values_df):
        
        # using machine data from newly created combined dataset same forward transformation is performed as discribled above
        
        # Most of the code in this section is for testing and visualization (not printed now), makin takaway is to obtain (X,Y,Z) from (x,y,z) in other words, forward transformation
        x_final = final_df['MachineX'].copy(deep=True)                
        y_final = final_df['MachineY'].copy(deep=True)
        z_final = final_df['MachineZ'].copy(deep=True)
        a_final = final_df['MachineA'].copy(deep=True)
        c_final = final_df['MachineC'].copy(deep=True)

        size1_final = x_final.shape[0]

        # converting pandas series to numpy array
        x_final = x_final.to_numpy()
        y_final = y_final.to_numpy()
        z_final = z_final.to_numpy()
        a_final = a_final.to_numpy(dtype =  np.float64)
        a_final = np.deg2rad(a_final)
        c_final = c_final.to_numpy(dtype =  np.float64)
        c_final = np.deg2rad(c_final)

        tool_tip_X_final = final_df['Tool Tip Point X'].to_numpy()        # using newly created combined dataset
        tool_tip_Y_final = final_df['Tool Tip Point Y'].to_numpy()
        tool_tip_Z_final = final_df['Tool Tip Point Z'].to_numpy()

        X_inv_final = tool_tip_X_final
        Y_inv_final = tool_tip_Y_final
        Z_inv_final = tool_tip_Z_final

        #compensation_values
        compensation_values = compensation_values_df.to_numpy()

        # Within each cube we have ranges defined in x,y,z for the machine position

        x_range = np.arange(-200,201,100)
        y_range = np.arange(-300,301,150)
        z_range = np.arange(-500,1,50)

        obj3_final = compensation.Compensation(compensation_values,x_range,y_range,z_range)

        # Caclculation of compensation error values based on machine positions obtained through inverse transformation
        deltaX1_final, deltaY1_final, deltaZ1_final,deltaI1_final,deltaJ1_final,deltaK1_final = obj3_final.calculate(x_final,y_final,z_final)  
        size3_final = x_final.shape[0]
        conc3_final = np.concatenate((deltaX1_final.reshape(size3_final,1),deltaY1_final.reshape(size3_final,1),deltaZ1_final.reshape(size3_final,1),deltaI1_final.reshape(size3_final,1),deltaJ1_final.reshape(size3_final,1),deltaK1_final.reshape(size3_final,1)),axis=1)

        # prininting the compensation error values
        df_obj3_final = pd.DataFrame(conc3_final, columns=['deltaX1 final','deltaY1 final','deltaZ1 final','deltaI1 final','deltaJ1 final','deltaK1 final'])
        df_obj3_final.head(5)  

        x_compensated_final = x_final + deltaX1_final*10**-3    # since given compensation is to be converted from microns to mm (10**-6 x 10**3 = 10**-3)
        y_compensated_final = y_final + deltaY1_final*10**-3
        z_compensated_final = z_final + deltaZ1_final*10**-3

        conc4_final = np.concatenate((x_compensated_final.reshape(size3_final,1),y_compensated_final.reshape(size3_final,1),z_compensated_final.reshape(size3_final,1)),axis=1)
        df_obj4_final = pd.DataFrame(conc4_final, columns=['x_compensated final','y_compensated final','z_compensated final'])

        obj_final = transformation.Transformation(size1_final,angle)
        # Forward Transformation fuction:
        # Input : Machine points in machine coordinate system
        # Output: returns too tip points and orientation in workpiece coordinate system
        tool_position_workpiece_CS_final, tool_orientation_workpiece_CS_final = obj_final.forward(x_compensated_final,y_compensated_final,z_compensated_final,a_final,c_final)

        X_final = tool_position_workpiece_CS_final[0,0,:]
        Y_final = tool_position_workpiece_CS_final[1,0,:]
        Z_final = tool_position_workpiece_CS_final[2,0,:]

        I_final = tool_orientation_workpiece_CS_final[0,0,:]
        J_final = tool_orientation_workpiece_CS_final[1,0,:]
        K_final = tool_orientation_workpiece_CS_final[2,0,:]

        # Verification the correctness of code, i.e on applying reverse transformation on the forward transformation we should get same values
        machine_points_xyz_final, machine_direction_ac_final = obj_final.backward(X_final,Y_final,Z_final,I_final,J_final,K_final)

        x_out_final = machine_points_xyz_final[0,0,:]
        y_out_final = machine_points_xyz_final[1,0,:]
        z_out_final = machine_points_xyz_final[2,0,:]

        a_out_final = machine_direction_ac_final[0,0,:]
        c_out_final = machine_direction_ac_final[1,0,:]

        conc1_final = np.concatenate((x_final.reshape(size1_final,1),x_out_final.reshape(size1_final,1),y_final.reshape(size1_final,1),y_out_final.reshape(size1_final,1),z_final.reshape(size1_final,1),z_out_final.reshape(size1_final,1),a_final.reshape(size1_final,1),a_out_final.reshape(size1_final,1),c_final.reshape(size1_final,1),c_out_final.reshape(size1_final,1)),axis=1)
        df_obj1_final = pd.DataFrame(conc1_final, columns=['x final','x_out final','y final','y_out final','z final','z_out final','a final','a_out final','c final','c_out'])

        conc11_final = np.concatenate((x_final.reshape(size1_final,1),y_final.reshape(size1_final,1),z_final.reshape(size1_final,1),a_final.reshape(size1_final,1),c_final.reshape(size1_final,1),X_final.reshape(size1_final,1),Y_final.reshape(size1_final,1),Z_final.reshape(size1_final,1),I_final.reshape(size1_final,1),J_final.reshape(size1_final,1),K_final.reshape(size1_final,1)),axis=1)
        dfout11_final = pd.DataFrame(conc11_final, columns=['x final','y final','z final','a final','c final','X final','Y final','Z final','I final','J final','K final'])

        return Y_final, final_df
    
    # Peak removal: by observing the data, it can be seen that the using just one of the components (in this case Y) all the outlers can be eliminated. 
    count = 0
    col = final_df.columns
    while True:
        Y_final, final_df = remove_peaks(final_df,compensation_values_df)
        
        # Two types of spikes are typicalled observed (check visualizaion by runnung  vizualization.ipnyb) one in positive direction and othes in negetive direction, so for robostness both are considered. At the same time, it is not wise to remove the original expected trajectory of the tool path in begining and the end, so especially for the positive direction they are intentionally not removed 
        peaks1, _ = find_peaks(-Y_final, height=(-175,None))
        peaks2,_ = find_peaks(Y_final[0:-1000], height=(240,None))
        peaks = list(peaks1) + list(peaks2)

        if(len(peaks)==0):
            break

        # if each peak on a plateu has to be removed its very expensive, so if two neighbous are consicutively classifies as peaks (given by count), it is assumed as plateu and 100 datapoint sare removed. This is a reasonable compromise in accuracy for significant speedup
        elif(len(peaks)==1 and count>2):
            arr = final_df.to_numpy()
            arr = np.delete(arr, np.arange(peaks[0],peaks[0]+100), 0)
            final_df = pd.DataFrame(arr,columns=col)
            count = 0

        else:
            arr = final_df.to_numpy()
            arr = np.delete(arr,peaks, 0)
            final_df = pd.DataFrame(arr,columns=col)
            count = count+1
            
    # saving data 
    final_df.to_excel(str(dir_final_save)+'finaldf_forward_with_compensation'+str(block)+'__'+str(angle)+'.xlsx',index=False) 

        





           
    
    
    
    
    
    
 
