# Author: Rohan Krishna Balaji
# Course : Simulation Science
# Date  : 25.06.2021
# Project : ICTM Analysis, Master's Thesis at Fraunhofer IPT
# Email :  rohan.balaji@rwth-aachen.de
# Verion : 1.01

# Script Description: To load the data needed from the given location

import numpy as np
import pandas as pd

def data_out(block,angle,path_parquet,path_planning):
    
    machine_data_file = str(path_parquet)+'ProgNumber-'+str(block)+'--Blade--'+str(angle)+'_DownsampledData.parquet'
    df_m = pd.read_parquet(machine_data_file)


    planning_file = str(path_planning)+'Planning_Data_OP'+str(block)+'.xlsx'
    df_p = pd.read_excel(planning_file)
    
    compensation_values_df = pd.read_csv('Compensation_Makino.csv')
    
    #print('Files loaded Succussfully!')
    
    return df_m,df_p,compensation_values_df