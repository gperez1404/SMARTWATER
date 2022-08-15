# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 13:08:30 2022

@author: uqgpere2
"""

###############################################################################
#%%                                                   PACKAGES YOU NEED TO LOAD
###############################################################################

import os
import time
import math

from itertools import product
import numpy as np
import pandas as pd
from sys import stdout as out
from mip import *
import gurobipy as grb
from gurobipy import GRB

import time
from datetime import datetime

###############################################################################
#  ^    ^    ^    ^    ^    ^    ^    ^    ^   EMD OF PACKAGES YOU NEED TO LOAD
###############################################################################

###############################################################################
#                                                                Define  inputs
###############################################################################

# Inputs :
#..............................................................................

Ecost=50 #US$_MWhr
Ccost=0.0013 #MillionUS$_m typical pipeline cost
Min=0
BigM=10000
v=1 #typical pipe vel m per s 
C=140 #Chezy coeff for smooth concrete pipe
R=0.3 #typical pipe radius in m
Min=0
BigM=1000

# location of input excel file with the network:
    
input_file_path  = r'D:\REPOSITORIES_SMI_828M6Z2\SMI-828M6Z2-Python-Scripts\Inputs_Scripts\Inputs_Gurobi_208_nodes.xls'

# Define location to save resutls:
location_output_files= r'R:\03_GISdata\data_Gabriel\04-Tools\02-Python-Scripts\03-Output-files-R'
suffix_results = r'test_desktop'

# This is the % of gap allowed in the solution of the model:
Opti_gap= 0.15

Number_of_threads = 6
Time_limit_hours  = 1

#optional paramters:

# list of nodes to force 0 in the solution
# This set defines the nodes that will be de-activated
Vz={}
    
#..............................................................................

#Load excel file with network data:
#..............................................................................

dist = pd.read_excel (input_file_path,sheet_name='distance')
dist=dist.values.tolist()
df = pd.DataFrame(dist)
df=df.replace(np.NaN, 0)
dist=df.values.tolist()

alt = pd.read_excel (input_file_path,sheet_name='altitude')
alt=alt.values.tolist()
df = pd.DataFrame(alt)
df=df.replace(np.NaN, 0)
alt=df.values.tolist()

cap = pd.read_excel (input_file_path,sheet_name='capacity')
cap=cap.values.tolist()
cap=cap[0]

dem = pd.read_excel (input_file_path,sheet_name='demand')
dem=dem.values.tolist()
dem=dem[0]
#..............................................................................

###############################################################################
#     ^     ^     ^     ^     ^     ^     ^     ^     ^        Define  inputs
###############################################################################

###############################################################################
#                                             Define variables and parameters
###############################################################################

time_before_execution = time.time()

n, V = len(dist[0]), set(range(len(dist[0])))
#print(n,V)

nc, Vc = len(cap), set(range(len(cap)))
#print(nc,Vc)

#nd, Vd = len(dem), set(i for i in range(len(dist[0])-len(cap)))
nd, Vd = len(dem), set(i for i in range(max(Vc)+1,max(V)+1))
#print(nd,Vd)


capdict=dict(zip(Vc,cap))

demdict=dict(zip(Vd,dem))

hff=((v*(3.1416*R**2))**1.852)/((0.849**1.852)*(C**1.852)*((R/2)**1.167))
cte=9.81*1000/1000000
hf = [[hff* j for j in sub] for i, sub in enumerate(dist)]

new_list=[]

for i in range(len(dist)):   
  list11=hf[i]
  list22=alt[i]
  
  sum_list = [a + b for a, b in zip(list11, list22)]  
  new_list.append(sum_list)
  
cost = [[cte* j for j in sub] for i, sub in enumerate(new_list)]
#print(cost)

df = pd.DataFrame(cost)
df=df.replace(np.NaN, 0)
cost=df.values.tolist()

capitalcost=dist

# Create TXT File to print results:
#..............................................................................   

FP_path_model_results = os.path.join(location_output_files,("01_Res_Gurobi_" + suffix_results + r'.txt'))
print(r'saving log at:' + FP_path_model_results)

TXT_file = open(FP_path_model_results, "w")
line= r'Scenario: ' + suffix_results
TXT_file.write(line + "\n")
line= r'Number of Threads to use: ' + str(Number_of_threads)
TXT_file.write(line + "\n")
line=  r'Minimizing Objective Function'
TXT_file.write(line + "\n")
line= r'Process started at:  ' + str( datetime.today().strftime('%Y-%m-%d %H:%M:%S'))
TXT_file.write(line + "\n")
TXT_file.flush()
#..............................................................................


###############################################################################
#      ^      ^      ^      ^      ^      ^     Define variables and parameters
###############################################################################

###############################################################################
#                                                          Optimisation problem
###############################################################################

# Define Optimisation problem:
#..............................................................................

#model = Model(solver_name=CBC)
model = Model(solver_name=GUROBI)
grb.setParam('MIPGap', Opti_gap)
grb.setParam(GRB.Param.Threads, Number_of_threads)
grb.setParam(GRB.Param.TimeLimit, (Time_limit_hours*3600))


y = [[model.add_var('y({},{})'.format(i, j),var_type=BINARY) for i in V] for j in V]


F = [[model.add_var('F({},{})'.format(i, j),lb=0,ub=5,var_type=CONTINUOUS) for i in V] for j in V]
   
model.objective = minimize(xsum((cost[i][j]*F[i][j]*20*365*24*Ecost/1000000+capitalcost[i][j]*y[i][j]*Ccost) for i in V for j in V))

for i in Vc:   
    model += capdict[i]>=xsum(F[i][j] for j in V-{i})-xsum(F[j][i] for j in V-{i}) 

for i in Vd:      
    model += xsum(F[j][i] for j in V-{i})==xsum(F[i][j] for j in V-{i})+demdict[i]  
    
for i in V:
    for j in V-{i}: 
        model+= F[i][j]<=BigM*y[i][j]
        model+= F[i][j]>=Min*y[i][j]       

if(len(Vz)> 0):
    for i in Vz:
        for j in V-{i}:
            model+= y[i][j]==0
            model+= y[j][i]==0
     
#..............................................................................




# Here you launch the optimsiation method
#..............................................................................
print(r'Minimizing Objective Function, please wait...')
print(r'Check the other console to see progress')

model.optimize()

elapsed_time = (time.time() - time_before_execution)
Fraction_of_Seconds, Seconds =math.modf(elapsed_time)
Fraction_of_hours, hours =math.modf(Seconds/3600)

print('Total execution time Optimsiation algorithm: ' + str(round(hours)) + ' hours ' + str(round(Seconds/60)%60)+ ' minutes ' + str(round(Seconds%60))+' seconds' )

line= r'Optimization process finished at:  ' + str( datetime.today().strftime('%Y-%m-%d %H:%M:%S'))
TXT_file.write(line + "\n")

#..............................................................................

###############################################################################
#     ^     ^     ^     ^     ^     ^     ^     ^     ^   Optimisation problem
###############################################################################


# Print results in console:
#..............................................................................

print("Solution {} found.".format(model.objective_value))
print("Solver status {} ".format(model.status))
print('model has {} vars, {} constraints and {} nzs'.format(model.num_cols, model.num_rows, model.num_nz))

for (i, j) in product(V, V):
    if y[i][j].x >= 0.99:
        print("Start {}: End {}: Flowrate {}".format(i, j, F[i][j].x))    
  
#..............................................................................

# Save results as txt file:
#..............................................................................

line = "Status: " + str(model.status)
TXT_file.write(line + "\n")
line = "Total Cost = " + str(model.objective_value)
TXT_file.write(line + "\n")
line = 'model has {} vars, {} constraints and {} nzs'.format(model.num_cols, model.num_rows, model.num_nz)
TXT_file.write(line + "\n")

for (i, j) in product(V, V):
    if y[i][j].x >= 0.99:
        strat_node= str(i)
        end_node  = str(j)
        Flow_Rate= str(round(F[i][j].x,3))  
        line =  r'Optimal flow rate from node:' + strat_node + r' to node: ' + end_node + " = " +  Flow_Rate  
        TXT_file.write( line  + "\n")


elapsed_time = (time.time() - time_before_execution)
Fraction_of_Seconds, Seconds =math.modf(elapsed_time)
Fraction_of_hours, hours =math.modf(Seconds/3600)

print('Total execution time : ' + str(round(hours)) + ' hours ' + str(round(Seconds/60)%60)+ ' minutes ' + str(round(Seconds%60))+' seconds' )
line = 'Total execution time : ' + str(round(hours)) + ' hours ' + str(round(Seconds/60)%60)+ ' minutes ' + str(round(Seconds%60))+' seconds'
TXT_file.write(line + "\n")
TXT_file.close()

#..............................................................................

# Save results as excel file:
#..............................................................................    
path_description_list =[]
Flow_rates_list =[]

for (i, j) in product(V, V):
    if y[i][j].x >= 0.99:
        Path_description = r'From-' + str(i) + r'-to-' + str(j)
        flow_rate= round(F[i][j].x,4)
        path_description_list.append(Path_description) 
        Flow_rates_list.append(flow_rate)

# create an excel file with results:

df_results = pd.DataFrame({'PathDes': path_description_list,'Qm3s': Flow_rates_list})    

excel_filename =  "00_Res_Gurobi_" + suffix_results + r'.xls'
excel_filepath = os.path.join(location_output_files,excel_filename)

with pd.ExcelWriter(excel_filepath) as writer:
    df_results.to_excel(writer)
    
#..............................................................................    

###############################################################################
#%%                                                                         End 
###############################################################################