# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 15:45:52 2022

This script solves a Mixed-Integer Linear Programming Problems (MIP) for the 
optimisation of a water distribution Network for water supply in Chile

More info about this at:
    https://docs.python-mip.com/en/latest/intro.html
    
@author: njame
"""

###############################################################################
#%%                                                   PACKAGES YOU NEED TO LOAD
###############################################################################

import os
import time
import math

from tqdm import tqdm

from itertools import product
import numpy as np
import pandas as pd
from sys import stdout as out
#from mip import Model, xsum, minimize, BINARY, CBC, OptimizationStatus
from mip import *

###############################################################################
#  ^    ^    ^    ^    ^    ^    ^    ^    ^   EMD OF PACKAGES YOU NEED TO LOAD
###############################################################################

###############################################################################
#                                                                List of inputs
###############################################################################

time_before_execution = time.time()

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
#..............................................................................

#load input file with network data:
#..............................................................................
location_input_file= r'C:\Pypy\Input_files'
location_output_files= r'C:\Pypy\Output_files'
suffix_results = r'41_nodes'

name_input_file = r'Inputs_41_nodos_scenario.xlsx' # 41 nodes
#name_input_file=r'Inputs_208_nodos_scenario.xlsx' # 208 nodes

input_file_path =os.path.join(location_input_file,name_input_file)

dist = pd.read_excel(input_file_path,sheet_name='distance')
alt  = pd.read_excel(input_file_path,sheet_name='altitude')
cap  = pd.read_excel(input_file_path,sheet_name='capacity')
dem  = pd.read_excel(input_file_path,sheet_name='demand')
#..............................................................................

###############################################################################
#    ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  List of inputs
###############################################################################

###############################################################################
#                                 Pre - run 
###############################################################################
dist=dist.values.tolist()
df_dist = pd.DataFrame(dist)
df_dist=df_dist.replace(np.NaN, 0)
dist=df_dist.values.tolist()

alt=alt.values.tolist()
df_alt = pd.DataFrame(alt)
df_alt=df_alt.replace(np.NaN, 0)
alt=df_alt.values.tolist()
numpy_array = np.array(alt)
transpose = numpy_array.T
alt = transpose.tolist()
#print(alt)

cap=cap.values.tolist()
#dict = cap.to_dict()
cap=cap[0]

dem=dem.values.tolist()
dem=dem[0]

###############################################################################
#         ^       ^      ^     pre- run procedures    ^       ^      ^    
###############################################################################

###############################################################################
#                                                                        Main() 
###############################################################################    

n, V = len(dist[0]), set(range(len(dist[0])))
#print(n,V)

nc, Vc = len(cap), set(range(len(cap)))
#print(nc,Vc)

#nd, Vd = len(dem), set(i for i in range(len(dist[0])-len(cap)))
nd, Vd = len(dem), set(i for i in range(max(Vc)+1,max(V)+1))
#print(nd,Vd)

Vz={24,33,35,36,37,38,39,40}

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

model = Model(solver_name=CBC)

y = [[model.add_var('y({},{})'.format(i, j),var_type=BINARY) for i in V] for j in V]


F = [[model.add_var('F({},{})'.format(i, j),lb=0,ub=5,var_type=CONTINUOUS) for i in V] for j in V]

print(r'Minimizing Objective Function, please wait...')
   
model.objective = minimize(xsum((cost[i][j]*F[i][j]*20*365*24*Ecost/1000000+capitalcost[i][j]*y[i][j]*Ccost) for i in V for j in V))

for i in Vc:   
    model += capdict[i]>=xsum(F[i][j] for j in V-{i})-xsum(F[j][i] for j in V-{i}) 

for i in Vd:      
    model += xsum(F[j][i] for j in V-{i})==xsum(F[i][j] for j in V-{i})+demdict[i]  
    
for i in V:
    for j in V-{i}: 
        model+= F[i][j]<=BigM*y[i][j]
        model+= F[i][j]>=Min*y[i][j]       

for i in Vz:   
    for j in V-{i}:         
        model+= y[i][j]==0
        model+= y[j][i]==0
     
model.optimize()

# Print results:
#..............................................................................

print("Solution {} found.".format(round(model.objective_value,4)))
print("Solver status {} ".format(model.status))
print('The model has {} vars, {} constraints and {} nzs'.format(model.num_cols, model.num_rows, model.num_nz))

for (i, j) in product(V, V):
    if y[i][j].x >= 0.99:
        print("Start {}: End {}: Flowrate {}".format(i, j, round(F[i][j].x,4)))

#..............................................................................

# Save resutls in txt files:
#..............................................................................

path_description_list =[]
Flow_rates_list =[]

for (i, j) in product(V, V):
    if y[i][j].x >= 0.99:
        Path_description = r'From-' + str(i) + r'-to-' + str(j)
        flow_rate= round(F[i][j].x,4)
        path_description_list.append(Path_description) 
        Flow_rates_list.append(flow_rate)

FP_path_descriptions = os.path.join(location_output_files,("Path_descriptions_" + suffix_results + r'.txt'))

TXT_file = open(FP_path_descriptions, "w")
for element in path_description_list:
    TXT_file.write(element + "\n")
TXT_file.close()

FP_flow_rates = os.path.join(location_output_files,("Flow_rates_" + suffix_results + r'.txt'))

TXT_file = open(FP_flow_rates, "w")
for element in Flow_rates_list:
    TXT_file.write(str(element) + "\n")
TXT_file.close()

#..............................................................................

elapsed_time = (time.time() - time_before_execution)
Fraction_of_Seconds, Seconds =math.modf(elapsed_time)
Fraction_of_hours, hours =math.modf(Seconds/3600)

print('Total execution time : ' + str(round(hours)) + ' hours ' + str(round(Seconds/60)%60)+ ' minutes ' + str(round(Seconds%60))+' seconds' )
 
###############################################################################
#%%                                                                      Main() 
###############################################################################    

