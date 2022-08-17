# SMARTWATER
This repository contains different scripts and tools that can be used to optimise and design water distribution networks in ATACAMA region (Chile), using Geospatial data with the followign dimensions:

  + Socio-enviromental
  + Topographic
  + Hydraulic
  + Economic

#..............................

02-Optimisation-Gurobi-MO.py 
This script solves a Mixed-Integer Linear Programming Problems (MIP) using a multi-objaective function (MO) with economic and environmental criterea.
This script uses as input an excel file with the right structure and the need Gurobi license.

 More info about how to install Gurobi and how to get a license at:   https://www.gurobi.com/academia/academic-program-and-licenses/
 
 Python:
 Dependencies: pandas
               gurobipy
               numpy                        
               os                        
               math                        
               time                        
               datetime                        
               itertools

#..............................
