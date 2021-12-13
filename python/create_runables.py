#########################################
# disclaimer: this script is used to run 
# several configs and to find a good 
# config to align the SciFi
# 
# for explicit use of the bash script to 
# run specific configs, comment out line 12
# where the random configs are created and instead
# MANUALLY write the constr_.txt and dofs_.txt files !
#########################################

import os
import numpy as np
import pandas as pd
import random

# these are the full data tables
dofs = pd.read_csv("dofs.txt", header=None)
dofs.columns = ['dof']
dof = dofs['dof']
length_dof = len(dof)

constr = pd.read_csv("constraints.txt", header=None)
constr.columns = ['constr']
con = constr['constr']
length_con = len(con)

n_dofs = os.environ['n_dofs']  # number of dofs i want to check
n_runfiles = os.environ['n_runfiles']  # number of runfile pairs i want to create for each n_dofs
n_constraints = os.environ['n_constraints']  # number of random constraints for each alignment job

# first iterate through the number of random dofs i want to use
rd_dof = random.sample(range(length_dof), int(n_dofs))  # chooses position in array of dofs

# save the dof names for the bash script to find files
file_dofnames = open("runlists/dofnames.txt", "a")
for i in rd_dof:
    file_dofnames.write(dof[i] + '\n')
file_dofnames.close()

for iii in rd_dof:
    use_dof = dof[iii]
    for ii in range(int(n_runfiles)):
        rd_c = random.sample(range(length_con), int(n_constraints))
        use_con = []
        for i in rd_c:
            use_con.append(con[i])

        # created dof_ii.txt and constraints_ii.txt and write the 
        # dofs and constr into it
        file_dof = open("runlists/dofs_{}_{}.txt".format(use_dof, ii), "a")
        file_dof.write(use_dof + "\n")
        file_dof.close()
        file_con = open("runlists/constr_{}_{}.txt".format(use_dof, ii), "a")
        for num in range(int(n_constraints)):
            file_con.write(use_con[num] + "\n")
        file_con.close()
