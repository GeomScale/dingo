#!/usr/bin/env python3
import os
import json
import numpy as np
from dingo.loading_models import read_json_file, read_mat_file


current_directory = os.getcwd()
input_file_json = current_directory +  '/ext_data/e_coli_core.json'
input_file_mat = current_directory +  '/ext_data/e_coli_core.mat'


met_net = read_json_file(input_file_json)

print("lb")
print(met_net[0])
print("------------------------")
print("ub")
print(met_net[1])
print("-------------------------")
print("S")
print(met_net[2])
print("-------------------------")
print("Metabolites are the following " + str(len(met_net[3])) + ":")
print(met_net[3])
print("-------------------------")
print("Reactions are the following " + str(len(met_net[4])) + ":")
print(met_net[4])

# -----------------------------------------------------------------

met_net = read_mat_file(input_file_mat)

print("lb")
print(met_net[0])
print("------------------------")
print("ub")
print(met_net[1])
print("-------------------------")
print("S")
print(met_net[2])
print("-------------------------")
print("Metabolites are the following " + str(len(met_net[3])) + ":")
print(met_net[3])
print("-------------------------")
print("Reactions are the following " + str(len(met_net[4])) + ":")
print(met_net[4])

