import os
import sys
import pandas as pd
sys.path.append('../../script/mu')
from Running_ABC_SMC import running_abc

# Parameter list with containing lines with format
# '<parameter name>,<low value>,<high value>
param_config = '../../input/param_config_mu.txt'
# dataset used to compare ABC simulations against
exp_data = "../../input/K12_scene2.csv"
relationship_file = "../../input/K12_scene2_Object relationships.csv"

output_directory = os.getcwd()
print(output_directory)

# configure ABC specific arguments
number_of_particles = 4
max_populations = 2
epsilon = 0.05

# read txt file
param_config = open(param_config).read().splitlines()

# read experimental data
# exp_data = pd.read_csv(exp_data)

running_abc(exp_data, relationship_file, param_config, output_directory, number_of_particles, max_populations, epsilon)
