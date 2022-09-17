import os
import sys
sys.path.append('../../script/')
from Running_ABC_SMC import running_abc


# parameter name: [low value, high value]
param_config = {'mu': [0, 2]}

# dataset used to compare ABC simulations against
exp_data_directory = "../../input/K12 Experimental data"
# save abc db
output_directory = os.getcwd()

# configure ABC specific arguments
number_of_particles = 4
max_populations = 2
epsilon = 0.05
summary_statistic_method_list = ["Aspect Ratio", "Anisotropy", 'Density', 'dist_vs_growth_rate']
# Simulation termination condition
max_num_cells_in_last_time_step = 409
dt = 0.025

num_cpu_cores = 8

running_abc(exp_data_directory, param_config, output_directory, number_of_particles, max_populations, epsilon,
            summary_statistic_method_list, max_num_cells_in_last_time_step, dt, num_cpu_cores)
