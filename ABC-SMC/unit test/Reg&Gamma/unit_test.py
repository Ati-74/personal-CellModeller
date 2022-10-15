import os
import sys
sys.path.append('../../script/')
from SampleScript2module import create_simulation_module_file

if __name__ == '__main__':

    # first step: create simulation module file
    # sample simulation script file
    sample_simulation_file = "../../input/SampleSimulation.py"
    params = ["gamma", "reg_param"]
    # Parameters to be calibrated pattern in script file
    params_for_calibration = ["gamma=10", "reg_param=0.1"]
    # By providing an example of the desired simulation file that includes the parameters to calibrate,
    #          the user can generate the module needed to run the simulation in ABC.
    create_simulation_module_file(sample_simulation_file, params, params_for_calibration)

    # second step: Running ABC
    from Running_ABC_SMC import running_abc
    # parameter name: [low value, high value]
    param_config = {'reg_param': [0, 1], 'gamma': [0, 100]}

    # dataset used to compare ABC simulations against
    exp_data_directory = "../../input/K12 Experimental data"
    # save abc db
    output_directory = os.getcwd()

    # configure ABC specific arguments
    number_of_particles = 50
    max_populations = 10
    epsilon = 0.05
    # "Aspect Ratio", "Anisotropy"
    summary_statistic_method_list = ["Aspect Ratio", "Anisotropy"]
    # Simulation termination condition
    max_num_cells_in_last_time_step = 409
    dt = 0.025

    num_cpu_cores = 6

    running_abc(exp_data_directory, param_config, output_directory, number_of_particles, max_populations, epsilon,
                summary_statistic_method_list, max_num_cells_in_last_time_step, dt, num_cpu_cores)
