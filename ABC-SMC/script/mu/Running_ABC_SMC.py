# Approximate Bayesian Computation with pyabc
import sys

# note: The addressing needs to be done from the unit test
sys.path.append('../../script/mu/summary statistics/simulation')
sys.path.append('../../script/mu/summary statistics/experiment')
from CellModeller.Simulator import Simulator
import os
import multiprocessing
import shutil
import subprocess
import numpy as np
import pyabc
import pandas
import random
import string
from CreateSimulationScript import simulation_script
from numpy.linalg import norm
from simulation_summary_statistics_analysis import process_sim_micro_colonies
from exp_summary_statistics_analysis import process_exp_micro_colonies
import tempfile
import uuid


def get_features(pickle_directory, summary_statistic_method_list):
    # cell_modeller_features list
    cell_modeller_features = process_sim_micro_colonies(pickle_directory, summary_statistic_method_list)

    return cell_modeller_features


# distance_function
def distance_calculation(cell_modeller_data_summary_statistics, exp_summary_statistics):
    print(cell_modeller_data_summary_statistics)
    print(exp_summary_statistics)
    n_features = len(cell_modeller_data_summary_statistics.keys())
    distance_vector = np.ones((n_features,))
    i = 0
    for feature in cell_modeller_data_summary_statistics.keys():
        diff = np.abs(np.mean(cell_modeller_data_summary_statistics[feature]) -
                      np.mean(exp_summary_statistics[feature]))
        distance_vector[i] = diff
        i += 1

    distance = norm(distance_vector)

    return distance


def simulate(modfilename, platform, device, output_name, steps=50):
    # total_simulation_time = 608
    last_time_step_bacteria_num = 409
    dt = 0.08

    (path, name) = os.path.split(modfilename)
    modname = str(name).split('.')[0]
    sys.path.append(path)
    sim = Simulator(modname, dt, outputDirName=output_name, clPlatformNum=platform, clDeviceNum=device,
                    saveOutput=True)
    while len(sim.cellStates) <= last_time_step_bacteria_num:
        sim.step()


def abc_simulation(params):
    # return loaded csv file of sim data
    # Execute CellModeller simulation with given parameters
    script_name = ''.join(
        random.SystemRandom().choice(string.ascii_uppercase + string.digits) for _ in range(6))
    print(script_name)

    parameter = ','.join([str(params[k]) for k in params])
    simulation_script(script_name, str(parameter))

    script_path = "scripts/" + script_name + ".py"
    # sys.stdout = open(os.devnull, 'w')
    simulate(script_path, 0, 0, script_name)
    print("Finished Simulation")

    # converting CellModeller Output Structure to CP Output Structure
    # Getting the location of simulations.
    pickle_directory = "data/" + script_name
    summary_statistic_method_list = ["Aspect Ratio", "Anisotropy"]

    # Start Processing
    print("start post-processing")
    return get_features(pickle_directory, summary_statistic_method_list)


def running_abc(exp_data, relationship_file, param_config, output_directory, n_particles, max_populations, epsilon):
    print(epsilon)
    # create directories
    simulation_result_path = "data"
    scripts_path = "scripts"
    if not os.path.exists(simulation_result_path):
        os.makedirs(simulation_result_path)
    if not os.path.exists(scripts_path):
        os.makedirs(scripts_path)

    print("Running ABC test")
    # parameter set for AstroABC package
    prior_distributions = {}

    for parameter_info in param_config:
        # split line in parameter file to get name, low, hi
        split_line = parameter_info.split(",")
        param_name = split_line[0]
        param_low = float(split_line[1])
        param_hi = float(split_line[2])
        width = abs(param_hi - param_low)

        # generate dictionary containing keywords for pyabc Distribution() class
        prior_distributions[param_name] = {"type": "uniform", "kwargs": {"loc": param_low + 1e-2, "scale": width}}

    parameter_prior = pyabc.Distribution()
    parameter_prior = parameter_prior.from_dictionary_of_dictionaries(prior_distributions)

    # extract features from experimental data
    summary_statistic_method_list = ["Aspect Ratio", "Anisotropy"]
    exp_features = process_exp_micro_colonies(exp_data, relationship_file, summary_statistic_method_list)

    # running ABC
    abc = pyabc.ABCSMC(models=abc_simulation, parameter_priors=parameter_prior, distance_function=distance_calculation,
                       population_size=int(n_particles))

    db_path = ("sqlite:///" + output_directory + "/spring_constants_run.db")

    abc.new(db_path, exp_features)

    # Let’s start the sampling now. We’ll sample until the acceptance threshold epsilon drops below 'minimum_epsilon'.
    # We also specify that we want a maximum number of 'max_populations' populations.
    # So whatever is reached first, minimum_epsilon or max_nr_populations, will stop further sampling.
    abc.run(max_nr_populations=int(max_populations), minimum_epsilon=epsilon)

    print('Finish')
