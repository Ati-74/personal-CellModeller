# Approximate Bayesian Computation with pyabc
# reference: https://buildmedia.readthedocs.org/media/pdf/pyabc/latest/pyabc.pdf
import sys

# note: The addressing needs to be done from the unit test
sys.path.append('../../script/SummaryStatistics/')
from CellModeller.Simulator import Simulator
import os
import multiprocessing
import shutil
import subprocess
import numpy as np
import pyabc
import random
import string
from CreateSimulationScript import simulation_script
from numpy.linalg import norm
from DataAnalysis import data_analysis


def calc_sim_summary_statistics(pickle_directory):

    simulation_summary_statistics = data_analysis(pickle_directory, summary_statistic_method_list, Exp_Sim_dt)

    return simulation_summary_statistics


# distance_function
def distance_calculation(simulation_summary_statistics, exp_summary_statistics):

    distance_list = []
    for summary_statistic in simulation_summary_statistics.keys():
        difference = np.abs(simulation_summary_statistics[summary_statistic] -
                            exp_summary_statistics[summary_statistic])
        distance_list.append(difference)

    distance = norm(distance_list)
    print("#################### distance ####################")
    print(distance)

    return distance


def simulate(modfilename, platform, device, output_name):
    (path, name) = os.path.split(modfilename)
    modname = str(name).split('.')[0]
    sys.path.append(path)
    sim = Simulator(modname, Exp_Sim_dt, outputDirName=output_name, clPlatformNum=platform, clDeviceNum=device,
                    saveOutput=True)
    while len(sim.cellStates) <= sim_max_num_cells_in_last_time_step:
        sim.step()


def abc_simulation(params):

    # Execute CellModeller simulation with given parameters
    script_name = ''.join(
        random.SystemRandom().choice(string.ascii_uppercase + string.digits) for _ in range(6))
    print(script_name)
    simulation_script(script_name, params)

    script_path = "scripts/" + script_name + ".py"
    simulate(script_path, 0, 0, script_name)
    print("Finished Simulation")

    # Start Processing
    print("Calculation of Summary Statistics")

    # Getting the location of simulation pickle files.
    pickle_directory = "data/" + script_name

    return calc_sim_summary_statistics(pickle_directory)


def running_abc(exp_data_directory, param_config, output_directory, n_particles, max_populations, epsilon,
                summary_statistic_method, max_num_cells_in_last_time_step, dt, num_cpu_cores):
    # create directories
    simulation_results_path = "data"
    scripts_path = "scripts"
    if not os.path.exists(simulation_results_path):
        os.makedirs(simulation_results_path)
    if not os.path.exists(scripts_path):
        os.makedirs(scripts_path)

    print("Running ABC")

    # parameter distribution
    prior_distributions = {}
    for parameter_name in param_config.keys():
        param_low = param_config[parameter_name][0]
        param_hi = param_config[parameter_name][1]
        width = abs(param_hi - param_low)

        # generate dictionary containing keywords for pyabc Distribution() class
        prior_distributions[parameter_name] = {"type": "uniform", "args": (param_low, width), "kwargs": {}}

    # create instance of `Distribution` class
    parameter_prior = pyabc.Distribution()
    parameter_prior = parameter_prior.from_dictionary_of_dictionaries(prior_distributions)

    # global variables
    global sim_max_num_cells_in_last_time_step
    global Exp_Sim_dt
    global summary_statistic_method_list

    # interval time
    Exp_Sim_dt = dt

    # calculation of summary statistics for experimental data
    exp_summary_statistics = data_analysis(exp_data_directory, summary_statistic_method, Exp_Sim_dt)

    # simulations configuration
    sim_max_num_cells_in_last_time_step = max_num_cells_in_last_time_step
    summary_statistic_method_list = summary_statistic_method

    # running ABC
    num_cpu_cores_for_abc = pyabc.sampler.MulticoreEvalParallelSampler(num_cpu_cores)
    abc = pyabc.ABCSMC(models=abc_simulation, parameter_priors=parameter_prior, distance_function=distance_calculation,
                       population_size=int(n_particles), sampler=num_cpu_cores_for_abc)

    db_path = ("sqlite:///" + output_directory + "/abc_results.db")

    abc.new(db_path, exp_summary_statistics)

    # Let’s start the sampling now. We’ll sample until the acceptance threshold epsilon drops below 'minimum_epsilon'.
    # We also specify that we want a maximum number of 'max_populations' populations.
    # So whatever is reached first, minimum_epsilon or max_nr_populations, will stop further sampling.
    abc.run(max_nr_populations=max_populations, minimum_epsilon=epsilon)

    print('Finish')
