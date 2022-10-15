import sys
from CellModeller.Simulator import Simulator
import os
import numpy as np
import random
import string
import pandas as pd
from SimulationModule import simulation_script
sys.path.append('../SummaryStatistics/')
from DataAnalysis import data_analysis


def calc_sim_summary_statistics(pickle_directory):
    simulation_summary_statistics = data_analysis(pickle_directory, summary_statistic_method_list)

    return simulation_summary_statistics


def simulate(modfilename, platform, device, output_name):
    (path, name) = os.path.split(modfilename)
    modname = str(name).split('.')[0]
    sys.path.append(path)
    sim = Simulator(modname, Sim_dt, outputDirName=output_name, clPlatformNum=platform, clDeviceNum=device,
                    saveOutput=True)
    while len(sim.cellStates) <= sim_max_num_cells_in_last_time_step:
        sim.step()


def run_simulation(param_name, param_val):
    # Execute CellModeller simulation with given parameters
    script_name = ''.join(
        random.SystemRandom().choice(string.ascii_uppercase + string.digits) for _ in range(6))
    print(script_name)
    simulation_script(script_name, {param_name: param_val})

    script_path = "scripts/" + script_name + ".py"
    simulate(script_path, 0, 0, script_name)
    print("Finished Simulation")

    # Start Processing
    print("Calculation of Summary Statistics")

    # Getting the location of simulation pickle files.
    pickle_directory = "data/" + script_name

    return calc_sim_summary_statistics(pickle_directory)


if __name__ == '__main__':

    # create directories
    simulation_results_path = "data"
    scripts_path = "scripts"
    if not os.path.exists(simulation_results_path):
        os.makedirs(simulation_results_path)
    if not os.path.exists(scripts_path):
        os.makedirs(scripts_path)

    global Sim_dt
    global sim_max_num_cells_in_last_time_step
    global summary_statistic_method_list
    Sim_dt = 0.025
    sim_max_num_cells_in_last_time_step = 1000
    summary_statistic_method_list = ['Aspect Ratio', 'Anisotropy']
    aspect_Ratio= []
    anisotropy = []

    print("Running Simulation")

    # parameter distribution
    param_distributions = {'reg_param': [0.01, 1]}
    num_simulation = 100

    for param_name in param_distributions.keys():
        step = (param_distributions[param_name][1] - param_distributions[param_name][0] + 0.01) / num_simulation
        param_values = np.arange(param_distributions[param_name][0], param_distributions[param_name][1], step)
        param_values = np.append(param_values, param_distributions[param_name][1])
        param_values_list = param_values.tolist()
        for param_val in param_values:
            summary_stat = run_simulation(param_name, param_val)
            aspect_Ratio.append(summary_stat[0])
            anisotropy.append(summary_stat[1])

        df = pd.DataFrame(zip(param_values_list, aspect_Ratio, anisotropy), columns=['parameter value', 'Aspect Ratio',
                                                                                     'Anisotropy'])
        df.to_csv(param_name + '.csv')

    print('Finish')
