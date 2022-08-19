from sim_MicroColonyAnalysis import micro_colony_analysis as sim_colony_analysis


def process_sim_micro_colonies(input_directory, summary_statistic_method_list):
    # start analysis
    summary_statistic = sim_colony_analysis(input_directory, summary_statistic_method_list)

    return summary_statistic
