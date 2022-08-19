from exp_MicroColonyAnalysis import micro_colony_analysis as exp_colony_analysis


def process_exp_micro_colonies(input_file, relationship_file, summary_statistic_method_list):

    # start analysis
    summary_statistic = exp_colony_analysis(input_file, relationship_file, summary_statistic_method_list)

    return summary_statistic

