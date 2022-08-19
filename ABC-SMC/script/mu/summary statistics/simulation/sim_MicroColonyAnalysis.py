import sys
import glob
import networkx as nx
import pandas as pd
import numpy as np
import pickle
from scipy.spatial import distance_matrix
from scipy.stats import lognorm
from calc_simulation_summary_statistics import fit_ellipse, aspect_ratio_calc, anisotropy_calc
from sim_ConnectedComponents import edge_list_to_adjacency_matrix

# set number of recursion
# https://stackoverflow.com/questions/14222416/recursion-in-python-runtimeerror-maximum-recursion-depth-exceeded-while-callin
sys.setrecursionlimit(10000)


def finding_sub_graphs(graph):
    """
    Goal: Finding the largest disconnected subparagraphs
    @param graph   graph  neighboring bacteria graph in current time-step
    """
    sub_graphs = sorted(nx.connected_components(graph), key=len, reverse=True)

    return sub_graphs


def make_adjacency_matrix(cs):
    """
    goal: make adjacency matrix
    @param cs dictionary  bacteria information
    @return: adjacency_matrix_df   dataframe  adjacency matrix dataframe
    """
    first_vertex_of_link = []
    second_vertex_of_link = []
    for it in cs:
        neighbours_list = cs[it].neighbours
        for neighbour_id in neighbours_list:
            first_vertex_of_link.append(cs[it].id)
            second_vertex_of_link.append(neighbour_id)
    adjacency_matrix = edge_list_to_adjacency_matrix(first_vertex_of_link, second_vertex_of_link)

    return adjacency_matrix


def make_distance_matrix(micro_colony_1, micro_colony_2, min_distance_between_micro_colonies):
    # each row shows distance of a bacterium in target micro colony to another micro colony bacteria
    distance_df = pd.DataFrame(distance_matrix(micro_colony_1, micro_colony_2))

    # find neighbour bacteria of two micro colonies
    neighbor_bacteria = distance_df.loc[(distance_df <= min_distance_between_micro_colonies).any(axis=1)]

    return neighbor_bacteria


def micro_colony_analysis(pickle_files_directory, summary_statistic_method):
    """
    Goal: this is the main function; the function calls other function that are described above and
    calculates aspect ratio for each micro colony in each time-step and store the outputs in png format and pickle.

    @param summary_statistic_method   str     the method that we apply on the micro-colonies
    
    """

    min_size_of_micro_colony = 2
    min_distance_between_micro_colonies = 2

    # store micro colonies summary statistics
    summary_statistic = {}
    # store aspect ratio of micro colonies
    aspect_ratio_list = []
    anisotropy_list = []

    # read pickle files
    path = pickle_files_directory + "/*.pickle"
    filename_list = [filename for filename in sorted(glob.glob(path))]

    for cnt, filename in enumerate(filename_list):

        # store ellipses
        ellipses = []
        # store anisotropies of micro colonies of this rime step
        local_anisotropy_list = []
        # store anisotropies of micro colonies of this rime step
        local_aspect_ratio_list = []
        micro_colonies_in_current_time_step = []

        timestep = cnt + 1
        print('time step:' + str(timestep))

        # read current pickle file
        current_bacteria_info = pickle.load(open(filename_list[cnt], 'rb'))
        cs = current_bacteria_info['cellStates']

        # find neighbours
        adjacency_matrix = make_adjacency_matrix(cs)

        # create graph
        graph = nx.from_pandas_adjacency(adjacency_matrix)

        # finding sub-graphs
        sub_graphs = finding_sub_graphs(graph)

        # get important features of bacteria to draw them
        bacteria_minor_in_current_timestep = [cs[it].radius for it in cs]
        # end points of bacteria
        bacteria_first_end_points_x_in_current_timestep = [cs[it].ends[0][0] for it in cs]
        bacteria_first_end_points_y_in_current_timestep = [cs[it].ends[0][1] for it in cs]
        bacteria_second_end_points_x_in_current_timestep = [cs[it].ends[1][0] for it in cs]
        bacteria_second_end_points_y_in_current_timestep = [cs[it].ends[1][1] for it in cs]

        # convert to dataframe
        bacteria_info = list(zip(bacteria_minor_in_current_timestep,
                                 bacteria_first_end_points_x_in_current_timestep,
                                 bacteria_first_end_points_y_in_current_timestep,
                                 bacteria_second_end_points_x_in_current_timestep,
                                 bacteria_second_end_points_y_in_current_timestep))
        columns_name = ['minor', 'x_end_point1', 'y_end_point1', 'x_end_point2', 'y_end_point2']
        df_current_time_step = pd.DataFrame(bacteria_info, columns=columns_name)

        for subgraph in sub_graphs:
            # check micro-colony
            micro_colony_is_merged = False
            if len(subgraph) >= min_size_of_micro_colony:
                subgraph_nodes = list(subgraph)
                bacteria_x_center_in_this_micro_colony = [cs[it].pos[0] for it in cs if cs[it].id in subgraph_nodes]
                bacteria_y_center_in_this_micro_colony = [cs[it].pos[1] for it in cs if cs[it].id in subgraph_nodes]
                for other_subgraph in sub_graphs:
                    if other_subgraph != subgraph and len(other_subgraph) >= min_size_of_micro_colony:
                        other_subgraph_nodes = list(other_subgraph)
                        bacteria_x_center_in_other_micro_colony = [cs[it].pos[0] for it in cs if
                                                                   cs[it].id in other_subgraph_nodes]
                        bacteria_y_center_in_other_micro_colony = [cs[it].pos[1] for it in cs if
                                                                   cs[it].id in other_subgraph_nodes]
                        # find neighbour bacteria of two micro colonies
                        neighbor_bacteria = make_distance_matrix(list(zip(bacteria_x_center_in_this_micro_colony,
                                                                          bacteria_y_center_in_this_micro_colony)),
                                                                 list(zip(bacteria_x_center_in_other_micro_colony,
                                                                          bacteria_y_center_in_other_micro_colony)),
                                                                 min_distance_between_micro_colonies)
                        if neighbor_bacteria.shape[0] >= 1:
                            # it means that two micro colonies are close together
                            micro_colony_is_merged = True
                            break
                        else:
                            pass
                if not micro_colony_is_merged:
                    # store micro colony bacteria id
                    micro_colonies_in_current_time_step.append(subgraph_nodes)

                    # get bacteria information in this micro-colony
                    bacteria_info_in_this_micro_colony = [cs[it] for it in cs if cs[it].id in subgraph_nodes]
                    # end points of bacteria
                    bacteria_first_end_points_x_in_this_micro_colony = \
                        [cs[it].ends[0][0] for it in cs if cs[it].id in subgraph]
                    bacteria_first_end_points_y_in_this_micro_colony = \
                        [cs[it].ends[0][1] for it in cs if cs[it].id in subgraph]
                    bacteria_second_end_points_x_in_this_micro_colony = \
                        [cs[it].ends[1][0] for it in cs if cs[it].id in subgraph]
                    bacteria_second_end_points_y_in_this_micro_colony = \
                        [cs[it].ends[1][1] for it in cs if cs[it].id in subgraph]
                    bacteria_orientation_in_this_micro_colony = \
                        [np.arctan2(cs[it].dir[1], cs[it].dir[0]) for it in cs if cs[it].id in subgraph]

                    bacteria_in_this_micro_info = list(zip(bacteria_first_end_points_x_in_this_micro_colony,
                                                           bacteria_first_end_points_y_in_this_micro_colony,
                                                           bacteria_second_end_points_x_in_this_micro_colony,
                                                           bacteria_second_end_points_y_in_this_micro_colony,
                                                           bacteria_orientation_in_this_micro_colony))
                    columns_name = ['x_end_point1', 'y_end_point1', 'x_end_point2', 'y_end_point2', 'orientation']
                    # create dataframe
                    bacteria_in_this_micro_colony_df = \
                        pd.DataFrame(bacteria_in_this_micro_info, columns=columns_name)
                    # print(bacteria_in_this_micro_colony_df)

                    # fit ellipse
                    ellipse_params = fit_ellipse(bacteria_in_this_micro_colony_df)
                    # append ellipse
                    ellipses.append(ellipse_params)

                    """
                                calculation of aspect ratio
                    """
                    if "Aspect Ratio" in summary_statistic_method:
                        aspect_ratio = aspect_ratio_calc(ellipse_params)
                        # store aspect ratio
                        aspect_ratio_list.append(aspect_ratio)
                        local_aspect_ratio_list.append(aspect_ratio)
                    """
                                Finish
                    """

                    """
                        calculation of Anisotropy
                    """
                    if "Anisotropy" in summary_statistic_method:
                        mean_anisotropy = anisotropy_calc(bacteria_in_this_micro_colony_df)
                        # store anisotropy
                        anisotropy_list.append(mean_anisotropy)
                        local_anisotropy_list.append(mean_anisotropy)
                    """
                        Finish
                    """

    if aspect_ratio_list:
        shape, loc, scale = lognorm.fit(aspect_ratio_list)
        summary_statistic['aspect ratio'] = round(np.log(scale), 2)
    if anisotropy_list:
        shape, loc, scale = lognorm.fit(anisotropy_list)
        summary_statistic['anisotropy'] = round(np.log(scale), 2)

    return summary_statistic
