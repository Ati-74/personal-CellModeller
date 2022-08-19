from lownerJohnEllipse import welzl, plot_ellipse
import numpy as np
import pandas as pd
import scipy.linalg as la
from scipy.spatial import distance_matrix
from sklearn.linear_model import LinearRegression


def fit_enclosing_ellipse(points):
    # convert dataframe to numpy array
    points = points.to_numpy()
    # print(points)
    # finds the smallest ellipse covering a finite set of points
    # https://github.com/dorshaviv/lowner-john-ellipse
    enclosing_ellipse = welzl(points)
    return enclosing_ellipse


def fit_ellipse(bac_in_micro_colony):
    # fit ellipse to micro colony
    # endpoints
    endpoint1_x = bac_in_micro_colony['x_end_point1'].values.tolist()
    endpoint1_y = bac_in_micro_colony['y_end_point1'].values.tolist()
    endpoint2_x = bac_in_micro_colony['x_end_point2'].values.tolist()
    endpoint2_y = bac_in_micro_colony['y_end_point2'].values.tolist()

    endpoints_x = endpoint1_x + endpoint2_x
    endpoints_y = endpoint1_y + endpoint2_y
    endpoints = pd.DataFrame(zip(endpoints_x, endpoints_y))

    # fit ellipse
    ellipse_params = fit_enclosing_ellipse(endpoints)

    return ellipse_params


def aspect_ratio_calc(ellipse_params):
    # calculate aspect ratio
    center_pos, major, minor, theta = ellipse_params
    aspect_ratio = round(minor / major, 3)

    return aspect_ratio


def anisotropy_calc(bac_in_micro_colony):
    # main idea: https://github.com/ingallslab/bsim-related/blob/main/bsim_related/data_processing/cell_data_processing.py#L184

    local_anisotropies = []

    # orientation of bacteria
    bacteria_orientation = bac_in_micro_colony["orientation"]

    for bacterium_index in range(bac_in_micro_colony.shape[0]):
        # Projection matrix
        projection_matrix = np.zeros(shape=(2, 2))
        for neighbor_index in range(bac_in_micro_colony.shape[0]):
            if neighbor_index != bacterium_index:
                # Compute the sum of the projection matrices on the orientation vectors of the neighbouring bacteria
                # projection matrix
                """
                cos(angle)                  cos(angle)*sin(angle)
                cos(angle)*sin(angle)       sin(angle)
                """
                projection_matrix += np.matrix([[np.cos(bacteria_orientation.iloc[neighbor_index]) ** 2,
                                                 np.cos(bacteria_orientation.iloc[neighbor_index]) * np.sin(
                                                     bacteria_orientation.iloc[neighbor_index])],
                                                [np.cos(bacteria_orientation.iloc[neighbor_index]) * np.sin(
                                                    bacteria_orientation.iloc[neighbor_index]),
                                                 np.sin(bacteria_orientation.iloc[neighbor_index]) ** 2]])

        # Compute the mean of the projection matrices on the orientation vectors of the neighbouring bacteria
        num_neighbours = bac_in_micro_colony.shape[0] - 1
        projection_matrix = projection_matrix / num_neighbours
        # Get the max real eigenvalues of the mean projection matrix; this is the local anisotropy
        local_anisotropies.append(max(la.eigvals(projection_matrix).real))

    # calculate mean anisotropy
    mean_anisotropy = round(np.mean(local_anisotropies), 3)

    return mean_anisotropy


