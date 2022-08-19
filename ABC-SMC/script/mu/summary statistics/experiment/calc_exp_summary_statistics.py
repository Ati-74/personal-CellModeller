import numpy as np
import pandas as pd
from lownerJohnEllipse import welzl
import scipy.linalg as la


def fit_enclosing_ellipse(points):
    # convert dataframe to numpy array
    points = points.to_numpy()
    # print(points)
    # finds the smallest ellipse covering a finite set of points
    # https://github.com/dorshaviv/lowner-john-ellipse
    enclosing_ellipse = welzl(points)
    return enclosing_ellipse


def bac_info(bac_in_micro_colony):
    """
    Goal: this function returns important features of bacteria in the micro-colony that we
    study (ex: orientation, minor/major axis length, etc.). we need this data for calculating endpoints for
    fitting ellipse and plotting bacteria.

    @param bac_in_micro_colony     data frame   information(orientation, major/minor axis length, etc.)
    about each bacterium in each micro-colony in current time-step.

    """
    # center coordinate
    bacteria_center_coord = bac_in_micro_colony[["AreaShape_Center_X", "AreaShape_Center_Y"]]
    # major axis length
    bacteria_major_axis = bac_in_micro_colony["AreaShape_MajorAxisLength"]
    # minor axis length
    bacteria_minor_axis = bac_in_micro_colony["AreaShape_MinorAxisLength"]
    # orientation
    bacteria_orientation = bac_in_micro_colony["AreaShape_Orientation"]

    return bacteria_center_coord, bacteria_major_axis, bacteria_minor_axis, bacteria_orientation


def fit_ellipse(bac_in_micro_colony):
    # bacteria info
    (bacteria_center_coord, bacteria_major, bacteria_minor, bacteria_orientation) = bac_info(bac_in_micro_colony)
    bac_angle = -(bacteria_orientation + 90) * np.pi / 180

    # endpoints
    endpoint1_x = bacteria_center_coord["AreaShape_Center_X"] + (bacteria_major / 2) * np.cos(bac_angle)
    endpoint1_y = bacteria_center_coord["AreaShape_Center_Y"] + (bacteria_major / 2) * np.sin(bac_angle)
    endpoint2_x = bacteria_center_coord["AreaShape_Center_X"] - (bacteria_major / 2) * np.cos(bac_angle)
    endpoint2_y = bacteria_center_coord["AreaShape_Center_Y"] - (bacteria_major / 2) * np.sin(bac_angle)
    endpoint1 = pd.concat([endpoint1_x, endpoint1_y], axis=1)
    endpoint2 = pd.concat([endpoint2_x, endpoint2_y], axis=1)
    endpoints = pd.concat([endpoint1, endpoint2], axis=0)
    # fit ellipse
    ellipse_params = fit_enclosing_ellipse(endpoints)

    return ellipse_params


def aspect_ratio_calc(ellipse_params):
    # calculate aspect ratio
    center_pos, major, minor, theta = ellipse_params
    aspect_ratio = round(minor / major, 3)

    return aspect_ratio


def anisotropy_calc(bac_in_micro_colony):
    # bacteria info
    (bacteria_center_coord, bacteria_major, bacteria_minor, bacteria_orientation) = bac_info(bac_in_micro_colony)

    # main idea: https://github.com/ingallslab/bsim-related/blob/main/bsim_related/data_processing/cell_data_processing.py#L184

    local_anisotropies = []

    # modification of orientation
    bac_angle = -(bacteria_orientation + 90) * np.pi / 180

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
                projection_matrix += np.matrix([[np.cos(bac_angle.iloc[neighbor_index]) ** 2,
                                                 np.cos(bac_angle.iloc[neighbor_index]) *
                                                 np.sin(bac_angle.iloc[neighbor_index])],
                                                [np.cos(bac_angle.iloc[neighbor_index]) *
                                                 np.sin(bac_angle.iloc[neighbor_index]),
                                                 np.sin(bac_angle.iloc[neighbor_index]) ** 2]])

        # Compute the mean of the projection matrices on the orientation vectors of the neighbouring bacteria
        num_neighbours = bac_in_micro_colony.shape[0] - 1
        projection_matrix = projection_matrix / num_neighbours
        # Get the max real eigenvalues of the mean projection matrix; this is the local anisotropy
        local_anisotropies.append(max(la.eigvals(projection_matrix).real))

    # calculate mean anisotropy
    mean_anisotropy = round(np.mean(local_anisotropies), 3)

    return mean_anisotropy
