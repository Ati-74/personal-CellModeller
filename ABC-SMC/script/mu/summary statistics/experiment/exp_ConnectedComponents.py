
"""
@author:  Atiyeh Ahmadi - Aaron Yip - Siddh
"""
import pandas as pd


def edge_list_to_adjacency_matrix(col1, col2):

    """
    Goal: neighbor cells' ID (of Object_relationship.csv) will be used to create adjacency data frame
    
    @param col1   data frame   First Object Number from Object_relationship.csv
    @param col2   data frame   Second Object Number from Object_relationship.csv
    """    
    df = pd.crosstab(col1, col2)
    idx = df.columns.union(df.index)
    df = df.reindex(index=idx, columns=idx, fill_value=0)
    return df
