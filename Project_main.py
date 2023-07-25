import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys

##Functions

def df_cleanup(df):
    """
    :param:df
    This Function receives a df of a freq file returns the same df with 3 columns extra:
    MOI - presenting the MOI of the experiment
    Passage - indicating the passage of each mutation
    Done_By - indicating the scientist responsible for the experiment
    :return: df_arranged
    """

    return df_arranged