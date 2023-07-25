#Libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys

##Functions
def df_cleanup(df):
    """
    This Function receives a df of a freq file returns the same df with 3 columns extra:
    MOI - presenting the MOI of the experiment
    Passage - indicating the passage of each mutation
    Done_By - indicating the scientist responsible for the experiment
    return: df_arranged
    """

    return df_arranged
def get_mut_list(merged_df):
    """
    This Function receives a merged df of freq files returns a new df containing the following:
    Position - the position of the nucleotide in query
    Ref_nuc - The nucleotide present in the reference sequence
    Mut_nuc - The nucleotide present in the mutative sequence
    Mutation - the mutation written in a formatted way :(N_ref+position+N_mut)
    Coverage - How many reads where mapped to that mutation
    Frequency - The frequency of the position.
    MOI - presenting the MOI of the experiment
    Passage - indicating the passage of each mutation
    Done_By - indicating the scientist responsible for the experiment
    return: df_arranged
    """

    return mut_df

def mut_cutoffs(mut_df,min_coverage,min_freq):
    """
    This Function receives a mutants df, a minimum coverage value and a minimum frequency value.
    It returns the same df where the mutants that does not meet the criteria are dropped.
    """

    return filtered_mut_df