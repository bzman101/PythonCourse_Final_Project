## Libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import re

## Functions

def df_cleanup(df,done_by):
    """
    This Function receives a df of a freq file returns the same df with 4 columns added:
    MOI - presenting the MOI of the experiment
    Passage - indicating the passage of each mutation
    Line - Which line of initial bacterial population
    Done_By - indicating the scientist responsible for the experiment
    return: df_arranged
    """
    df_arranged = df.copy()
    # Extract the number after the first "p" to a new column 'Passage'
    df_arranged['Passage'] = df_arranged['File'].str.extract(r'p(\d+)').astype(int)
    # Extract the letter after the number to a new column 'Line'
    df_arranged['Line'] = df_arranged['File'].str.extract(r'p\d+([A-Za-z])')
    # Extract the number after the string "moi" to a new column 'MOI'
    df_arranged['MOI'] = df_arranged['File'].str.extract(r'moi(\d+)').astype(int)
    df_arranged['Done_by']= done_by
    return df_arranged
def get_mut_list(merged_df):
    """
    This Function receives a merged df of freq files returns a new df containing the following:
    Position - the position of the nucleotide in query
    Ref_nuc - The nucleotide present in the reference sequence
    Mut_nuc - The nucleotide present in the mutative sequence
    Mutation - the mutation written in a formatted way :(N_ref+position+N_mut)
    Coverage - How many reads where mapped to that mutation
    Avg_Score, Frequency, MOI, Passage, Done_By
    The new df will be cleaned up from all Null
    return: df_arranged
    """

    #return mut_df

def mut_cutoffs(mut_df,min_coverage,min_freq):
    """
    This Function receives a mutants df, a minimum coverage value and a minimum frequency value.
    It returns the same df where the mutants that does not meet the criteria are dropped.
    """

    #return filtered_mut_df

## Constants and Parameters
#maria_path = "/Users/adibnzv/Desktop/DESKTOP/SCHOOL/PhD/Year 1/Semester 2/Python for Biologists/Final Project/PythonCourse_Final_Project/DATA/Maria/Maria1.csv"
carmel_path = "/Users/adibnzv/Desktop/DESKTOP/SCHOOL/PhD/Year 1/Semester 2/Python for Biologists/Final Project/PythonCourse_Final_Project/DATA/Carmel/Carmel1.csv"
shir_path = "/Users/adibnzv/Desktop/DESKTOP/SCHOOL/PhD/Year 1/Semester 2/Python for Biologists/Final Project/PythonCourse_Final_Project/DATA/Shir/Shir1.csv"

## Main Code

#Load the intial freq.csv files into a pandas df
#df_maria =pd.read_csv(maria_path)
df_carmel =pd.read_csv(carmel_path)
df_shir =pd.read_csv(shir_path)

#Fix Shir's 'File' column to match Carmel's
# Remove the hyphen after the number
df_shir['File'] = df_shir['File'].str.replace('-', '', n=1)
# Insert "_moi10_" after the capitalized letter
df_shir['File'] = df_shir['File'].str.replace('([A-Z])', r'\1_moi10_', n=1)

#Clean up the freq files and add relevant data
df_arr_carmel = df_cleanup(df_carmel, "Carmel")
df_arr_shir = df_cleanup(df_shir, "Shir")

# Create a Joined Freq File and remove lines where base_count is zero
joined_freq = pd.concat([df_arr_carmel, df_arr_shir])
joined_freq = joined_freq[joined_freq['base_count'] != 0]

