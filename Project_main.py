## Libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as ticker
import sys
import re
from MS2_COLORS import COLORS

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
    This Function receives a merged df of freq files adds a Mutation field:
    Mutation - the mutation written in a formatted way :(N_ref+position+N_mut)
    In addition, the function will lose mutation from known problematic regions or primer region.
    return: mut_df
    """
    # Add Mutation column
    mut_df = merged_df.copy()
    mut_df['Full Mutation'] = mut_df['ref_base'] + mut_df['ref_pos'].astype(str) + mut_df['read_base']
    # Lose problematic lines (Primers and Region of known sequencing error)
    primers = list(range(1, 20)) + list(range(1291, 1304)) + list(range(1179, 1200)) + list(range(2270, 2288)) + \
              list(range(2167, 2188)) + list(range(3548, 3570))
    problematic = [17, 18, 19, 20, 21, 22, 23, 183, 188, 189, 190, 196, 274, 317, 364, 452, 762, 2719, 3117, 3133, 3139,
                   3143, 3146, 3150, 3401, 3539, 3542]
    remove_list = primers + problematic
    mut_df = mut_df[~joined_freq['ref_pos'].isin(remove_list)]
    # Remove Non Mutations
    mut_df = mut_df[mut_df['read_base'] != mut_df['ref_base']]
    return mut_df

def mut_cutoffs(mut_df,min_coverage,min_frequency):
    """
    This Function receives a mutants df, a minimum coverage value and a minimum frequency value.
    It returns the same df where the mutants that does not meet the criteria are dropped.
    """
    filtered_mut_df = mut_df.copy()
    filtered_mut_df = filtered_mut_df[filtered_mut_df['base_count'] >= min_coverage]
    filtered_mut_df = filtered_mut_df[filtered_mut_df['frequency'] >= min_frequency]

    return filtered_mut_df

def create_mutation_figure(df, mutations_list, output_path):
    """
    This function gets a df of freq files, a list of mutations and a path to save graph to.
    """
    plt.style.use('ggplot')
    # Filter the mutations_list to only include mutations that appear in more than one passage
    mutations_list = [m for m in mutations_list if df[df['Full Mutation'] == m]['Passage'].nunique() > 1]
    # Create a subplot for each mutation with a smaller size
    fig, axes = plt.subplots(nrows=len(mutations_list), ncols=1, figsize=(5, 2.5 * len(mutations_list)))
    # If there's only one mutation, axes will not be an array, so we convert it to an array
    if len(mutations_list) == 1:
        axes = [axes]
    for mutation, ax in zip(mutations_list, axes):
        # Filter the DataFrame for rows where 'Full Mutation' is equal to the mutation and create a copy
        df_mutation = df[df['Full Mutation'] == mutation].copy()
        # Create a new column that concatenates 'Line', 'MOI', and 'Done_by'
        df_mutation['Experiment'] = df_mutation['Line'] + ' - ' + df_mutation['MOI'].astype(str) + ' - ' + df_mutation['Done_by']
        # Create a line plot for each unique combination of 'Line', 'MOI', and 'Done_by'
        sns.lineplot(x='Passage', y='frequency', hue='Experiment', data=df_mutation, ax=ax)
        # Set the title of the subplot to the mutation
        ax.set_title(mutation)
        # Set the y-limit to [0, 1]
        ax.set_ylim(0, 1)
        # Set x-axis to only contain integers
        ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
        # Set legend font size
        ax.legend(fontsize='small')
    # Adjust the layout
    plt.tight_layout()
    # Save the figure with a lower dpi
    plt.savefig(output_path, dpi=300)
    # Show the figure
    plt.show()


## Constants and Parameters
maria_path = "/Users/adibnzv/Desktop/DESKTOP/SCHOOL/PhD/Year 1/Semester 2/Python for Biologists/Final Project/PythonCourse_Final_Project/DATA/Maria/Maria1.csv"
carmel_path = "/Users/adibnzv/Desktop/DESKTOP/SCHOOL/PhD/Year 1/Semester 2/Python for Biologists/Final Project/PythonCourse_Final_Project/DATA/Carmel/Carmel1.csv"
shir_path = "/Users/adibnzv/Desktop/DESKTOP/SCHOOL/PhD/Year 1/Semester 2/Python for Biologists/Final Project/PythonCourse_Final_Project/DATA/Shir/Shir1.csv"
export_path = "/Users/adibnzv/Desktop/DESKTOP/SCHOOL/PhD/Year 1/Semester 2/Python for Biologists/Final Project/PythonCourse_Final_Project/Export/"
min_freq = 0.01
min_cov = 100


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

# Create a Mutation DF File
Mutation_df = get_mut_list(joined_freq)

# Fileter out using the cutoffs (can be found in the parameters section)
Mutation_df_filtered = mut_cutoffs(Mutation_df,min_cov,min_freq)

#Get a list of Muataion that passed the filteration
mut_lst = Mutation_df_filtered['Full Mutation'].tolist()
mut_lst = list(set(mut_lst))

# Create Graph per Mutation and save them to the Export folder:
create_mutation_figure(Mutation_df_filtered, mut_lst, export_path+'Figure1')

