## Libraries
import pandas as pd
import random
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as ticker
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
from MS2_COLORS import COLORS
import math
from matplotlib.lines import Line2D
import numpy as np

## Functions
def df_cleanup(df, done_by):
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
    df_arranged['Done_by'] = done_by
    return df_arranged
def get_mut_column(merged_df):
    """
    This Function receives a merged df of freq files adds a Mutation field:
    Mutation - the mutation written in a formatted way :(N_ref + position + N_mut)
    return: mut_df
    """
    # Add Mutation column
    mut_df = merged_df.copy()
    mut_df['Full Mutation'] = mut_df['ref_base'] + mut_df['ref_pos'].astype(str) + mut_df['read_base']
    # Remove Non Mutations
    mut_df = mut_df[mut_df['read_base'] != mut_df['ref_base']]
    return mut_df
def mut_cutoffs(mut_df,min_coverage,min_frequency):
    """
    This Function receives a mutants df, a minimum coverage value and a minimum frequency value.
    It returns the df when it is filtered by these cutoffs
    In addition, the function will lose mutation from known problematic regions or primer region.
    """
    filtered_mut_df = mut_df.copy()
    # Lose problematic lines (Primers and Region of known sequencing error)
    primers = list(range(1, 20)) + list(range(1291, 1304)) + list(range(1179, 1200)) + list(range(2270, 2288)) + \
              list(range(2167, 2188)) + list(range(3548, 3570))
    problematic = [17, 18, 19, 20, 21, 22, 23, 183, 188, 189, 190, 196, 274, 317, 364, 452, 762, 2719, 3117, 3133, 3139,
                   3143, 3146, 3150, 3401, 3539, 3542, 3547]
    remove_list = primers + problematic
    filtered_mut_df = filtered_mut_df[~filtered_mut_df['ref_pos'].isin(remove_list)]
    # Lose only the mutations that has not met the cutoffs in any of the passages.
    filtered_mut_df = filtered_mut_df[filtered_mut_df['base_count'] >= min_coverage]
    filtered_mut_df = filtered_mut_df[filtered_mut_df['frequency'] >= min_frequency]
    relevant_mutations = filtered_mut_df['Full Mutation'].tolist()
    relevant_mutations = list(set(relevant_mutations))
    relevant_mutations.sort()
    return relevant_mutations
def get_color_for_new_mutation(mutation):
    possible_colors = list(mcolors.CSS4_COLORS.keys())
    relevant_colors = [x for x in possible_colors if x not in COLORS.items()]
    picked_color = random.choice(relevant_colors)
    relevant_mut_colors[mutation]=picked_color
    return picked_color
def create_per_line_figure(df, mutation_lst, output_path):
    """
    This function gets a df of freq files, a list of mutations and a path to save graph to.
    """
    # Create a list of the different experiments
    experiments = list(set(df['Experiment']))
    experiments.sort()
    # Create a dictionary to store unique labels and their corresponding colors
    legend_elements = {}
    # Create ggplot alike plot
    plt.style.use('ggplot')
    fig, axes = plt.subplots(layout='constrained',ncols=3, nrows=math.ceil(len(experiments)/3),figsize=(20, 10))
    axes = axes.flatten()
    # Adding the different graphs looping over experiments
    for experiment, ax in zip(experiments, axes):
        df_exp = df[df['Experiment'] == experiment].copy()
        df_exp = df_exp.sort_values('Passage')
        for m in mutation_lst:
            df_exp_mutation = df_exp[df_exp['Full Mutation'] == m].sort_values('Passage').copy()
            # Adding a unique color to the mutation according to a pre-made dictionary - 'COLORS'
            if m in COLORS:
                ax.plot('Passage', 'frequency', data=df_exp_mutation, linestyle='-', marker='.', label=m, color=relevant_mut_colors[m])
                color = COLORS[m]
            else:
                color = get_color_for_new_mutation(m)
                ax.plot('Passage', 'frequency', data=df_exp_mutation, linestyle='-', marker='.', label=m, color=color)
            # Add the label and color to the legend_elements dictionary
            legend_elements[m] = color
        # Tweaking and fixing the title and axis
        ax.set_title(experiment, fontsize='large')
        ax.set_xlabel('Passage', fontsize='small', color='black')
        ax.set_ylim(0, 1)
    # Cover the bottom left and bottom center subplots with a white rectangle
    axes[-2].axis('off')
    axes[-1].axis('off')
    # Making sure the y-axis presented only once
    axes[3].set_ylabel('Frequency', fontsize='small', color='black')
    # Create custom legend elements using Line2D
    custom_lines = [Line2D([0], [0], color=color, lw=2) for label, color in legend_elements.items()]
    # Add the custom legend to the figure
    fig.legend(custom_lines, legend_elements.keys(), loc='lower right', ncol=4, fontsize='large')
    # Saving and presenting the graph
    plt.savefig(output_path, dpi=800)
    plt.show()
    return
def create_per_mutation_figure(df, mutations_list, output_path):
    """
    This function gets a df of freq files, a list of mutations and a path to save graph to.
    """
    # Filter the mutations_list to only include mutations that appear in more than one passage
    mutations_list = [m for m in mutations_list if df[df['Full Mutation'] == m]['Passage'].nunique() > 1]
    num_of_pages = math.ceil(len(mutations_list) / 9)
    # Divide the mutation into sublists of 9 graphs a page
    n = len(mutations_list)
    sublists = []
    for i in range(0, n, 9):
        sublist = mutations_list[i:i + 9]
        sublists.append(sublist)
    for page in range(num_of_pages):
        plt.style.use('ggplot')
        # Create a subplot for each mutation
        fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(20, 10))
        axes = np.array(axes).flatten()  # Flatten the axes array
        # If there's only one mutation, axes will not be an array, so I convert it to an array
        if len(sublists[page]) == 1:
            axes = [axes]
        for mutation, ax in zip(sublists[page], axes):
            # Filter the DataFrame for rows where 'Full Mutation' is equal to the mutation and create a copy
            df_mutation = df[df['Full Mutation'] == mutation].reset_index(drop=True)
            # Create a line plot for each unique combination of 'Line', 'MOI', and 'Done_by'
            sns.lineplot(x='Passage', y='frequency', hue='Experiment', ax=ax, data=df_mutation)
            # Set the title of the subplot to the mutation
            ax.set_title(mutation)
            # Set the y-limit to [0, 1]
            ax.set_ylim(0, 1)
            # Set x-axis to only contain integers
            ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
            # Set legend font size
            ax.legend(fontsize='small')
        plt.tight_layout()
        # Save the figure with a lower dpi
        plt.savefig(output_path+"_"+str(page), dpi=800)
        # Show the figure
        plt.show()
    return
def create_heatmap_figure(df, mutation_lst, output_path):
    # Filter the dataframe
    #mutatin lsr
    df = df[(df['Passage'] == 0) & (df['Full Mutation'].isin(mutation_lst))]
    # Pivot the dataframe
    df_pivot = df.pivot(index='Experiment', columns='Full Mutation', values='frequency')
    # Create the heatmap
    plt.figure(figsize=(15, 12))
    sns.heatmap(df_pivot, cmap=sns.cubehelix_palette(as_cmap=True), vmin=0, vmax=0.1)
    # Save and show the figure
    plt.savefig(output_path)
    plt.show()
    return
def create_genome_map_figure(df, mutation_lst, expe_col, output_path):
    # Filter the dataframe
    df = df[df['Full Mutation'].isin(mutation_lst)]
    # Group the dataframe and take the maximum frequency for each group
    df_grouped = df.groupby(['Experiment', 'ref_pos'])['frequency'].max().reset_index()
    # Create a list of unique colors for each experiment
    unique_experiments = df_grouped['Experiment'].unique()
    # Create the scatter plot
    fig, ax = plt.subplots(figsize=(15, 4))
    # Define the x-axis ranges for the different genes
    gene_ranges = [(1, 1311), (1311, 1727), (1678, 1905), (1761, 3650)]
    # Define the colors for the different genes
    gene_colors = ['lightsteelblue', 'darkorchid', 'royalblue', 'mediumaquamarine']
    # Define the labels for the different genes
    gene_labels = ['mat', 'cp', 'lys', 'rep']
    # Fill the background colors for the different genes and create custom legend handles
    gene_legend_handles = []
    for gene_range, gene_color, gene_label in zip(gene_ranges, gene_colors, gene_labels):
        ax.axvspan(gene_range[0], gene_range[1], facecolor=gene_color, alpha=0.2)
        gene_legend_handles.append(Patch(facecolor=gene_color, alpha=0.2, label=gene_label))
    for experiment in unique_experiments:
        df_experiment = df_grouped[df_grouped['Experiment'] == experiment]
        ax.scatter(df_experiment['ref_pos'], df_experiment['frequency'], color=expe_col[experiment], label=experiment)
    ax.set_xlim(1, 3650)
    ax.set_xlabel('ref_pos')
    ax.set_ylabel('frequency')
    ax.set_title('Genome Map')
    # Create a legend for the experiments on the left side
    experiment_legend_handles = [Patch(facecolor=exp_col[experiment], label=experiment) for experiment in unique_experiments]
    legend1 = ax.legend(handles=experiment_legend_handles, title='Experiment', loc='center left', bbox_to_anchor=(-0.17, 0.5))
    # Create a legend for the background colors on the right side
    legend2 = ax.legend(handles=gene_legend_handles, title='Genes', loc='center left', bbox_to_anchor=(1.01, 0.5))
    # Add both legends to the plot
    ax.add_artist(legend1)
    ax.add_artist(legend2)
    # Save the figure
    plt.savefig(output_path, dpi=800)
    # Show the figure
    plt.show()
    return


## Constants and Parameters
carmel_path = "/Users/adibnzv/Desktop/DESKTOP/SCHOOL/PhD/Year 1/Semester 2/Python for Biologists/Final Project/PythonCourse_Final_Project/DATA/Carmel/Carmel1.csv"
shir_path = "/Users/adibnzv/Desktop/DESKTOP/SCHOOL/PhD/Year 1/Semester 2/Python for Biologists/Final Project/PythonCourse_Final_Project/DATA/Shir/Shir1.csv"
export_path = "/Users/adibnzv/Desktop/DESKTOP/SCHOOL/PhD/Year 1/Semester 2/Python for Biologists/Final Project/PythonCourse_Final_Project/Export/"
min_freq = 0.05
min_cov = 100
relevant_mut_colors = COLORS

## Main Code
# Load the initial freq.csv files into a pandas df
df_carmel = pd.read_csv(carmel_path)
df_shir = pd.read_csv(shir_path)

# Fix Shir's 'File' column to match Carmel's
# Remove the hyphen after the number
df_shir['File'] = df_shir['File'].str.replace('-', '', n=1)
# Insert "_moi10_" after the capitalized letter
df_shir['File'] = df_shir['File'].str.replace('([A-Z])', r'\1_moi10_', n=1)

# Clean up the freq files and add relevant data
df_arr_carmel = df_cleanup(df_carmel, "Carmel")
df_arr_shir = df_cleanup(df_shir, "Shir")

# Create a Joined Freq File and remove lines where base_count is zero
joined_freq = pd.concat([df_arr_carmel, df_arr_shir])
joined_freq = joined_freq[joined_freq['base_count'] != 0].reset_index(drop=True)

# Create Uniq Identifier for each Experiment
joined_freq['Experiment'] = joined_freq['Done_by'] + "-" + \
                            joined_freq['MOI'].astype(str) + "-" \
                            + joined_freq['Line'].astype(str)

# Add a Mutation column
Mutation_df = get_mut_column(joined_freq)

# create a list of mutation that met the cutoffs (can be found in the parameters section)
mut_lst = mut_cutoffs(Mutation_df,min_cov,min_freq)

# create a dictionary to color-code the different experiments
exp_col = {'Carmel-1-A': 'brown', 'Carmel-1-B': 'rosybrown', 'Carmel-10-A': 'darkgreen',
           'Carmel-10-B': 'seagreen', 'Shir-10-A': 'silver',
           'Shir-10-B': 'gray', 'Shir-10-C': 'black'}

# Create Graph per Line and save them to the Export folder:
#create_per_line_figure(Mutation_df, mut_lst, export_path + 'Figure1')

# Create Graph per Mutation and save them to the Export folder:
#create_per_mutation_figure(Mutation_df, mut_lst, export_path + 'Figure2')

# Create Heatmap for passage 0
create_heatmap_figure(Mutation_df, mut_lst, export_path + 'Figure3')

# Create a graph of position of mutation along the genome of MS2
#create_genome_map_figure(Mutation_df, mut_lst,exp_col, export_path + 'Figure4')