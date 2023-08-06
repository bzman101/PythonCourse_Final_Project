## Libraries
import pandas as pd
import random
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as ticker
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
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
def get_mut_column(merged_df):
    """
    This Function receives a merged df of freq files adds a Mutation field:
    Mutation - the mutation written in a formatted way :(N_ref+position+N_mut)
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
    return relevant_mutations
def assign_colors_to_experiment(df):
    """
    Assigns random colors to each experiment in the df and creates and returns a dictionary.
    """
    exp_lst = df['Experiment'].tolist()
    unique_items = set(exp_lst)
    num_items = len(unique_items)
    # Get the list of named colors from Matplotlib
    named_colors = list(mcolors.CSS4_COLORS.keys())
    # Generate a list of random colors from the named colors
    random_colors = [named_colors[random.randint(0, len(named_colors)-1)] for _ in range(num_items)]
    color_assignment = {}
    for index, item in enumerate(unique_items):
        color_assignment[item] = random_colors[index]
    return color_assignment
def create_per_Line_figure(df, mutation_lst ,output_path):
    """
    This function gets a df of freq files, a list of mutations and a path to save graph to.
    """
    # Create a list of the different experiments
    experiments = list(set(df['Experiment']))
    # Create ggplot alike plot
    plt.style.use('ggplot')
    fig, axes = plt.subplots(layout='constrained',nrows=(len(experiments)),figsize=(6, 4 *len(experiments)))
    axes = axes.flatten()
    # Adding the different graphs looping over experiments
    for experiment, ax in zip(experiments, axes):
        df_exp = df[df['Experiment'] == experiment].copy()
        df_exp = df_exp.sort_values('Passage')
        for m in mutation_lst:
            df_exp_mutation = df_exp[df_exp['Full Mutation'] == m].sort_values('Passage').copy()
            # Adding a unique color to the mutation according to a pre-made dictionary - 'COLORS'
            if m in COLORS:
                ax.plot('Passage', 'frequency', data=df_exp_mutation, linestyle='-', marker='.', label=m, color=COLORS[m])
            else:
                ax.plot('Passage', 'frequency', data=df_exp_mutation, linestyle='-', marker='.', label=m)
        # Tweaking and fixing the title and axis
        ax.set_title(experiment, fontsize='small')
        ax.set_xlabel('Passage', fontsize='small', color='black')
        ax.set_ylim(0, 1)
    # Making sure the y axis presented only once
    axes[0].set_ylabel('Frequency', fontsize='small', color='black')
    # Adjust the layout
    plt.tight_layout()
    # Adding legend
    plt.legend(ncol=1, loc="center right", borderaxespad=0, facecolor='white', edgecolor='white')
    # Saving and presenting the graph
    plt.savefig(output_path, bbox_inches='tight', dpi=800)
    plt.show()
    return
def create_per_mutation_figure(df, mutations_list, output_path):
    """
    This function gets a df of freq files, a list of mutations and a path to save graph to.
    """
    plt.style.use('ggplot')
    # Filter the mutations_list to only include mutations that appear in more than one passage
    mutations_list = [m for m in mutations_list if df[df['Full Mutation'] == m]['Passage'].nunique() > 1]
    # Create a subplot for each mutation with a smaller size
    fig, axes = plt.subplots(nrows=len(mutations_list), ncols=3, figsize=(5, 2.5 * len(mutations_list)))
    # If there's only one mutation, axes will not be an array, so we convert it to an array
    if len(mutations_list) == 1:
        axes = [axes]
    for mutation, ax in zip(mutations_list, axes):
        # Filter the DataFrame for rows where 'Full Mutation' is equal to the mutation and create a copy
        df_mutation = df[df['Full Mutation'] == mutation].reset_index(drop=True)
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
    plt.tight_layout()
    # Save the figure with a lower dpi
    plt.savefig(output_path, dpi=800)
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
    sns.heatmap(df_pivot, cmap=sns.cubehelix_palette(as_cmap=True))
    # Save and show the figure
    plt.savefig(output_path)
    plt.show()
    return
def create_genome_map_figure(df, mutation_lst, exp_col, output_path):
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
        ax.scatter(df_experiment['ref_pos'], df_experiment['frequency'], color=exp_col[experiment], label=experiment)
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
    plt.savefig(output_path, bbox_inches='tight', dpi=800)
    # Show the figure
    plt.show()
    return


## Constants and Parameters
#maria_path = "/Users/adibnzv/Desktop/DESKTOP/SCHOOL/PhD/Year 1/Semester 2/Python for Biologists/Final Project/PythonCourse_Final_Project/DATA/Maria/Maria1.csv"
carmel_path = "/Users/adibnzv/Desktop/DESKTOP/SCHOOL/PhD/Year 1/Semester 2/Python for Biologists/Final Project/PythonCourse_Final_Project/DATA/Carmel/Carmel1.csv"
shir_path = "/Users/adibnzv/Desktop/DESKTOP/SCHOOL/PhD/Year 1/Semester 2/Python for Biologists/Final Project/PythonCourse_Final_Project/DATA/Shir/Shir1.csv"
export_path = "/Users/adibnzv/Desktop/DESKTOP/SCHOOL/PhD/Year 1/Semester 2/Python for Biologists/Final Project/PythonCourse_Final_Project/Export/"
min_freq = 0.2
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
joined_freq = joined_freq[joined_freq['base_count'] != 0].reset_index(drop=True)

#Create Uniq Identifier for each Experiment
joined_freq['Experiment'] = joined_freq['Done_by'] + "-" + \
                            joined_freq['MOI'].astype(str) + "-" \
                            + joined_freq['Line'].astype(str)

# Add a Mutation column
Mutation_df = get_mut_column(joined_freq)

# create a list of mutation that met the cutoffs (can be found in the parameters section)
mut_lst = mut_cutoffs(Mutation_df,min_cov,min_freq)

# create a dictionary to color-code the different experiments
exp_col = assign_colors_to_experiment(Mutation_df)

# Create Graph per Line and save them to the Export folder:
create_per_Line_figure(Mutation_df, mut_lst, export_path + 'Figure1')

# Create Graph per Mutation and save them to the Export folder:
#create_per_mutation_figure(Mutation_df, mut_lst, export_path + 'Figure2')

# Create Heatmap for passage 0
create_heatmap_figure(Mutation_df, mut_lst, export_path + 'Figure3')

# Create a graph of positon of mutation along the genome of MS2
create_genome_map_figure(Mutation_df, mut_lst,exp_col, export_path + 'Figure4')