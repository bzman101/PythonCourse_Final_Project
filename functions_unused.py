def create_per_Line_figure2(df, mutation_lst ,output_path):
    """
     This function gets a df of freq files and a path to save graph to.
     """
    experiments = list(set(df['Experiment']))
    plt.style.use('ggplot')
    # Create a subplot for each of the Lines
    fig, axes = plt.subplots(nrows=len(experiments), ncols=1, figsize=(5, 2.5 * len(experiments)))
    # If there's only one mutation, axes will not be an array, so we convert it to an array
    if len(experiments) == 1:
        axes = [axes]
    for experiment, ax in zip(experiments, axes):
        # Filter the DataFrame for rows where 'Experiment' is equal to the experiment and create a copy
        df_exp = df[df['Experiment'] == experiment].copy()
        # Create a line plot for each experiment
        for mutation in mutation_lst:
            df_exp_mut = df_exp[df_exp['Full Mutation'] == mutation].copy()
            if mutation in COLORS.keys():
                sns.lineplot(x='Passage', y='frequency', data=df_exp_mut, ax=ax, palette = COLORS, hue = 'Full Mutation')
            else:
                sns.lineplot(x='Passage', y='frequency', data=df_exp_mut, ax=ax, hue='Full Mutation')
            # Set the title of the subplot to the mutation
            ax.set_title(experiment)
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


def create_per_mutation2_figure(df, mutations_list, output_path):
    """
    This function gets a df of freq files, a list of mutations, and a path to save the graph to.
    """
    plt.style.use('ggplot')
    # Filter the mutations_list to only include mutations that appear in more than one passage
    mutations_list = [m for m in mutations_list if df[df['Full Mutation'] == m]['Passage'].nunique() > 1]

    # Create a figure and axis objects for the subplots
    fig, axes = plt.subplots(nrows=len(mutations_list), ncols=3, figsize=(15, 2.5 * len(mutations_list)))

    # If there's only one mutation, axes will not be an array, so we convert it to an array
    if len(mutations_list) == 1:
        axes = [axes]

    for i, mutation in enumerate(mutations_list):
        # Filter the DataFrame for rows where 'Full Mutation' is equal to the mutation and create a copy
        df_mutation = df[df['Full Mutation'] == mutation].reset_index(drop=True)

        # Check if the DataFrame has enough data points to create the line plot
        if len(df_mutation) >= 2:
            # Create a line plot for each unique combination of 'Line', 'MOI', and 'Done_by'
            for experiment, group_data in df_mutation.groupby('Experiment'):
                axes[i, 0].plot(group_data['Passage'], group_data['frequency'], label=experiment)

            # Set the title of the subplot to the mutation
            axes[i, 0].set_title(mutation)
            # Set the y-limit to [0, 1]
            axes[i, 0].set_ylim(0, 1)
            # Set x-axis to only contain integers
            axes[i, 0].xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
            # Set legend font size
            axes[i, 0].legend(fontsize='small')

            # Hide other two empty subplots
            axes[i, 1].axis('off')
            axes[i, 2].axis('off')
        else:
            # If not enough data points, hide all three subplots
            axes[i, 0].axis('off')
            axes[i, 1].axis('off')
            axes[i, 2].axis('off')

    # Adjust the layout
    plt.tight_layout()

    # Save the figure with a lower dpi
    plt.savefig(output_path, dpi=800)

    # Show the figure
    plt.show()

