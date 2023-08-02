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
