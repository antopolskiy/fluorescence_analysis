from os.path import join
import matplotlib.pyplot as plt
import numpy as np
import logging

def plot_traces_for_all_cells(datasets, dataset_names, conf, save_path=None):
    logging.info(f'''Plotting cells from datasets with names {dataset_names} 
and saving the figures to {save_path}''')
    if type(datasets)!=list:
        datasets = [datasets]
        dataset_names = [dataset_names]

    for cell_id, _ in datasets[0].iterrows():
        ax = plt.figure(figsize=(10,8)).add_subplot(1,1,1)
        for dataset, name in zip(datasets, dataset_names):
            ax.plot(dataset.loc[cell_id].values, label=name)
        
        ## stimuli begin lines and names
        for stim_name, stim_ts in zip(conf.stim_names, conf.stim_begins):
            ax.axvline(stim_ts, color='r', linestyle=':', zorder=0.5)
            ax.text(stim_ts+ax.get_xlim()[1]*0.005, 
                     ax.get_ylim()[1] - np.diff(ax.get_ylim())*0.03, 
                     stim_name)
        
        ax.set_title(f'Cell {cell_id}')
        plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)

        if save_path:
            plt.savefig(join(save_path,f'Cell {cell_id}.png'), dpi=300, bbox_inches='tight')
            plt.close()
        else:
            raise NotImplementedError('Interactive mode for inspecting the figures is not yet implemented')
        break

import inspect

def heatmap(dataset, order_function=np.mean, save_path=None):

    logging.info(f'Plotting heatmap, ordered by {order_function.__name__}, saving into {save_path}')

    if order_function:    
        order = np.argsort(dataset.apply(order_function, axis=1))
        order = dataset.index[order]

    fig, ax = plt.subplots(figsize=(8,14))
    im = ax.imshow(dataset.loc[order], aspect=dataset.shape[1]/dataset.shape[0]*2, cmap='jet')

    ax.set_yticks(np.arange(dataset.shape[0]))
    ax.set_yticklabels(order.values)

    ax.set_ylabel('Cell #')
    ax.set_xlabel('Time, frame #')

    ## to create gridlines on the y axis, need to create minor ticks and then
    ## set their length to 0 so that they don't appear as ticks on the axis
    ax.set_yticks(np.arange(0.5, dataset.shape[0]-1, 1), minor=True);
    ax.grid(which='minor', color='grey', linestyle='-', linewidth=1)
    ax.tick_params(which='minor', length=0)

    fig.tight_layout()

    cbar = ax.figure.colorbar(im, ax=ax, aspect=50)

    if save_path:
        plt.savefig(join(save_path,'heatmap.png'), dpi=300, bbox_inches='tight')