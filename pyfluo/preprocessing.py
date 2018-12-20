# Functions and classes for working with fluorescence timeseries data
#
# Sergey Antopolskiy (s.antopolsky@gmail.com), 2018-2019

import pandas as pd
import numpy as np
import logging

def read_xlsx(fname, header=None, usecols='B:ZZZ', read_params=True):
    dF = pd.read_excel(fname, header=header, usecols=usecols).T

    if read_params:
        params = pd.read_excel(fname, sheet_name=1, header=None, 
                               usecols='A:B', dtype=object)
        # transform parameters from DataFrame to dict
        params = params.set_index(0).to_dict()[1]
        # execute strings, if possible (helps parse lists 
        # and may work with other data types)
        params = {k:eval_str(v) for k,v in params.items()}

        return dF, params
    else:
        return dF

def normalize_on_background(dF, n_frames_in_background):
    '''
    Normalizes a trace `dF` on the first `n_frames_in_background` samples,
    using the following formula: (dF - MEAN(bg)) / MEAN(bg) * 100.
    Thus the result is in % of MEAN deviations from the MEAN of bg.
    
    Returns
    dF_norm : normalized trace
    bg_mean : mean of background period
    bg_std : STD of background period
    '''
    bg_mean = dF.loc[:,:n_frames_in_background].mean(axis=1)
    dF_norm = dF.subtract(bg_mean, axis='rows').divide(bg_mean, axis='rows')*100
    return dF_norm, bg_mean

# def clip_backgroud(dF_norm_bg, std_clip=2):
#     #TODO: Shouls I use median here instead of the mean and MAD instead of STD?
#     # The data are not normally distributed
#     dF_norm_bg = dF_norm_bg.copy()
#     for cell_id, cell_signal in dF_norm_bg.iterrows():
#         dF_norm_bg.loc[cell_id,:] = np.clip(
#             cell_signal, 
#             a_min=cell_signal.mean() - std_clip*cell_signal.std(),
#             a_max=cell_signal.mean() + std_clip*cell_signal.std())
#     return dF_norm_bg

def clip_backgroud(dF, ix_to_clip, multiplier=2, agg_function=np.mean, var_function=np.std):
    logging.info(
        f'Clipping signals in range {ix_to_clip} using {agg_function.__name__}+/-{var_function.__name__}*{multiplier}')

    dF_clipped = dF.copy()
    aggregates = dF_clipped.loc[:,ix_to_clip[0]:ix_to_clip[1]].apply(agg_function, axis=1)
    variances = dF_clipped.loc[:,ix_to_clip[0]:ix_to_clip[1]].apply(var_function, axis=1)
    part_clipped = dF_clipped.loc[:,ix_to_clip[0]:ix_to_clip[1]]\
        .clip(lower=aggregates-variances*multiplier, 
              upper=aggregates+variances*multiplier, axis=0)
    dF_clipped.loc[:,ix_to_clip[0]:ix_to_clip[1]] = part_clipped
    return dF_clipped

def check_similar_traces(dF, threshold=0.5):
    order = np.argsort(dF.mean(axis=1))
    df_sorted = dF.iloc[order,:]
    for (i1,s1),(i2,s2) in zip(df_sorted.iloc[:-1].iterrows(), 
                               df_sorted.iloc[1:].iterrows()):
        if (s1==s2).mean() > threshold:
            fig, ax = plt.subplots(1,3,figsize=(14,5))
            ax[0].plot(s1)
            ax[0].set_title('Cell number {}'.format(i1))
            ax[1].plot(s2)
            ax[1].set_title('Cell number {}'.format(i2))
            ax[2].plot(s1)
            ax[2].plot(s2)
            ax[2].set_title('Cells {} and {} overlayed'.format(i1,i2))
            raise ValueError('Cells with indexed {} and {} share over {:.0f}% of values'
                             .format(i1,i2,threshold*100))

def exclude_cells(dF, ixs_exclude=[]):
    return dF.drop(ixs_exclude, axis='rows')

#######
# Utils
#######

def eval_str(inp, verbose=False):
    if type(inp)==str:
        try:
            inp = eval(inp)
        except Exception as e:
            if verbose:
                print('Did not execute the string, because:',str(e))
    return inp

#####
# WIP
#####

class Config():
    '''
    A class containing methods for interacting with the configuration
    parameters for the script.
    '''
    
    def __init__(self, params):
        logging.info('Creating config')
        self.params = params
        self.__dict__.update(params)
#         self.infer_missing_values()
#         self.check_configuration_consistency()
        self.extract_stimuli_info()
    
    def extract_stimuli_info(self):
        # get a list of stimuli names, defined in the config
        self.stim_names = ['_'.join(s.split('_')[1:-1]) for s in self.params.keys() if 'stim_' in s]
        self.stim_names = np.unique(self.stim_names)
        self.stim_n = len(self.stim_names)
        
        # get the list of stimuli initiation frames
        self.stim_begins = [self.params['stim_{}_begin'.format(s)] for s in self.stim_names]
        self.stim_ends = [self.params['stim_{}_end'.format(s)] for s in self.stim_names]
        
#     def __call__(self, name, what):
#         return params['stim_{}_{}'.format(name, what)]
    
#     def __repr__(self):
#         return '(Stim, begin, end): {}'.format(list(zip(self.names, self.begins, self.ends)))
        
# config = Config(params)

## WIP: Version of remove drift (old algorithm from Madina's script) for dataframes. Incomplete.

# std_factor=3 
# n_frames_baseline=30 
# n_frames_skip=20
    
# no_spikes_data = pd.DataFrame(0, index=dF_norm.index, columns=dF_norm.columns)
# movingav = pd.DataFrame(0, index=dF_norm.index, columns=dF_norm.columns)
# result = pd.DataFrame(0, index=dF_norm.index, columns=dF_norm.columns)

# no_spikes_data.iloc[:,:n_frames_skip] = dF_norm.loc[:,:n_frames_skip]

# for c in dF_norm.index:

#     threshold = std_factor * np.std(dF_norm.loc[c,:n_frames_baseline])
#     for t in dF_norm.columns[n_frames_skip:]:
# #         import pdb; pdb.set_trace()
#         abs1=abs(dF_norm.loc[c,t] - dF_norm.loc[c,t-1])
#         abs2=abs(dF_norm.loc[c,t] - no_spikes_data.loc[c,t-1])
#         if (abs1 > threshold) or (abs2 > threshold):
#             no_spikes_data.loc[c,t] = np.mean(no_spikes_data.loc[c,(t-4):t])
#         else:
#             no_spikes_data.loc[c,t] = dF_norm.loc[c,t]

#     movingav.loc[c,:]=convolve(no_spikes_data.loc[c,:],[0.1]*10,mode='same')
#     result.loc[c,:]=dF_norm.loc[c,:]-movingav.loc[c,:]

