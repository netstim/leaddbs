import matplotlib.pyplot as plt
import numpy as np
import os
import h5py
#import seaborn as sns
#sns.set()
import sys

os.environ['PATIENTDIR'] = '/opt/Patient'
sys.path.insert(1, os.environ['PATIENTDIR'])


# those inputs would come from Lead-DBS
#index_side instead of stim_folder
import GUI_inp_dict
from GUI_inp_dict import d

# choose the output folder depending on the simulated hemisphere
if d['Stim_side'] == 0:
    stim_folder = 'Results_rh/'
else:
    stim_folder = 'Results_lh/'


# later we should store connectome name in GUI_inp_dict.py
# load parameters from the file prepared by Lead-DBS
file_inp = h5py.File(os.environ['PATIENTDIR'] + '/oss-dbs_parameters.mat', mode='r')
array_ascii = file_inp['settings']['connectome'][:]
list_ascii = []
for i in range(array_ascii.shape[0]):
    list_ascii.append(array_ascii[i][0])
# list_ascii = map(lambda s: s.strip(), list_ascii)
Connectome_name = ''.join(chr(i) for i in list_ascii)

Projections = []
for i in range(len(file_inp['settings']['connectomeTractNames'][0])):
    ext_string = file_inp[file_inp['settings']['connectomeTractNames'][0][i]]
    list_ascii = []
    for i in range(ext_string.shape[0]):
        list_ascii.append(ext_string[i][0])
    # list_ascii = map(lambda s: s.strip(), list_ascii)
    projection_name = ''.join(chr(i) for i in list_ascii)
    # print(projection_name)
    Projections.append(projection_name)

# It will always be Multi-Tract
# 'Multi-tract' connectomes contain multiple pathways in separate .mat files
if 'Multi-Tract' in Connectome_name:
    Full_paths = [
        os.environ['PATIENTDIR'] + '/' + Connectome_name.rsplit(' ', 1)[1] + '/data' + str(d['Stim_side'] + 1) + '.mat']
    file = h5py.File(Full_paths[0], mode='r')
    number_original = []
    name_original = []
    # run this in loop, later you will check if you actually have this pathways in Summary_status.h5
    for projection_name in Projections:
        if 'origNum' in file[projection_name]:
            number_original.append(file[projection_name]['origNum'][:][0][0])
            name_original.append(projection_name)
        #else: # that means no fibers of the pathway were passed from Lead-DBS
        #    number_original.append(None)

else:
    print('This connectome has only one pathway, exiting')
    raise SystemExit


#print(number_original)
#print(name_original)

# now get the results
hf = h5py.File(os.environ['PATIENTDIR'] + '/' + stim_folder + '/Summary_status.h5', 'r')
lst = list(hf.keys())

Summary_status_all = np.zeros((len(lst), 5), float)
number_original_filtered = []

for pop_i in range(len(lst)):

    #hf = h5py.File(os.environ['PATIENTDIR'] + '/' + stim_folder + 'Summary_status.h5', 'r')
    Summary_status = hf.get(lst[pop_i])
    Summary_status_all[pop_i, :] = np.array(Summary_status)
    print(str(lst[pop_i]))

    index_pathway = name_original.index(lst[pop_i])
    number_original_filtered.append(number_original[index_pathway])

hf.close()

number_original_filtered = np.array(number_original_filtered)
pos = np.arange(len(lst))  # the x locations for the groups
pos_adjusted = pos - 0.25

#fig, ax = plt.subplots(figsize=(len(lst) / 2.0, len(lst) / 8.0))
fig, ax = plt.subplots(figsize=(12, 3))
width = 0.125
colors_opt = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']
labels_opt = ['Activated', 'Not Activated', 'Encapsulation', 'CSF', 'Other']
#bars_status =
for i in range(Summary_status_all.shape[1]):
    ax.bar(pos_adjusted - 0.250 + i*width, Summary_status_all[:, i] / number_original_filtered[:], width, color=colors_opt[i], label=labels_opt[i])
ax.legend(loc='upper right',bbox_to_anchor=(0, 1.25, 1, 0), ncol=5,mode="expand", borderaxespad=0.)
ax.set_xticks(pos_adjusted)
ax.set_xticklabels((lst))
plt.xticks(rotation = 45)
fig.tight_layout()
plt.savefig(os.environ['PATIENTDIR'] + '/Images/Percent_activation_profile.png', format='png',
                    dpi=1000)