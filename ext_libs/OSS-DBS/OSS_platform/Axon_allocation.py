"""
@author: Konstantin Butenko

This script uses fiber trajectories to allocate axon models

"""

# import tables
import numpy as np
# import nibabel as nib
import os
import h5py
import sys



def fibers_to_axons(name_of_combined_file, name_of_fiber_file, projection_name, axon_model, diam_fib, axon_length,
                    active_contact_coordinates, multiple_projections):

    name_of_directory = name_of_fiber_file.rsplit('/', 1)[0]

    file = h5py.File(name_of_fiber_file, mode='r')
    if multiple_projections == False:
        fiber_array = file['fibers'][:]        # fiber_array has 4 columns (x,y,z,fiber_index), raws - all points
    else:
        fiber_array = file[projection_name]['fibers'][:]

    if fiber_array.ndim == 1:
        print(projection_name, 'projection is empty, check settings for fib. diameter and axon length')
        return 0

    # from nibabel.streamlines.array_sequence import ArraySequence
    from nibabel_SequenceArray import ArraySequence
    streamlines = ArraySequence()
    # streamlines=[]

    N_streamlines = int(fiber_array[3, :].max())  # yes, indexing starts with one in those .mat files

    k = 0
    i_previous = 0
    for i in range(N_streamlines):
        loc_counter = 0

        while ((i + 1) == fiber_array[3, k]):  # this is very slow, you need to extract a pack by np.count?
            k += 1
            loc_counter += 1
            if (k == fiber_array[3, :].shape[0]):
                break

        stream_line = fiber_array[:3, i_previous:i_previous + loc_counter].T
        i_previous = k
        streamlines.append(stream_line)

    streamlines_filtered_ROIs = streamlines  # pre-filtered in Lead-DBS

    if axon_model == 'McIntyre2002':
        from Axon_files.axon import Axon

        param_ax = {
            'centered': True,
            'diameter': diam_fib
        }
        a = Axon(param_ax)
        nr = Axon.get_axonparams(a)
        n_comp = ((nr["ranvier_nodes"] - 1) + nr["inter_nodes"] + nr["para1_nodes"] + nr["para2_nodes"]) / (
                    nr["ranvier_nodes"] - 1)

        ranvier_length, para1_length, para2_length, node_step = (
        nr["ranvier_length"] * 1e-3, nr["para1_length"] * 1e-3, nr["para2_length"] * 1e-3, nr["deltax"] * 1e-3)
        if diam_fib >= 5.7:
            inter_length = (node_step - para1_length * 2 - para2_length * 2) / 6
        else:
            inter_length = (node_step - para1_length * 2 - para2_length * 2) / 3

        n_Ranviers = int(axon_length / node_step)

        n_segments = int((n_Ranviers - 1) * n_comp + 1)  # overall number of points on Axon incl. internodal

        n_total = (n_Ranviers - 1) * n_comp + 1  # total incl. Ranvier
    elif axon_model == 'McIntyre2002_ds':
        from Axon_files.axon import Axon

        param_ax = {
            'centered': True,
            'diameter': diam_fib
        }
        a = Axon(param_ax)
        nr = Axon.get_axonparams(a)
        if diam_fib >= 5.7:
            n_comp = 3  # node -- -- internodal -- -- -- -- internodal -- -- node
        else:
            n_comp = 2  # mode -- -- -- internodal -- -- -- node

        ranvier_length, para1_length, para2_length, node_step = (
        nr["ranvier_length"] * 1e-3, nr["para1_length"] * 1e-3, nr["para2_length"] * 1e-3, nr["deltax"] * 1e-3)
        if diam_fib >= 5.7:
            inter_length = (node_step - para1_length * 2 - para2_length * 2) / 6
        else:
            inter_length = (node_step - para1_length * 2 - para2_length * 2) / 3

        n_Ranviers = int(axon_length / node_step)
        n_segments = int((n_Ranviers - 1) * n_comp + 1)  # overall number of points on Axon incl. internodal
        n_total = (n_Ranviers - 1) * n_comp + 1  # total incl. Ranvier
    elif axon_model == 'Reilly2016':
        n_comp = 2  # only nodes and one internodal per segment

        node_step = diam_fib * 0.2  # from 1 to 2 mm
        n_Ranviers = int(axon_length / node_step)

        n_segments = int((n_Ranviers - 1) * n_comp + 1)
        n_total = (n_Ranviers - 1) * 2 + 1  # one internodal per segment + the last Ranvier
    else:
        print("The neuron model is not implemented")
        raise SystemExit

    # resampling with different number of nodes of Ranvier
    from Arbitrary_streamline_to_Ranviers import length_fiber
    lengths_streamlines_filtered = list(length_fiber(streamlines_filtered_ROIs))

    streamlines_resampled = ArraySequence()
    # streamlines_resampled=[]

    from Arbitrary_streamline_to_Ranviers import resample_streamline_for_Ranvier
    streamline_index = 0
    excluded_streamlines = []
    total_points = 0
    for single_stream_lime in streamlines_filtered_ROIs:
        n_Ranvier_this_axon = int(lengths_streamlines_filtered[streamline_index] / node_step)
        streamline_resampled = resample_streamline_for_Ranvier(streamlines_filtered_ROIs[streamline_index],
                                                               n_Ranvier_this_axon * node_step, n_Ranvier_this_axon)
        if len(streamline_resampled) < n_Ranviers:
            print("streamline ", streamline_index, " is too short")
            excluded_streamlines.append(streamline_index)
        else:
            streamlines_resampled.append(streamline_resampled)
            total_points = total_points + len(streamline_resampled)
        streamline_index = streamline_index + 1
    streamline_index = 0

    # print(len(streamlines_resampled),len(streamlines_filtered_ROIs),n_total)

    #    streamlines_resampled=streamlines_filtered_ROIs # the length check was already done in Lead-DBS

    #    from dipy.viz import ui, window
    #    from dipy.viz import window, actor, colormap as cmap
    #    from nibabel.streamlines.array_sequence import ArraySequence
    #
    #
    #    stream_actor_dm = actor.line(streamlines_resampled,[0.0,0.4470,0.7410])        #particular color
    #    ren = window.Renderer()
    #    #hdp
    #    ren.add(stream_actor_dm)
    #    #ren.add(seedroi_actor)
    #    window.show(ren)

    #    # now we transform it back to .mat for visualization
    #    #point_fibers_array=np.zeros((4,total_points),float)
    #    point_fibers_array=np.zeros((total_points,4),float)
    #    glob_counter=0
    #    for i_steamline in range(len(streamlines_resampled)):
    #
    #        streamline_length = len(streamlines_resampled[i_steamline])
    #
    #
    #        point_fibers_array[glob_counter:glob_counter+streamline_length,:3]=streamlines_resampled[i_steamline][:]
    #        point_fibers_array[glob_counter:glob_counter+streamline_length,3]=i_steamline+1
    #
    #        glob_counter=glob_counter+streamline_length
    #
    #    from scipy.io import savemat
    #    mdic = {"fibers": point_fibers_array, "ea_fibformat": "1.0"}
    #    savemat(name_of_directory+'/'+projection_name+"_resampled_to_nodes.mat", mdic)

    center_of_mass = active_contact_coordinates  # it is a list

    from scipy import spatial

    streamlines_ROI_centered = ArraySequence()
    # streamlines_ROI_centered=[]

    max_A_len = 0
    overall_shape = 0

    for inx_axn in range(len(streamlines_resampled)):

        single_streamline_ROI_centered = np.zeros((n_Ranviers, 3), float)

        A = streamlines_resampled[inx_axn]

        if inx_axn > 600:
            for k in range(A.shape[0]):
                if A[k, 2] == 0.0:
                    print(inx_axn)

        overall_shape = overall_shape + A.shape[0]
        if max_A_len < A.shape[0]:
            max_A_len = A.shape[0]
        # A[spatial.KDTree(A).query(pt)[1]]
        # distance,index = spatial.KDTree(A).query(center_of_mass)        #distance is a local index of closest node of Ranvier on the axon. You should check for several center_of_mass as you do already when seeding for the network

        distance_list = []
        index_list = []
        for j in range(len(center_of_mass)):
            distance, index = spatial.KDTree(A).query(
                center_of_mass[j])  # distance is a local index of closest node of Ranvier on the axon
            distance_list.append(distance)
            index_list.append(index)

        index = index_list[distance_list.index(min(distance_list))]  # index of the closest point as assigned as index
        distance = min(distance_list)

        loc_index = 0

        if index < int(n_Ranviers / 2):

            bias_to_right = int(n_Ranviers / 2) - index
            # for i in xrange(0,index+int(n_Ranviers/2)+bias_to_right):
            for i in range(0, int(n_Ranviers)):
                single_streamline_ROI_centered[loc_index, :] = A[i]
                loc_index = loc_index + 1
        elif index + int(n_Ranviers / 2) + 1 > A.shape[0]:
            # bias_to_left=index+int(n_Ranviers/2)-A.shape[0]
            # for i in xrange(index-int(n_Ranviers/2)-bias_to_left,index+int(n_Ranviers/2)-bias_to_left):
            # for i in range(0,A.shape[0]-n_Ranviers,A.shape[0]):
            for i in range(A.shape[0] - n_Ranviers, A.shape[0]):
                single_streamline_ROI_centered[loc_index, :] = A[i]
                loc_index = loc_index + 1
        else:
            if n_Ranviers % 2 == 0:
                for i in range(index - int(n_Ranviers / 2), index + int(n_Ranviers / 2)):
                    # for i in xrange(index-int(n_Ranviers/2),index+int(n_Ranviers/2)):
                    single_streamline_ROI_centered[loc_index, :] = A[i]
                    loc_index = loc_index + 1
            else:
                for i in range(index - int(n_Ranviers / 2), index + int(n_Ranviers / 2) + 1):
                    # for i in xrange(index-int(n_Ranviers/2),index+int(n_Ranviers/2)):
                    single_streamline_ROI_centered[loc_index, :] = A[i]
                    loc_index = loc_index + 1

        streamlines_ROI_centered.append(single_streamline_ROI_centered)

    if len(streamlines_ROI_centered) != len(streamlines_resampled):
        print("Failed to sample some axons!")
        raise SystemExit

    #    from dipy.viz import ui, window
    #    from dipy.viz import window, actor, colormap as cmap
    #    from nibabel.streamlines.array_sequence import ArraySequence
    #
    #
    #    stream_actor_dm = actor.line(streamlines_ROI_centered,[0.0,0.4470,0.7410])        #particular color
    #    ren = window.Renderer()
    #    #hdp
    #    ren.add(stream_actor_dm)
    #    #ren.add(seedroi_actor)
    #    window.show(ren)

    '''streamlines_ROI_centered should already contain the position of Ranvier nodes. Now we get internodal compartments'''

    def normalized(a, axis=-1, order=2):
        l2 = np.atleast_1d(np.linalg.norm(a, order, axis))
        l2[l2 == 0] = 1
        return a / np.expand_dims(l2, axis)

    n_total = int(n_total)
    Array_coord = np.zeros((n_total, 3, len(streamlines_ROI_centered)), float)

    for inx_axn in range(len(streamlines_ROI_centered)):

        if axon_model == 'McIntyre2002':
            inx_comp = 0
            for inx in range(n_Ranviers - 1):  # the last node is not included
                Array_coord[inx_comp, :, inx_axn] = streamlines_ROI_centered[inx_axn][inx]

                if diam_fib >= 5.7:
                    Array_coord[inx_comp + 11, :, inx_axn] = streamlines_ROI_centered[inx_axn][inx + 1]
                    internodal_vector = Array_coord[inx_comp + 11, :, inx_axn] - Array_coord[inx_comp, :, inx_axn]
                else:
                    Array_coord[inx_comp + 8, :, inx_axn] = streamlines_ROI_centered[inx_axn][inx + 1]
                    internodal_vector = Array_coord[inx_comp + 8, :, inx_axn] - Array_coord[inx_comp, :, inx_axn]

                # internodal_vector_normalized=preprocessing.normalize(internodal_vector,norm='l2')
                internodal_vector_normalized = normalized(internodal_vector)

                loc_pos = 0.0
                if diam_fib >= 5.7:
                    for inx_loc in np.arange(1,
                                             n_comp):  # only internodal compartments. The distances will be computed from the node of Ranvier using loc_pos
                        inx_loc = int(inx_loc)
                        if inx_loc == 1:
                            loc_pos = (ranvier_length + para1_length) / 2
                        if inx_loc == 2 or inx_loc == 10:
                            loc_pos = loc_pos + (para1_length + para2_length) / 2
                        if inx_loc == 3 or inx_loc == 9:
                            loc_pos = loc_pos + (para2_length + inter_length) / 2
                        if inx_loc == 4 or inx_loc == 5 or inx_loc == 6 or inx_loc == 7 or inx_loc == 8:
                            loc_pos = loc_pos + (inter_length) / 1  # switch to mm from µm

                        # print(inx_comp,inx_loc,inx_axn)
                        Array_coord[inx_comp + inx_loc, 0, inx_axn] = Array_coord[inx_comp, 0, inx_axn] + loc_pos * \
                                                                      internodal_vector_normalized[0][0]
                        Array_coord[inx_comp + inx_loc, 1, inx_axn] = Array_coord[inx_comp, 1, inx_axn] + loc_pos * \
                                                                      internodal_vector_normalized[0][1]
                        Array_coord[inx_comp + inx_loc, 2, inx_axn] = Array_coord[inx_comp, 2, inx_axn] + loc_pos * \
                                                                      internodal_vector_normalized[0][2]

                    inx_comp = inx_comp + 11

                if diam_fib < 5.7:
                    for inx_loc in np.arange(1,
                                             n_comp):  # only internodal compartments. The distances will be computed from the node of Ranvier using loc_pos
                        if inx_loc == 1:
                            loc_pos = (ranvier_length + para1_length) / 2
                        if inx_loc == 2 or inx_loc == 7:
                            loc_pos = loc_pos + (para1_length + para2_length) / 2
                        if inx_loc == 3 or inx_loc == 6:
                            loc_pos = loc_pos + (para2_length + inter_length) / 2
                        if inx_loc == 4 or inx_loc == 5:
                            loc_pos = loc_pos + (inter_length) / 1  # switch to mm from µm

                        inx_loc = int(inx_loc)
                        Array_coord[inx_comp + inx_loc, 0, inx_axn] = Array_coord[inx_comp, 0, inx_axn] + loc_pos * \
                                                                      internodal_vector_normalized[0][0]
                        Array_coord[inx_comp + inx_loc, 1, inx_axn] = Array_coord[inx_comp, 1, inx_axn] + loc_pos * \
                                                                      internodal_vector_normalized[0][1]
                        Array_coord[inx_comp + inx_loc, 2, inx_axn] = Array_coord[inx_comp, 2, inx_axn] + loc_pos * \
                                                                      internodal_vector_normalized[0][2]

                    inx_comp = inx_comp + 8

        elif axon_model == 'McIntyre2002_ds': # downsampled version
            inx_comp = 0
            for inx in range(n_Ranviers - 1):  # the last node is not included
                Array_coord[inx_comp, :, inx_axn] = streamlines_ROI_centered[inx_axn][inx]

                if diam_fib >= 5.7:
                    Array_coord[inx_comp + 3, :, inx_axn] = streamlines_ROI_centered[inx_axn][inx + 1]
                    internodal_vector = Array_coord[inx_comp + 3, :, inx_axn] - Array_coord[inx_comp, :, inx_axn]
                else:
                    Array_coord[inx_comp + 2, :, inx_axn] = streamlines_ROI_centered[inx_axn][inx + 1]
                    internodal_vector = Array_coord[inx_comp + 2, :, inx_axn] - Array_coord[inx_comp, :, inx_axn]

                # internodal_vector_normalized=preprocessing.normalize(internodal_vector,norm='l2')
                internodal_vector_normalized = normalized(internodal_vector)

                loc_pos = 0.0
                if diam_fib >= 5.7:   # node -- -- internodal -- -- -- -- internodal -- -- node
                    for inx_loc in np.arange(1,
                                             n_comp):  # only internodal compartments. The distances will be computed from the node of Ranvier using loc_pos
                        inx_loc = int(inx_loc)
                        if inx_loc == 1:
                            loc_pos = (ranvier_length + inter_length) / 2 + para1_length + para2_length
                        elif inx_loc == 2:
                            loc_pos = loc_pos + 5 * inter_length
                        else:
                            print('wrong number of compartments')

                        Array_coord[inx_comp + inx_loc, 0, inx_axn] = Array_coord[inx_comp, 0, inx_axn] + loc_pos * \
                                                                      internodal_vector_normalized[0][0]
                        Array_coord[inx_comp + inx_loc, 1, inx_axn] = Array_coord[inx_comp, 1, inx_axn] + loc_pos * \
                                                                      internodal_vector_normalized[0][1]
                        Array_coord[inx_comp + inx_loc, 2, inx_axn] = Array_coord[inx_comp, 2, inx_axn] + loc_pos * \
                                                                      internodal_vector_normalized[0][2]

                    inx_comp = inx_comp + 3

                if diam_fib < 5.7:   # mode -- -- -- internodal -- -- -- node

                    loc_pos = 0.5 * ranvier_length + 1.5 * inter_length + para1_length + para2_length

                    Array_coord[inx_comp + 1, 0, inx_axn] = Array_coord[inx_comp, 0, inx_axn] + loc_pos * \
                                                                  internodal_vector_normalized[0][0]
                    Array_coord[inx_comp + 1, 1, inx_axn] = Array_coord[inx_comp, 1, inx_axn] + loc_pos * \
                                                                  internodal_vector_normalized[0][1]
                    Array_coord[inx_comp + 1, 2, inx_axn] = Array_coord[inx_comp, 2, inx_axn] + loc_pos * \
                                                                  internodal_vector_normalized[0][2]
                    inx_comp = inx_comp + 2

        elif axon_model == 'Reilly2016':
            inx_comp = 0
            for inx in range(n_Ranviers - 1):  # the last node is not included
                Array_coord[inx_comp, :, inx_axn] = streamlines_ROI_centered[inx_axn][inx]

                Array_coord[inx_comp + 2, :, inx_axn] = streamlines_ROI_centered[inx_axn][inx + 1]
                internodal_vector = Array_coord[inx_comp + 2, :, inx_axn] - Array_coord[inx_comp, :, inx_axn]

                internodal_vector_normalized = normalized(internodal_vector)

                # we need only
                loc_pos = node_step * 0.5
                Array_coord[inx_comp + 1, 0, inx_axn] = Array_coord[inx_comp, 0, inx_axn] + loc_pos * \
                                                        internodal_vector_normalized[0][0]
                Array_coord[inx_comp + 1, 1, inx_axn] = Array_coord[inx_comp, 1, inx_axn] + loc_pos * \
                                                        internodal_vector_normalized[0][1]
                Array_coord[inx_comp + 1, 2, inx_axn] = Array_coord[inx_comp, 2, inx_axn] + loc_pos * \
                                                        internodal_vector_normalized[0][2]

                inx_comp = inx_comp + 2
    #

    # np.savetxt('Array_coord.csv', Array_coord, delimiter=" ")

    Array_coord_platform = np.zeros((n_total * len(streamlines_ROI_centered), 3), float)
    Array_coord_colored = np.zeros((n_total * len(streamlines_ROI_centered), 4), float)

    glob_ind = 0
    for axon_index in range(len(streamlines_ROI_centered)):
        Array_coord_platform[glob_ind:glob_ind + n_total, :] = Array_coord[:, :, axon_index]
        Array_coord_colored[glob_ind:glob_ind + n_total, :3] = Array_coord[:, :, axon_index]
        Array_coord_colored[glob_ind:glob_ind + n_total, 3] = axon_index + 1  # because in Matlab they start from 1

        glob_ind = glob_ind + n_total

        # name_of_directory+'/'+projection_name+

    np.savetxt(name_of_directory + '/' + 'Array_coord_colored_' + projection_name + '.csv', Array_coord_colored,
               delimiter=" ")
    np.savetxt(name_of_directory + '/' + 'Array_coord_' + projection_name + '.csv', Array_coord_platform, delimiter=" ")

    # np.savetxt('Array_coord_colored_'+name_of_fiber_file+'.csv', Array_coord_colored, delimiter=" ")
    # np.savetxt('Array_coord_'+name_of_fiber_file+'.csv', Array_coord_platform, delimiter=" ")

    # hf = h5py.File(name_of_combined_file + '.h5', 'a')

    # hf = h5py.File(name_of_directory+'/'+name_of_combined_file + '.h5', 'a')
    # hf.create_dataset(projection_name, data=Array_coord_platform)
    # hf.close()

    hf = h5py.File(os.environ['PATIENTDIR'] + '/' + name_of_combined_file + '.h5', 'a')
    hf.create_dataset(projection_name, data=Array_coord_platform)
    hf.close()

    from scipy.io import savemat
    mdic = {"fibers": Array_coord_colored, "ea_fibformat": "1.0"}
    savemat(name_of_directory + '/' + name_of_combined_file + '_' + projection_name + "_axons.mat", mdic)

    return int(n_Ranviers)


if __name__ == '__main__':

    """ To ensure the proper input, use Lead-DBS to call this function
        args - stim. (patient) directory, hemisphere index (0-rh, 1-lh), 
        
        optional: NEURON model ('McIntyre2002', 'McIntyre2002_ds' (downsampled), 'Reilly2016' (classic McNeil's)) 
                  cluster call (1-yes, 0 - no) """

    import json

    # make sure we are in OSS_platform/
    oss_dbs_folder = os.path.dirname(os.path.realpath(sys.argv[0]))
    os.chdir(oss_dbs_folder)

    patient_dir = sys.argv[1:][0]   # patient dir refers to the stim. folder here
    index_side = int(sys.argv[1:][1])

    # not necessary to use env.variable, but just for consistency
    os.environ['PATIENTDIR'] = patient_dir

    if len(sys.argv[1:]) > 3:
        cluster_run = int(sys.argv[1:][3])

    # load parameters from the file prepared by Lead-DBS
    file_inp = h5py.File(os.environ['PATIENTDIR'] + '/oss-dbs_parameters.mat', mode='r')

    if len(sys.argv[1:]) > 2:
        Axon_model = str(sys.argv[1:][2])
    elif 'neuronModel' in file_inp['settings']:
        array_ascii = file_inp['settings']['neuronModel'][:]
        list_ascii = []
        for i in range(array_ascii.shape[0]):
            list_ascii.append(array_ascii[i][0])
        # list_ascii = map(lambda s: s.strip(), list_ascii)
        Axon_Model = ''.join(chr(i) for i in list_ascii)
        if Axon_Model not in ['McIntyre2002', 'McIntyre2002_ds', 'Reilly2016']:
            print('The selected NEURON models is not recognized, check oss-dbs_parameters.mat')
            raise SystemExit
    else:
        Axon_model = 'McIntyre2002'  # Assume McIntyre by default
            
    array_ascii = file_inp['settings']['connectomePath'][:]
    list_ascii = []
    for i in range(array_ascii.shape[0]):
        list_ascii.append(array_ascii[i][0])
    # list_ascii = map(lambda s: s.strip(), list_ascii)
    Path_to_files = ''.join(chr(i) for i in list_ascii)
    # print('Path_to_files',Path_to_files)

    array_ascii = file_inp['settings']['connectome'][:]
    list_ascii = []
    for i in range(array_ascii.shape[0]):
        list_ascii.append(array_ascii[i][0])
    # list_ascii = map(lambda s: s.strip(), list_ascii)
    Connectome_name = ''.join(chr(i) for i in list_ascii)

    # 'Multi-tract' connectomes contain multiple pathways in separate .mat files
    if 'Multi-Tract' in Connectome_name:
        Full_paths = [
            os.environ['PATIENTDIR'] + '/' + Connectome_name.rsplit(' ', 1)[1] + '/data' + str(index_side + 1) + '.mat']
    else:
        Full_paths = [os.environ['PATIENTDIR'] + '/' + Connectome_name + '/data' + str(index_side + 1) + '.mat']

    if 'Multi-Tract' in Connectome_name:
        Projections = []
        for i in range(len(file_inp['settings']['connectomeTractNames'][0])):
            ext_string = file_inp[file_inp['settings']['connectomeTractNames'][0][i]]
            list_ascii = []
            for j in range(ext_string.shape[0]):
                list_ascii.append(ext_string[j][0])
            # list_ascii = map(lambda s: s.strip(), list_ascii)
            projection_name = ''.join(chr(i) for i in list_ascii)
            # print(projection_name)
            Projections.append(projection_name)

    Name_to_save = 'Allocated_axons'  # hardcoded name for axons pre-filtered by Lead-DBS

    axon_length = list(file_inp['settings']['axonLength'][:][0][:])
    diams_fib = list(file_inp['settings']['fiberDiameter'][:][0][:])

    # check which contacts are active
    ampl_vector = file_inp['settings']['Phi_vector'][:, index_side]
    ampl_vector = list(ampl_vector)
    Active_contact_coordinates = []
    import math

    stimSets = bool(file_inp['settings']['stimSetMode'][0][0])  # if StimSets are used, create a dummy ampl_vector

    if stimSets:  # create ampl_vector from maximum values of the protocols
        stim_protocols = np.genfromtxt(os.environ['PATIENTDIR'] + '/Current_protocols_' + str(index_side) + '.csv',
                                       dtype=float, delimiter=',', names=True)

        total_contacts = len(list(stim_protocols[0]))
        total_protocols = stim_protocols.shape[0]

        protocols_array = np.zeros((total_protocols, total_contacts), float)
        ampl_vector = np.zeros(total_contacts, float)
        ampl_vector = list(stim_protocols[0])  # just initialize

        for j in range(total_protocols):
            protocols_array[j, :] = list(stim_protocols[j])
            for i in range(total_contacts):
                if not math.isnan(protocols_array[j, i]):
                    ampl_vector[i] = 1.0  # you do not need a value, just substitute NaN

    for i in range(len(ampl_vector)):
        if not (math.isnan(ampl_vector[i])):
            a_ref = file_inp['settings']['contactLocation'][index_side][0]
            b = file_inp[a_ref]
            Active_contact_coordinates.append(b[:, i])


    Fiber_names = []
    for fiber_file in Full_paths:  # obsolete way, we store mutiple fibers in one data.mat

        name_of_fiber_file = fiber_file.rsplit('/', 1)[1]
        Fiber_names.append(name_of_fiber_file[:-4])  # cut .mat

    # name_of_fiber_file can be something like ['ansa_lenticularis','cerebellothalamic','medial_lemniscus']

    name_of_directory = fiber_file.rsplit('/', 1)[0]  # they should be in the same directory

    # just an example of parameters that will be filled in
    axon_dict = {
        'Axon_Model_Type': 'No_model',
        'Neuron_model_array_prepared': 1,  # now it is known
        'Name_prepared_neuron_array': '0',
        'diam_fib': [10],
        'n_Ranvier': [17],
    }

    axon_dict['Axon_Model_Type'] = Axon_model
    axon_dict['Name_prepared_neuron_array'] = Name_to_save + '.h5'

    n_Ranviers_per_projection = np.zeros(len(axon_length), int)
    for i in range(len(diams_fib)):
        if 'Multi-Tract' in Connectome_name:
            # print(Projections[i])
            n_Ranviers_per_projection[i] = fibers_to_axons(Name_to_save, Full_paths[0], Projections[i], Axon_model,
                                                           diams_fib[i], axon_length[i], Active_contact_coordinates,
                                                           True)
            print("Projection", Projections[i], "seeded with", n_Ranviers_per_projection[i], "nodes of Ranvier\n")
        else:
            n_Ranviers_per_projection[i] = fibers_to_axons(Name_to_save, Full_paths[i], Fiber_names[i], Axon_model,
                                                           diams_fib[i], axon_length[i], Active_contact_coordinates,
                                                           False)
            print("Projection", Fiber_names[i], "seeded with", n_Ranviers_per_projection[i], "nodes of Ranvier\n")
    np.savetxt(os.environ['PATIENTDIR'] + '/' + Name_to_save + '_N_nodes.csv', n_Ranviers_per_projection, delimiter=" ")

    if len(diams_fib) == 1 and not os.path.exists(os.environ['PATIENTDIR'] + '/' + Name_to_save + '.h5'):  # stupid way
        n_Ranviers_per_projection_list = [n_Ranviers_per_projection[0]]
    elif len(diams_fib) == 1 and os.path.exists(os.environ['PATIENTDIR'] + '/' + Name_to_save + '.h5'):  # stupid way
        n_Ranviers_per_projection_list = [int(n_Ranviers_per_projection + 0)]  # stupid trick
    else:
        n_Ranviers_per_projection_list = list(n_Ranviers_per_projection)

    diams_fib_true = []
    n_Ranviers_per_projection_true = []
    for i in range(len(diams_fib)):
        if n_Ranviers_per_projection[i] != 0:
            diams_fib_true.append(float(diams_fib[i]))
            n_Ranviers_per_projection_true.append(int(n_Ranviers_per_projection[i]))

    axon_dict['n_Ranvier'] = n_Ranviers_per_projection_true
    axon_dict['diam_fib'] = diams_fib_true


    with open(os.environ['PATIENTDIR'] + '/Allocated_axons_parameters.json', 'w') as save_as_dict:
        json.dump(axon_dict, save_as_dict)


