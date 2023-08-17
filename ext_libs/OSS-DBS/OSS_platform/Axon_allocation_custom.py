'''
    By K. Butenko
    Allocate pathways for each trk file and stores in a combined file in .h5 format
'''

# import tables
import numpy as np
# import nibabel as nib
import os
import h5py
import sys

#from dipy.io.streamline import load_tractogram, save_tractogram
#from dipy.io.stateful_tractogram import Space, StatefulTractogram
#from dipy.io.image import load_nifti_data, load_nifti, save_nifti
#from dipy.viz import window, actor, colormap as cmap

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx,array[idx]


def length(xyz, along=False):                   # from Dipy (https://dipy.org/)
    ''' Euclidean length of track line
    This will give length in mm if tracks are expressed in world coordinates.
    Parameters
    ------------
    xyz : array-like shape (N,3)
       array representing x,y,z of N points in a track
    along : bool, optional
       If True, return array giving cumulative length along track,
       otherwise (default) return scalar giving total length.
    Returns
    ---------
    L : scalar or array shape (N-1,)
       scalar in case of `along` == False, giving total length, array if
       `along` == True, giving cumulative lengths.
    Examples
    ----------
    >>> from dipy.tracking.metrics import length
    >>> xyz = np.array([[1,1,1],[2,3,4],[0,0,0]])
    >>> expected_lens = np.sqrt([1+2**2+3**2, 2**2+3**2+4**2])
    >>> length(xyz) == expected_lens.sum()
    True
    >>> len_along = length(xyz, along=True)
    >>> np.allclose(len_along, expected_lens.cumsum())
    True
    >>> length([])
    0
    >>> length([[1, 2, 3]])
    0
    >>> length([], along=True)
    array([0])
    '''
    xyz = np.asarray(xyz)
    if xyz.shape[0] < 2:
        if along:
            return np.array([0])
        return 0
    dists = np.sqrt((np.diff(xyz, axis=0)**2).sum(axis=1))
    if along:
        return np.cumsum(dists)
    return np.sum(dists)

def length_fiber(streamlines, affine=None):   # from Dipy (https://dipy.org/)
    """
    Calculate the lengths of many streamlines in a bundle.
    Parameters
    ----------
    streamlines : list
        Each item in the list is an array with 3D coordinates of a streamline.
    affine : 4 x 4 array
        An affine transformation to move the fibers by, before computing their
        lengths.
    Returns
    -------
    Iterator object which then computes the length of each
    streamline in the bundle, upon iteration.
    """

    return map(length, streamlines)


def index_for_length(xyz,req_length, along=True): #from Dipy (https://dipy.org/)
    """ Euclidean length of track line
    This will give length in mm if tracks are expressed in world coordinates.
    Parameters
    ------------
    xyz : array-like shape (N,3)
       array representing x,y,z of N points in a track
    along : bool, optional
       If True, return array giving cumulative length along track,
       otherwise (default) return scalar giving total length.
    Returns
    ---------
    L : scalar or array shape (N-1,)
       scalar in case of `along` == False, giving total length, array if
       `along` == True, giving cumulative lengths.
    Examples
    ----------
    >>> from dipy.tracking.metrics import length
    >>> xyz = np.array([[1,1,1],[2,3,4],[0,0,0]])
    >>> expected_lens = np.sqrt([1+2**2+3**2, 2**2+3**2+4**2])
    >>> length(xyz) == expected_lens.sum()
    True
    >>> len_along = length(xyz, along=True)
    >>> np.allclose(len_along, expected_lens.cumsum())
    True
    >>> length([])
    0
    >>> length([[1, 2, 3]])
    0
    >>> length([], along=True)
    array([0])
    """
    xyz = np.asarray(xyz)
    if xyz.shape[0] < 2:
        if along:
            return np.array([0])
        return 0

    dists = np.sqrt((np.diff(xyz, axis=0)**2).sum(axis=1))

    if along:
        cummulated_lengths=np.cumsum(dists)
        idx,value=find_nearest(cummulated_lengths, req_length)
        if value>req_length:
            idx=idx-1

        return idx,cummulated_lengths[idx]
    return idx,cummulated_lengths[idx]

def resample_streamline_for_Ranvier(streamline_array,axon_length,n_Ranviers):
    cut_index,cummulated_length=index_for_length(streamline_array,axon_length)        #after this index we cut the streamline

    streamline_array_Ranvier=np.zeros((cut_index+1+1+1,3),float)         #check notes in Cicero! Don't mix up sums and positions. +1 for the last Ranvier node, +1 for the sum, +1 for index
    last_segment_length=axon_length-cummulated_length

    #print "Last_point_from_the_streamline: ",streamline_array[cut_index+1,0],streamline_array[cut_index+1,1],streamline_array[cut_index+1,2]

    x_vect=streamline_array[cut_index+1+1,0]-streamline_array[cut_index+1,0]    #check in between the last taken and the next one
    y_vect=streamline_array[cut_index+1+1,1]-streamline_array[cut_index+1,1]
    z_vect=streamline_array[cut_index+1+1,2]-streamline_array[cut_index+1,2]
    v=np.array([x_vect,y_vect,z_vect])

    v_hat = v / (v**2).sum()**0.5

    streamline_array_Ranvier[:cut_index+1+1,:]=streamline_array[:cut_index+1+1,:]

    streamline_array_Ranvier[cut_index+1+1,:]=last_segment_length*v_hat+streamline_array[cut_index+1,:]

    #print streamline_array_Ranvier

    from dipy.tracking.streamline import set_number_of_points
    #from streamlinespeed import set_number_of_points

    streamline_resampled = set_number_of_points(streamline_array_Ranvier, nb_points=n_Ranviers)
    #streamline_resampled =seq_sampling(streamline_array_Ranvier, res=n_Ranviers)

    return streamline_resampled



def fibers_to_axons(combined_h5_file, pathway_mat_file, axon_model, diam_fib, axon_length,
                    centering_coordinates, projection_name=None):

    if projection_name == None:
        # extract it from trk file name
        projection_name = pathway_mat_file.rsplit('/', 1)[1][0:-4]

    name_of_directory = combined_h5_file.rsplit('/', 1)[0]
    print(name_of_directory)


    file = h5py.File(pathway_mat_file, mode='r')
    fiber_array = file['fibers'][:]        # fiber_array has 4 columns (x,y,z,fiber_index), raws - all points

    if fiber_array.ndim == 1:
        print(projection_name, 'projection is empty')
        return 0

    from nibabel.streamlines.array_sequence import ArraySequence
    #from nibabel_SequenceArray import ArraySequence
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
    #from Arbitrary_streamline_to_Ranviers import length_fiber
    lengths_streamlines_filtered = list(length_fiber(streamlines_filtered_ROIs))

    streamlines_resampled = ArraySequence()
    # streamlines_resampled=[]

    #from Arbitrary_streamline_to_Ranviers import resample_streamline_for_Ranvier
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

    center_of_mass = centering_coordinates  # it is a list of lists (multiple 3-D coordinates allowed)

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

        index = index_list[distance_list.index(min(distance_list))]   # index of the closest point as assigned as index
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

    hf = h5py.File(combined_h5_file + '.h5', 'a')
    hf.create_dataset(projection_name, data=Array_coord_platform)
    hf.close()

    #from scipy.io import savemat
    #mdic = {"fibers": Array_coord_colored, "ea_fibformat": "1.0"}
    #savemat(name_of_directory + '/' + name_of_combined_file + '_' + projection_name + "_axons.mat", mdic)

    return int(n_Ranviers), projection_name


if __name__ == '__main__':

    # Inputs
    combined_h5_file = '/home/konstantin/Documents/Codes/Useful_Py_scripts/Custom_mat_to_h5_pathways/all_tracts'
    pathway_mat_files = ['/home/konstantin/Documents/GitHub/leaddbs/connectomes/dMRI_MultiTract/PetersenOldNew/SMA_hdp_left.mat',
                         '/home/konstantin/Documents/GitHub/leaddbs/connectomes/dMRI_MultiTract/PetersenOldNew/SMA_hdp_right.mat',
                         '/home/konstantin/Documents/GitHub/leaddbs/connectomes/dMRI_MultiTract/PetersenOldNew/gpe2stn_sm_left.mat']

    pathway_names = ['tract1','tract2','tract3']
    diam_fib = [5.7,5.7,3.0]
    axon_length = [20.0,20.0,10.0]

    centering_coordinates = [[7.5838, -18.3984, 1.8932],[-7.5838, -18.3984, 1.8932]]  # in this case, we just have some STN coordinates for left and right
    axon_model = 'McIntyre2002_ds'



    n_Ranviers_per_projection = np.zeros(len(axon_length), int)
    for i in range(len(pathway_mat_files)):
        n_Ranviers_per_projection[i], projection_name = fibers_to_axons(combined_h5_file, pathway_mat_files[i], axon_model, diam_fib[i], axon_length[i], centering_coordinates)
        print("Projection", projection_name, "seeded with", n_Ranviers_per_projection[i], "nodes of Ranvier\n")

    # just an example of parameters that will be filled in
    axon_dict = {
        'Axon_Model_Type': 'No_model',
        'Neuron_model_array_prepared': 1,  # now it is known
        'Name_prepared_neuron_array': '0',
        'diam_fib': [10],
        'n_Ranvier': [17],
    }

    if len(diam_fib) == 1 and not os.path.exists(combined_h5_file + '.h5'):  # stupid way
        n_Ranviers_per_projection_list = [n_Ranviers_per_projection[0]]
    elif len(diam_fib) == 1 and os.path.exists(combined_h5_file + '.h5'):  # stupid way
        n_Ranviers_per_projection_list = [int(n_Ranviers_per_projection + 0)]  # stupid trick
    else:
        n_Ranviers_per_projection_list = list(n_Ranviers_per_projection)

    diams_fib_true = []
    n_Ranviers_per_projection_true = []
    for i in range(len(diam_fib)):
        if n_Ranviers_per_projection[i] != 0:
            diams_fib_true.append(float(diam_fib[i]))
            n_Ranviers_per_projection_true.append(int(n_Ranviers_per_projection[i]))

    axon_dict['n_Ranvier'] = n_Ranviers_per_projection_true
    axon_dict['diam_fib'] = diams_fib_true
    axon_dict['Axon_Model_Type'] = axon_model
    axon_dict['Name_prepared_neuron_array'] = combined_h5_file + '.h5'

    name_of_directory = combined_h5_file.rsplit('/', 1)[0]

    import json
    with open(name_of_directory + '/Allocated_axons_parameters.json', 'w') as save_as_dict:
        json.dump(axon_dict, save_as_dict)