"""
@author: Konstantin Butenko

get_callmap() is obsolete and should be substituted in the code for the more general get_cellmap_tensors()
The former iterates over the mesh and marks FEM elements with tissue indices according to the segmented MRI data
and the electrode parameters. The latter does the same, but also maps the DTI data onto the elements.

IMPORTANT: the dielectric properties are not assigned here, only the tissue/material indices.
"""

import logging

logging.getLogger('UFL').setLevel(logging.WARNING)
logging.getLogger('FFC').setLevel(logging.WARNING)

from dolfin import MeshFunction, cells, HDF5File, File, set_log_active
import numpy as np
import time
import os

set_log_active(False)  # turns off debugging info

# this function is obsolete and should be rotated out (use get_cellmap_tensors() instead)
def get_cellmap(mesh, subdomains_assigned, Domains, MRI_param, default_material):

    """ mesh and subdomains_assigned are .xml files of the mesh and the corresponding mesh_physical_region.xml
        class instance Domain is defined in CAD_Salome.py
        class instance MRI_param is defined in MRI_DTI_processing.py
        default_material: index of the tissue assigned as default (1-CSF, 2-WM, 3-GM)

        Returns: a cell function that contacts material(tissue) indices. """

    Tissue_array = np.load(os.environ['PATIENTDIR'] + '/MRI_DTI_derived_data/Tissue_array_MRI.npy')

    subdomains = MeshFunction('size_t', mesh, 3)
    subdomains.set_all(0)

    for cell in cells(mesh):  # iterate over midpoints of cells and find the corresponding voxels of segmented tissue

        x_coord, y_coord, z_coord = (cell.midpoint().x(), cell.midpoint().y(), cell.midpoint().z())

        # find indices of these coordinates in the 3D nifti file
        coords = np.array([x_coord, y_coord, z_coord, 1.0]) - np.array(
            [MRI_param.MRI_shift[0], MRI_param.MRI_shift[1], MRI_param.MRI_shift[2], 0.0])  # we need to shift back to MRI coordinates
        basis_coords = np.dot(np.linalg.inv(MRI_param.affine_MRI), coords)  # map back to the 'orthonormal space' to compute the indices
        x_ind_basis, y_ind_basis, z_ind_basis, one = basis_coords.astype(int)[:]

        xv_mri = x_ind_basis  # defines number of steps to get to the voxels containing x[0] coordinate
        yv_mri = y_ind_basis * MRI_param.N_voxels[0]  # defines number of steps to get to the voxels containing x[0] and x[1] coordinates
        zv_mri = z_ind_basis * MRI_param.N_voxels[0] * MRI_param.N_voxels[1]  # defines number of steps to get to the voxels containing x[0], x[1] and x[2] coordinates
        k_mri = int(xv_mri + yv_mri + zv_mri)  # 1-D (flat) index

        # if x_ind_basis >= MRI_param.N_voxels[0] or y_ind_basis >= MRI_param.N_voxels[1] or z_ind_basis >= MRI_param.N_voxels[2] or x_coord < 0.0 or y_coord < 0.0 or z_coord < 0.0:
        #     print("Detected a cell outside of the segmented MRI")

        if k_mri >= MRI_param.N_voxels[0] * MRI_param.N_voxels[1] * MRI_param.N_voxels[2] or k_mri < 0:
            subdomains[cell] = default_material
            # skip the cell
            continue

        if int(Tissue_array[k_mri]) == 0:
            subdomains[cell] = default_material
        else:
            subdomains[cell] = int(Tissue_array[k_mri])   # subdomains have the same ordering as Tissue array [default, csf, wm, gm] + [encap, float]

    # reassign cells that are in the floating contacts
    if Domains.Float_contacts != -1:
        if isinstance(Domains.Float_contacts, int):
            subdomains.array()[subdomains_assigned.array() == Domains.Float_contacts] = 5
        else:
            for i in range(len(Domains.Float_contacts)):
                subdomains.array()[subdomains_assigned.array() == Domains.Float_contacts[i]] = 5  # 5 is index for float contact (very high cond and perm)

    # reassign cells that are in the encapsulation layer
    for i in range(len(Domains.Encup_index)):
        subdomains.array()[subdomains_assigned.array() == Domains.Encup_index[i]] = 4  # 4 is index of encap

    del Tissue_array

    return subdomains


def get_cellmap_tensors(mesh, subdomains_assigned, Domains, MRI_param, DTI_param, default_material):

    """ mesh and subdomains_assigned are .xml files of the mesh and the corresponding mesh_physical_region.xml
        class instance Domain is defined in CAD_Salome.py
        class instances MRI_param and DTI_param are defined in MRI_DTI_processing.py
        default_material: index of the tissue assigned as default (1-CSF, 2-WM, 3-GM)

        Returns: a cell function that contacts material(tissue) indices.
        If tensor data were provided, also stores their mapping on the mesh
         in /Results_adaptive/Tensors_to_solve_num_el_ """

    Tissue_array = np.load(os.environ['PATIENTDIR'] + '/MRI_DTI_derived_data/Tissue_array_MRI.npy')

    if DTI_param != 0:
        DTI_array = np.load(os.environ['PATIENTDIR'] + '/MRI_DTI_derived_data/Tensor_array_DTI.npy')
        # Create mesh functions for c00, c01, c11
        c00 = MeshFunction("double", mesh, 3, 1.0)  #  FSL tensor ordering (xx,xy,xz,yy,yz,zz)
        c01 = MeshFunction("double", mesh, 3, 0.0)
        c02 = MeshFunction("double", mesh, 3, 0.0)
        c11 = MeshFunction("double", mesh, 3, 1.0)
        c12 = MeshFunction("double", mesh, 3, 0.0)
        c22 = MeshFunction("double", mesh, 3, 1.0)

    subdomains = MeshFunction('size_t', mesh, 3)  #  will contain indices of the tissue for each FEM element
    subdomains.set_all(0)

    for cell in cells(mesh): # iterate over midpoints of cells and find the corresponding voxels of segmented tissue

        x_coord, y_coord, z_coord = (cell.midpoint().x(), cell.midpoint().y(), cell.midpoint().z())

        # find indices of these coordinates in the 3D nifti file
        coords = np.array([x_coord, y_coord, z_coord, 1.0]) - np.array(
            [MRI_param.MRI_shift[0], MRI_param.MRI_shift[1], MRI_param.MRI_shift[2], 0.0])  # we need to shift back to MRI coordinates
        basis_coords = np.dot(np.linalg.inv(MRI_param.affine_MRI), coords)  # map back to the 'orthonormal space' to compute the indices
        x_ind_basis, y_ind_basis, z_ind_basis, one = basis_coords.astype(int)[:]

        xv_mri = x_ind_basis  # defines number of steps to get to the voxels containing x[0] coordinate
        yv_mri = y_ind_basis * MRI_param.N_voxels[0]  # defines number of steps to get to the voxels containing x[0] and x[1] coordinates
        zv_mri = z_ind_basis * MRI_param.N_voxels[0] * MRI_param.N_voxels[1]  # defines number of steps to get to the voxels containing x[0], x[1] and x[2] coordinates
        k_mri = int(xv_mri + yv_mri + zv_mri)  # 1-D (flat) index

        # if x_ind_basis >= MRI_param.N_voxels[0] or y_ind_basis >= MRI_param.N_voxels[1] or z_ind_basis >= MRI_param.N_voxels[2] or x_coord < 0.0 or y_coord < 0.0 or z_coord < 0.0:
        #     print("Detected a cell outside of the segmented MRI")

        # check if the flat index is outside of the bounds
        if k_mri >= MRI_param.N_voxels[0] * MRI_param.N_voxels[1] * MRI_param.N_voxels[2] or k_mri < 0:
            subdomains[cell] = default_material
            # skip the cell, it will have the default material assigned with isotropic properties
            continue

        if x_ind_basis >= MRI_param.N_voxels[0] or y_ind_basis >= MRI_param.N_voxels[1] or z_ind_basis >= MRI_param.N_voxels[2] or x_coord < 0.0 or y_coord < 0.0 or z_coord < 0.0:
            subdomains[cell] = default_material
            # skip the cell, it will have the default material assigned with isotropic properties
            continue

        # assign the tissue indices
        if int(Tissue_array[k_mri]) == 0:
            subdomains[cell] = default_material
        else:
            subdomains[cell] = int(Tissue_array[k_mri])   # subdomains have the same ordering as Tissue array [default, csf, wm, gm] + [encap, float]

        # map tensors to elements if DTI data was provided (and processed)
        if DTI_param != 0:
            basis_coords_DTI = np.dot(np.linalg.inv(DTI_param.affine_DTI), coords)  # map back to the 'orthonormal space' to compute the indices
            x_ind_basis_DTI, y_ind_basis_DTI, z_ind_basis_DTI, one = basis_coords_DTI.astype(int)[:]

            xv_dti = x_ind_basis_DTI  # defines number of steps to get to the voxels containing x[0] coordinate
            yv_dti = y_ind_basis_DTI * DTI_param.N_voxels[0]  # defines number of steps to get to the voxels containing x[0] and x[1] coordinates
            zv_dti = z_ind_basis_DTI * DTI_param.N_voxels[0] * DTI_param.N_voxels[1]  # defines number of steps to get to the voxels containing x[0], x[1] and x[2] coordinates
            k_dti = int(xv_dti + yv_dti + zv_dti)

            if k_dti >= DTI_param.N_voxels[0] * DTI_param.N_voxels[1] * DTI_param.N_voxels[2] or k_dti < 0 or x_ind_basis_DTI >= DTI_param.N_voxels[0] or y_ind_basis_DTI >= DTI_param.N_voxels[1] or z_ind_basis_DTI >= DTI_param.N_voxels[2] or x_ind_basis_DTI < 0 or y_ind_basis_DTI < 0 or z_ind_basis_DTI < 0:
                # cell is outside DTI data, assign first voxel of DTI, which will not fullfil the following condition, and the cell will be considered as isotropic
                continue
            else:
                c00[cell] = DTI_array[k_dti, 0]  # we don't scale the tensor with cond. values here, because they depend on frequency
                c01[cell] = DTI_array[k_dti, 1]
                c02[cell] = DTI_array[k_dti, 2]
                c11[cell] = DTI_array[k_dti, 3]
                c12[cell] = DTI_array[k_dti, 4]
                c22[cell] = DTI_array[k_dti, 5]

    # reassign cells that are in the encapsulation layer (isotropic conductivity!)
    for i in range(len(Domains.Encup_index)):
        subdomains.array()[subdomains_assigned.array() == Domains.Encup_index[i]] = 4
        if DTI_param != 0:
            c00.array()[subdomains_assigned.array() == Domains.Encup_index[i]] = 1.0
            c01.array()[subdomains_assigned.array() == Domains.Encup_index[i]] = 0.0
            c02.array()[subdomains_assigned.array() == Domains.Encup_index[i]] = 0.0
            c11.array()[subdomains_assigned.array() == Domains.Encup_index[i]] = 1.0
            c12.array()[subdomains_assigned.array() == Domains.Encup_index[i]] = 0.0
            c22.array()[subdomains_assigned.array() == Domains.Encup_index[i]] = 1.0

    # reassign cells that are in the floating contacts (isotropic conductivity!)
    if Domains.Float_contacts != -1:
        if isinstance(Domains.Float_contacts, int):
            subdomains.array()[subdomains_assigned.array() == Domains.Float_contacts] = 5
            if DTI_param != 0:
                c00.array()[subdomains_assigned.array() == Domains.Float_contacts] = 1.0
                c01.array()[subdomains_assigned.array() == Domains.Float_contacts] = 0.0
                c02.array()[subdomains_assigned.array() == Domains.Float_contacts] = 0.0
                c11.array()[subdomains_assigned.array() == Domains.Float_contacts] = 1.0
                c12.array()[subdomains_assigned.array() == Domains.Float_contacts] = 0.0
                c22.array()[subdomains_assigned.array() == Domains.Float_contacts] = 1.0
        else:
            for i in range(len(Domains.Float_contacts)):
                subdomains.array()[subdomains_assigned.array() == Domains.Float_contacts[i]] = 5
                if DTI_param != 0:
                    c00.array()[subdomains_assigned.array() == Domains.Float_contacts[i]] = 1.0
                    c01.array()[subdomains_assigned.array() == Domains.Float_contacts[i]] = 0.0
                    c02.array()[subdomains_assigned.array() == Domains.Float_contacts[i]] = 0.0
                    c11.array()[subdomains_assigned.array() == Domains.Float_contacts[i]] = 1.0
                    c12.array()[subdomains_assigned.array() == Domains.Float_contacts[i]] = 0.0
                    c22.array()[subdomains_assigned.array() == Domains.Float_contacts[i]] = 1.0

    if DTI_param != 0:
        hdf = HDF5File(mesh.mpi_comm(), os.environ['PATIENTDIR'] + '/Results_adaptive/Tensors_to_solve_num_el_' + str(
            mesh.num_cells()) + '.h5', 'w')
        hdf.write(c00, "/c00")
        hdf.write(c01, "/c01")
        hdf.write(c02, "/c02")
        hdf.write(c11, "/c11")
        hdf.write(c12, "/c12")
        hdf.write(c22, "/c22")
        hdf.close()

        file = File(os.environ['PATIENTDIR'] + '/Tensors/c00_unscaled.pvd')
        file << c00, mesh
        file = File(os.environ['PATIENTDIR'] + '/Tensors/c01_unscaled.pvd')
        file << c01, mesh
        file = File(os.environ['PATIENTDIR'] + '/Tensors/c02_unscaled.pvd')
        file << c02, mesh
        file = File(os.environ['PATIENTDIR'] + '/Tensors/c11_unscaled.pvd')
        file << c11, mesh
        file = File(os.environ['PATIENTDIR'] + '/Tensors/c12_unscaled.pvd')
        file << c12, mesh
        file = File(os.environ['PATIENTDIR'] + '/Tensors/c22_unscaled.pvd')
        file << c22, mesh

        # Store to file
        # mesh_file = File("mesh.xml.gz")
        # c00_file = File("Tensors/c00.xml.gz")
        # c01_file = File("Tensors/c01.xml.gz")
        # c02_file = File("Tensors/c02.xml.gz")
        # c11_file = File("Tensors/c11.xml.gz")
        # c12_file = File("Tensors/c12.xml.gz")
        # c22_file = File("Tensors/c22.xml.gz")

        ##mesh_file << mesh
        # c00_file << c00
        # c01_file << c01
        # c02_file << c02
        # c11_file << c11
        # c12_file << c12
        # c22_file << c22

    return subdomains
