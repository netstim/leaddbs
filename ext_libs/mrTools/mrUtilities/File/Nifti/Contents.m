% CBI - NIFTI tools for Matlab
%
% (c) Jonas Larsson 2005
%
% NOTE: these tools currently do not support compressed Nifti files. 
%
% TOOLS FOR READING NIFTI FILES
% =============================
% cbiReadNiftiHeader.m
% - Reads in the header of a Nifti-1 file
% cbiReadNifti.m
% - Reads in the data and header of a Nifti-1 file
% 
% TOOLS FOR WRITING NIFTI FILES
% =============================
% cbiWriteNiftiHeader.m
% - Writes a Nifti-1 header
% cbiWriteNifti.m
% - Writes a Nifti-1 file (data and header)
%
% TOOLS FOR CONVERTING NIFTI FILES
% ================================
% cbiSwapNiftiDimensions.m
% - Swaps the dimensions of a Nifti data set or file. Analogous to FSL avwswapdim, but correctly sets qform & sform data.
%
% NIFTI UTILITIES
% ===============
% cbiCreateNiftiHeader.m 
% - Creates a new header and/or checks an existing header for consistency and valid fields
% cbiHomogeneousToQuaternion.m 
% - Converts a 4x4 or 3x3 matrix to quaternions
% cbiQuaternionToHomogeneous.m
% - Converts quaternions and qoffset data into a 4x4 homogeneous matrix
% cbiSetNiftiQform.m
% - Sets the qform and quaternion information of a Nifti header
% cbiSetNiftiSform.m
% - Sets the sform and srow information of a Nifti header
% cbiMatlabDatatype2Nifti.m
% - Converts matlab data type in string format to a Nifti-1 integer code; see nifti1.h for details 
% cbiNiftiDatatype2Matlab.m
% - Converts a Nifti-1 integer data tyoe code into matlab data type in string format
% cbiSizeofNifti.m
% - Returns the size in bytes per data point of a given Nifti-1 format.
% nifti1.h
% - Header file of the nifti-1 library written by Robert Cox; describes the Nifti-1 format
