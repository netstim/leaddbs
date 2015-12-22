function impexp_NiftiMrStruct = impexp_cfg_NiftiMrStruct
%_______________________________________________________________________
%
% This toolbox contains various helper functions to convert NifTI images 
% to and from mrStructs.
%
% This toolbox is free but copyright software, distributed under the 
% terms of the GNU General Public Licence as published by the Free 
% Software Foundation (either version 2, as given in file 
% spm_LICENCE.man, or at your option, any later version). Further details 
% on "copyleft" can be found at http://www.gnu.org/copyleft/.
% The toolbox consists of the files listed in its Contents.m file.
%_______________________________________________________________________
%
% @(#) $Id: impexp_cfg_NiftiMrStruct.m,v 1.14 2014/11/03 15:33:05 kellnere Exp $


rev='$Revision: 1.14 $';

% MATLABBATCH Configuration file for toolbox 'MedPhysConvert'
% This code has been automatically generated.
% ---------------------------------------------------------------------
% srcimgs Source Images
% ---------------------------------------------------------------------
srcimgs         = cfg_files;
srcimgs.tag     = 'srcimgs';
srcimgs.name    = 'Source Images';
srcimgs.help    = {'Specify here files to be included in the mrStruct. The number of files must match the expected values: 1 for image/volume, any number for series2D/3D and image/volumeEchos, #series-by-#echos for series2D/3DEchos. The decision on whether 2D or 3D data are imported will be made based on the NifTI header of the images. All images need to have same orientation and voxel sizes.'};
srcimgs.filter = 'image';
srcimgs.ufilter = '.*';
srcimgs.num     = [1 Inf];
% ---------------------------------------------------------------------
% volume image/volume
% ---------------------------------------------------------------------
volume         = cfg_const;
volume.tag     = 'volume';
volume.name    = 'image/volume';
volume.val = {'volume'};
volume.help    = {''};
% ---------------------------------------------------------------------
% series3D series2D/3D
% ---------------------------------------------------------------------
series3D         = cfg_const;
series3D.tag     = 'series3D';
series3D.name    = 'series2D/3D';
series3D.val = {'series3D'};
series3D.help    = {''};
% ---------------------------------------------------------------------
% volumeEchos image/volumeEchos
% ---------------------------------------------------------------------
volumeEchos         = cfg_const;
volumeEchos.tag     = 'volumeEchos';
volumeEchos.name    = 'image/volumeEchos';
volumeEchos.val = {'volumeEchos'};
volumeEchos.help    = {''};
% ---------------------------------------------------------------------
% inorder Input Order
% ---------------------------------------------------------------------
inorder         = cfg_menu;
inorder.tag     = 'inorder';
inorder.name    = 'Input Order';
inorder.help    = {'Specify the order of your images with respect to echos and series.'};
inorder.labels = {
                  'series1 all echos | series2 all echos ... seriesX all echos'
                  'echos1 all series | echos2 all series ... echosX all series'
}';
inorder.values = {
                  'series'
                  'echos'
}';
% ---------------------------------------------------------------------
% nseries #images in series
% ---------------------------------------------------------------------
nseries         = cfg_entry;
nseries.tag     = 'nseries';
nseries.name    = '#images in series';
nseries.help    = {'The number of images must be equal to #images-by-#echos.'};
nseries.strtype = 'n';
nseries.num     = [1 1];
% ---------------------------------------------------------------------
% nechos #echos
% ---------------------------------------------------------------------
nechos         = cfg_entry;
nechos.tag     = 'nechos';
nechos.name    = '#echos';
nechos.help    = {'The number of images must be equal to #images-by-#echos.'};
nechos.strtype = 'n';
nechos.num     = [1 1];
% ---------------------------------------------------------------------
% series3DEchos series2D/3DEchos
% ---------------------------------------------------------------------
series3DEchos         = cfg_branch;
series3DEchos.tag     = 'series3DEchos';
series3DEchos.name    = 'series2D/3DEchos';
series3DEchos.val     = {inorder nseries nechos };
series3DEchos.help    = {'Specify image order, #echos and #images.'};
% ---------------------------------------------------------------------
% mrtype mrStruct type
% ---------------------------------------------------------------------
mrtype         = cfg_choice;
mrtype.tag     = 'mrtype';
mrtype.name    = 'mrStruct type';
mrtype.help    = {'Specify the appropriate type of mrStruct. The decision whether the data is 2D or 3D will be made based on the image header information.'};
mrtype.values  = {volume series3D volumeEchos series3DEchos };
% ---------------------------------------------------------------------
% outdir Output directory
% ---------------------------------------------------------------------
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output directory';
outdir.help    = {'Files produced by this function will be written into this output directory.'};
outdir.filter = 'dir';
outdir.ufilter = '.*';
outdir.num     = [1 1];
% ---------------------------------------------------------------------
% fname Output Filename
% ---------------------------------------------------------------------
fname         = cfg_entry;
fname.tag     = 'fname';
fname.name    = 'Output Filename';
fname.help    = {'The output file specified here is written to the selected output directory. If the input mrStruct contains multiple volumes (series, echos), then the filename part will be extended with running indices.'};
fname.strtype = 's';
fname.num     = [1 Inf];
% ---------------------------------------------------------------------
% outimg Output File & Directory
% ---------------------------------------------------------------------
outimg         = cfg_branch;
outimg.tag     = 'outimg';
outimg.name    = 'Output File & Directory';
outimg.val     = {outdir fname };
outimg.help    = {'Specify a output filename and target directory.'};
% ---------------------------------------------------------------------
% outvar Output variable
% ---------------------------------------------------------------------
outvar         = cfg_const;
outvar.tag     = 'outvar';
outvar.name    = 'Output variable';
outvar.val{1}  = true;
% ---------------------------------------------------------------------
% outchoice Output destination
% ---------------------------------------------------------------------
outchoice         = cfg_choice;
outchoice.tag     = 'outchoice';
outchoice.name    = 'Output destination';
outchoice.help    = {'Output can be saved to disk or into a MATLAB variable in the ''base'' workspace. No check is performed whether the specified file or variable already exists and any previous contents will be overwritten.'};
outchoice.values  = {outimg outvar};
% ---------------------------------------------------------------------
% nifti2mrstruct Convert Nifti images to mrStruct
% ---------------------------------------------------------------------
nifti2mrstruct         = cfg_exbranch;
nifti2mrstruct.tag     = 'nifti2mrstruct';
nifti2mrstruct.name    = 'Convert Nifti images to mrStruct';
nifti2mrstruct.val     = {srcimgs mrtype outchoice };
nifti2mrstruct.help    = {
                        'Import NifTI images into mrstruct'
                        'FORMAT [mrStruct, errStr] = nifti_to_mrstruct(mrType, fNameCell)'
                        '======'
                        'Input arguments'
                        'mrType    - a valid mrstruct type (all image, volume, series types are supported, but no spectral types)'
                        'fNameCell - cell array of filenames. The meaning of cell array dimensions depends on the mrType:'
                        '''image'',''volume''      - a 1x1 cell array'
                        '''(image|volume)Echos'' - #Echos-by-1 cell array'
                        '''series2D'',''series3D'' - #Timepoints-by-1 cell array'
                        '''(series2D|series3D)Echos'' - #Echos-by-#Timepoints cell array'
                        'Output arguments'
                        'mrStruct - mrStruct with all data and position information'
                        'errStr   - empty, if successful. Otherwise error description.'
                        'All files that have to be imported into a single mrstruct need to have identical voxel sizes, image dimensions and slice orientation. If this is not the case, they need to be resliced within SPM before they can be imported.'
                        'Data are rearranged from NifTI array order to mrStruct array order. A coordinate transformation matrix from mrStruct voxel coordinates to mm coordinates is computed based on the NifTI transformation matrix.'
                        '_______________________________________________________________________'
                        'Bjoern W. Kreher'
                        '08/05'
                        ''
                        'UNIX'
                        '_______________________________________________________________________'
                        'nifti_to_mrstruct.m,v 1.1 2006/07/18 13:49:39 bkreher Exp'
}';
nifti2mrstruct.prog = @(job)impexp_run_nifti2mrstruct('run',job);
nifti2mrstruct.vout = @(job)impexp_run_nifti2mrstruct('vout',job);
% ---------------------------------------------------------------------
% srcstruct File containing mrStruct
% ---------------------------------------------------------------------
srcstruct         = cfg_files;
srcstruct.tag     = 'srcstruct';
srcstruct.name    = 'File containing mrStruct';
srcstruct.help    = {'Select a .mat file containing the mrStruct.'};
srcstruct.filter = 'mat';
srcstruct.ufilter = '.*';
srcstruct.num     = [1 Inf];
% ---------------------------------------------------------------------
% srcvar Directly passed variable
% ---------------------------------------------------------------------
srcvar         = cfg_entry;
srcvar.tag     = 'srcvar';
srcvar.name    = 'mrStruct';
srcvar.help    = {'Specify a mrStruct variable.'};
srcvar.strtype = 'e';
srcvar.num     = [1 Inf];
srcvar.check   = @(job)impexp_run_mrstruct2nifti('check','ismrstruct',job);
% ---------------------------------------------------------------------
% srcchoice Input source
% ---------------------------------------------------------------------
srcchoice         = cfg_choice;
srcchoice.tag     = 'srcchoice';
srcchoice.name    = 'Input source';
srcchoice.help    = {'Input mrStructs can be either loaded from disk or passed as a variable.'};
srcchoice.values  = {srcstruct srcvar};
% ---------------------------------------------------------------------
% outdir Output directory
% ---------------------------------------------------------------------
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output directory';
outdir.help    = {'Files produced by this function will be written into this output directory.'};
outdir.filter = 'dir';
outdir.ufilter = '.*';
outdir.num     = [1 1];
% ---------------------------------------------------------------------
% fname Output Filename
% ---------------------------------------------------------------------
fname         = cfg_entry;
fname.tag     = 'fname';
fname.name    = 'Output Filename';
fname.help    = {'The output file specified here is written to the selected output directory. '...
    'If there is more than one input mrStruct or an mrStruct contains '...
    'multiple volumes (series, echos), then the filename part will be extended with running indices.'};
fname.strtype = 's';
fname.num     = [1 Inf];
% ---------------------------------------------------------------------
% outimg Output File & Directory
% ---------------------------------------------------------------------
outimg         = cfg_branch;
outimg.tag     = 'outimg';
outimg.name    = 'Specify Output File & Directory';
outimg.val     = {outdir fname};
outimg.help    = {'Specify a output filename and target directory.'};
% ---------------------------------------------------------------------
% autoimg Automatically generate Output Filename
% ---------------------------------------------------------------------
autoimg        = cfg_const;
autoimg.tag    = 'autoimg';
autoimg.name   = 'Automatically generate Output Filename';
autoimg.val    = {true};
autoimg.help   = {'Generate image filenames automatically. If the sources '...
    'are mrStruct files, then their filenames will be used as basename. '...
    'Otherwise, files will be saved to the current directory using running '...
    'indices.'};
% ---------------------------------------------------------------------
% outname Output Naming Scheme
% ---------------------------------------------------------------------
outname         = cfg_choice;
outname.tag     = 'outname';
outname.name    = 'Output Naming Scheme';
outname.values  = {outimg autoimg};
outname.help    = {'Specify a output filename and target directory or use some heuristics.'};
% ---------------------------------------------------------------------
% dtype Data Type
% ---------------------------------------------------------------------
dtype         = cfg_menu;
dtype.tag     = 'dtype';
dtype.name    = 'Data Type';
dtype.help    = {'Data type of output images.'};
dtype.labels = cellstr(spm_type(spm_type));
dtype.values = num2cell(spm_type);
% ---------------------------------------------------------------------
% output Output Options
% ---------------------------------------------------------------------
output         = cfg_branch;
output.tag     = 'output';
output.name    = 'Output Options';
output.val     = {outname dtype};
output.help    = {'Specify a output filename and data type.'};
% ---------------------------------------------------------------------
% mrstruct2nifti Convert mrStruct to Nifti image(s)
% ---------------------------------------------------------------------
mrstruct2nifti         = cfg_exbranch;
mrstruct2nifti.tag     = 'mrstruct2nifti';
mrstruct2nifti.name    = 'Convert mrStruct to Nifti image(s)';
mrstruct2nifti.val     = {srcchoice output };
mrstruct2nifti.help    = {
                        'Convert mrStruct image/volume/series to NifTI compatible files'
                        'FORMAT [res, errStr]= mrstruct_to_nifti(mrStruct, fName)'
                        '======'
                        'This routine tries to convert mrStruct data to NifTI compatible files using SPM5 functions spm_type, spm_platform, spm_write_vol. Image data will be saved in user given format (dataType). Default format is int16. If available, a NifTI compatible coordinate transformation matrix will be'
                        'computed from mrStruct.edges. If this is not possible, axial orientation'
                        'is assumed and voxel sizes are read from mrStruct.vox. If no voxel sizes'
                        'are given, a default of 1x1x1mm will be assumed. Image data are'
                        'rearranged in Nifti/Analyze format and slices are re-ordered depending on'
                        'the handedness of the voxel-to-world coordinate mapping.'
                        ''
                        'Input arguments'
                        'mrStruct - mrstruct to be converted'
                        'fName    - location and filename for created file(s). If mrStruct contains series or multi-echo data, separate volumes will be created for each echo or series image/volume by appending a running series/echo number to the file name.'
                        'dataType - data type for output (''uint8'',''int16'',''int32'',''float32'','
                        '''float64'',''int8'',''uint16'',''uint32'')'
                        'Output arguments'
                        'res      - a cell array containing volume handles (see spm_vol). The array is shaped according to the series/echo dimensions of the mrStruct data array.'
                        'errStr   - empty, if successful operation. Otherwise, an error message.'
                        '_______________________________________________________________________'
                        'Bjoern W. Kreher'
                        '08/05'
                        ''
                        'UNIX'
                        '_______________________________________________________________________'
                        '$Id: impexp_cfg_NiftiMrStruct.m,v 1.14 2014/11/03 15:33:05 kellnere Exp $'
}';
mrstruct2nifti.prog = @(job)impexp_run_mrstruct2nifti('run',job);
mrstruct2nifti.vout = @(job)impexp_run_mrstruct2nifti('vout',job);

% ---------------------------------------------------------------------
% srcstruct File containing mrStruct
% ---------------------------------------------------------------------
srcdtdstruct         = cfg_files;
srcdtdstruct.tag     = 'srcdtdstruct';
srcdtdstruct.name    = 'File containing tensors and B0-image';
srcdtdstruct.help    = {'Select a _DTD.mat file containing the tensor data.'};
srcdtdstruct.filter = 'mat';
srcdtdstruct.ufilter = '_DTD.mat';
srcdtdstruct.num     = [1 Inf];
% ---------------------------------------------------------------------
% srcvar Directly passed variable
% ---------------------------------------------------------------------
srcdtdvar         = cfg_entry;
srcdtdvar.tag     = 'srcdtdvar';
srcdtdvar.name    = 'dtdStruct';
srcdtdvar.help    = {'Specify a dtdStruct variable.'};
srcdtdvar.strtype = 'e';
srcdtdvar.num     = [1 Inf];
srcdtdvar.check   = @(job)impexp_run_bo2nifti('check','isdtdstruct',job);
% ---------------------------------------------------------------------
% srcchoice Input source
% ---------------------------------------------------------------------
srcdtdchoice         = cfg_choice;
srcdtdchoice.tag     = 'srcdtdchoice';
srcdtdchoice.name    = 'Input source';
srcdtdchoice.help    = {'Input dtdStructs can be either loaded from disk or passed as a variable.'};
srcdtdchoice.values  = {srcdtdstruct srcdtdvar};
% ---------------------------------------------------------------------
% outdir Output directory
% ---------------------------------------------------------------------
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output directory';
outdir.help    = {'Files produced by this function will be written into this output directory.'};
outdir.filter = 'dir';
outdir.ufilter = '.*';
outdir.num     = [1 1];
% ---------------------------------------------------------------------
% fname Output Filename
% ---------------------------------------------------------------------
fname         = cfg_entry;
fname.tag     = 'fname';
fname.name    = 'Output Filename';
fname.help    = {'The output file specified here is written to the selected output directory. '...
    'If there is more than one input mrStruct or an mrStruct contains '...
    'multiple volumes (series, echos), then the filename part will be extended with running indices.'};
fname.strtype = 's';
fname.num     = [1 Inf];
% ---------------------------------------------------------------------
% outimg Output File & Directory
% ---------------------------------------------------------------------
outimg         = cfg_branch;
outimg.tag     = 'outimg';
outimg.name    = 'Specify Output File & Directory';
outimg.val     = {outdir fname};
outimg.help    = {'Specify a output filename and target directory.'};
% ---------------------------------------------------------------------
% autoimg Automatically generate Output Filename
% ---------------------------------------------------------------------
autoimg        = cfg_const;
autoimg.tag    = 'autoimg';
autoimg.name   = 'Automatically generate Output Filename';
autoimg.val    = {true};
autoimg.help   = {'Generate image filenames automatically. If the sources '...
    'are mrStruct files, then their filenames will be used as basename. '...
    'Otherwise, files will be saved to the current directory using running '...
    'indices.'};
% ---------------------------------------------------------------------
% outname Output Naming Scheme
% ---------------------------------------------------------------------
outname         = cfg_choice;
outname.tag     = 'outname';
outname.name    = 'Output Naming Scheme';
outname.values  = {outimg autoimg};
outname.help    = {'Specify a output filename and target directory or use some heuristics.'};

% ---------------------------------------------------------------------
% eigVal3_2_nifti Extract eigVal1 map and convert mrStruct to Nifti image
% ---------------------------------------------------------------------
eigVal3_2_nifti         = cfg_exbranch;
eigVal3_2_nifti.tag     = 'eigVal3_2_nifti';
eigVal3_2_nifti.name    = 'Eigval 3 map (DTD) to Nifti';
eigVal3_2_nifti.val     = {srcdtdchoice outname};
eigVal3_2_nifti.help    = {
                        'Extract the eigenvalue 3 from the dtdstruct and convert this'
                        'mrStruct volume to NifTI compatible files'
                        '======'
                        'This routine tries to convert mrStruct data to NifTI compatible files using SPM5 functions spm_type, spm_platform, spm_write_vol. Image data will be saved in user given format (dataType). Default format is int16. If available, a NifTI compatible coordinate transformation matrix will be'
                        'computed from mrStruct.edges. If this is not possible, axial orientation'
                        'is assumed and voxel sizes are read from mrStruct.vox. If no voxel sizes'
                        'are given, a default of 1x1x1mm will be assumed. Image data are'
                        'rearranged in Nifti/Analyze format and slices are re-ordered depending on'
                        'the handedness of the voxel-to-world coordinate mapping.'
                        ''
}';
eigVal3_2_nifti.prog = @(job)impexp_run_eigVal3_2_nifti('run',job);
eigVal3_2_nifti.vout = @(job)impexp_run_eigVal3_2_nifti('vout',job);

% ---------------------------------------------------------------------
% eigVal2_2_nifti Extract eigVal1 map and convert mrStruct to Nifti image
% ---------------------------------------------------------------------
eigVal2_2_nifti         = cfg_exbranch;
eigVal2_2_nifti.tag     = 'eigVal2_2_nifti';
eigVal2_2_nifti.name    = 'Eigval 2 map (DTD) to Nifti';
eigVal2_2_nifti.val     = {srcdtdchoice outname};
eigVal2_2_nifti.help    = {
                        'Extract the eigenvalue 2 from the dtdstruct and convert this'
                        'mrStruct volume to NifTI compatible files'
                        '======'
                        'This routine tries to convert mrStruct data to NifTI compatible files using SPM5 functions spm_type, spm_platform, spm_write_vol. Image data will be saved in user given format (dataType). Default format is int16. If available, a NifTI compatible coordinate transformation matrix will be'
                        'computed from mrStruct.edges. If this is not possible, axial orientation'
                        'is assumed and voxel sizes are read from mrStruct.vox. If no voxel sizes'
                        'are given, a default of 1x1x1mm will be assumed. Image data are'
                        'rearranged in Nifti/Analyze format and slices are re-ordered depending on'
                        'the handedness of the voxel-to-world coordinate mapping.'
                        ''
}';
eigVal2_2_nifti.prog = @(job)impexp_run_eigVal2_2_nifti('run',job);
eigVal2_2_nifti.vout = @(job)impexp_run_eigVal2_2_nifti('vout',job);

% ---------------------------------------------------------------------
% eigVal1_2_nifti Extract eigVal1 map and convert mrStruct to Nifti image
% ---------------------------------------------------------------------
eigVal1_2_nifti         = cfg_exbranch;
eigVal1_2_nifti.tag     = 'eigVal1_2_nifti';
eigVal1_2_nifti.name    = 'Eigval 1 map (DTD) to Nifti';
eigVal1_2_nifti.val     = {srcdtdchoice outname};
eigVal1_2_nifti.help    = {
                        'Extract the eigenvalue 1 from the dtdstruct and convert this'
                        'mrStruct volume to NifTI compatible files'
                        '======'
                        'This routine tries to convert mrStruct data to NifTI compatible files using SPM5 functions spm_type, spm_platform, spm_write_vol. Image data will be saved in user given format (dataType). Default format is int16. If available, a NifTI compatible coordinate transformation matrix will be'
                        'computed from mrStruct.edges. If this is not possible, axial orientation'
                        'is assumed and voxel sizes are read from mrStruct.vox. If no voxel sizes'
                        'are given, a default of 1x1x1mm will be assumed. Image data are'
                        'rearranged in Nifti/Analyze format and slices are re-ordered depending on'
                        'the handedness of the voxel-to-world coordinate mapping.'
                        ''
}';
eigVal1_2_nifti.prog = @(job)impexp_run_eigVal1_2_nifti('run',job);
eigVal1_2_nifti.vout = @(job)impexp_run_eigVal1_2_nifti('vout',job);

% ---------------------------------------------------------------------
% trace2nifti Extract trace map and convert mrStruct to Nifti image
% ---------------------------------------------------------------------
trace2nifti         = cfg_exbranch;
trace2nifti.tag     = 'trace2nifti';
trace2nifti.name    = 'Trace map (DTD) to Nifti';
trace2nifti.val     = {srcdtdchoice outname};
trace2nifti.help    = {
                        'Extract the Trace map from the dtdstruct and convert this'
                        'mrStruct volume to NifTI compatible files'
                        '======'
                        'This routine tries to convert mrStruct data to NifTI compatible files using SPM5 functions spm_type, spm_platform, spm_write_vol. Image data will be saved in user given format (dataType). Default format is int16. If available, a NifTI compatible coordinate transformation matrix will be'
                        'computed from mrStruct.edges. If this is not possible, axial orientation'
                        'is assumed and voxel sizes are read from mrStruct.vox. If no voxel sizes'
                        'are given, a default of 1x1x1mm will be assumed. Image data are'
                        'rearranged in Nifti/Analyze format and slices are re-ordered depending on'
                        'the handedness of the voxel-to-world coordinate mapping.'
                        ''
}';
trace2nifti.prog = @(job)impexp_run_trace2nifti('run',job);
trace2nifti.vout = @(job)impexp_run_trace2nifti('vout',job);

% ---------------------------------------------------------------------
% fa2nifti Extract FA map and convert mrStruct to Nifti image
% ---------------------------------------------------------------------
fa2nifti         = cfg_exbranch;
fa2nifti.tag     = 'fa2nifti';
fa2nifti.name    = 'FA map (DTD) to Nifti';
fa2nifti.val     = {srcdtdchoice outname};
fa2nifti.help    = {
                        'Extract the FA map from the dtdstruct and convert this'
                        'mrStruct volume to NifTI compatible files'
                        '======'
                        'This routine tries to convert mrStruct data to NifTI compatible files using SPM5 functions spm_type, spm_platform, spm_write_vol. Image data will be saved in user given format (dataType). Default format is int16. If available, a NifTI compatible coordinate transformation matrix will be'
                        'computed from mrStruct.edges. If this is not possible, axial orientation'
                        'is assumed and voxel sizes are read from mrStruct.vox. If no voxel sizes'
                        'are given, a default of 1x1x1mm will be assumed. Image data are'
                        'rearranged in Nifti/Analyze format and slices are re-ordered depending on'
                        'the handedness of the voxel-to-world coordinate mapping.'
                        ''
}';
fa2nifti.prog = @(job)impexp_run_fa2nifti('run',job);
fa2nifti.vout = @(job)impexp_run_fa2nifti('vout',job);

% ---------------------------------------------------------------------
% bo2nifti Extract b0-image and convert mrStruct to Nifti image
% ---------------------------------------------------------------------
bo2nifti         = cfg_exbranch;
bo2nifti.tag     = 'bo2nifti';
bo2nifti.name    = 'B0 image (DTD) to Nifti';
bo2nifti.val     = {srcdtdchoice outname};
bo2nifti.help    = {
                        'Extract the mean b0-image from the dtdstruct and convert this'
                        'mrStruct volume to NifTI compatible files'
                        '======'
                        'This routine tries to convert mrStruct data to NifTI compatible files using SPM5 functions spm_type, spm_platform, spm_write_vol. Image data will be saved in user given format (dataType). Default format is int16. If available, a NifTI compatible coordinate transformation matrix will be'
                        'computed from mrStruct.edges. If this is not possible, axial orientation'
                        'is assumed and voxel sizes are read from mrStruct.vox. If no voxel sizes'
                        'are given, a default of 1x1x1mm will be assumed. Image data are'
                        'rearranged in Nifti/Analyze format and slices are re-ordered depending on'
                        'the handedness of the voxel-to-world coordinate mapping.'
                        ''
                        '_______________________________________________________________________'
                        '$Id: impexp_cfg_NiftiMrStruct.m,v 1.14 2014/11/03 15:33:05 kellnere Exp $'
}';
bo2nifti.prog = @(job)impexp_run_bo2nifti('run',job);
bo2nifti.vout = @(job)impexp_run_bo2nifti('vout',job);

% ---------------------------------------------------------------------
% srcstruct File containing probStruct
% ---------------------------------------------------------------------
srcstruct         = cfg_files;
srcstruct.tag     = 'srcstruct';
srcstruct.name    = 'File containing probStruct';
srcstruct.help    = {'Select a .mat file containing the probstruct.'};
srcstruct.filter = 'mat';
srcstruct.ufilter = '.*';
srcstruct.num     = [1 Inf];
% ---------------------------------------------------------------------
% srcvar Directly passed variable
% ---------------------------------------------------------------------
srcvar         = cfg_entry;
srcvar.tag     = 'srcvar';
srcvar.name    = 'probStruct';
srcvar.help    = {'Specify a probStruct variable.'};
srcvar.strtype = 'e';
srcvar.num     = [1 Inf];
srcvar.check   = @(job)impexp_run_probstruct2nifti('check','isprobstruct',job);
% ---------------------------------------------------------------------
% srcchoice Input source
% ---------------------------------------------------------------------
srcchoice         = cfg_choice;
srcchoice.tag     = 'srcchoice';
srcchoice.name    = 'Input source';
srcchoice.help    = {'Input probstructs can be either loaded from disk or passed as a variable.'};
srcchoice.values  = {srcstruct srcvar};
% ---------------------------------------------------------------------
% outdir Output directory
% ---------------------------------------------------------------------
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output directory';
outdir.help    = {'Files produced by this function will be written into this output directory.'};
outdir.filter = 'dir';
outdir.ufilter = '.*';
outdir.num     = [1 1];
% ---------------------------------------------------------------------
% fname Output Filename
% ---------------------------------------------------------------------
fname         = cfg_entry;
fname.tag     = 'fname';
fname.name    = 'Output Filename';
fname.help    = {'The output file specified here is written to the selected output directory. '...
    'If there is more than one input probstruct or an probstruct contains '...
    'multiple volumes (series, echos), then the filename part will be extended with running indices.'};
fname.strtype = 's';
fname.num     = [1 Inf];
% ---------------------------------------------------------------------
% outimg Output File & Directory
% ---------------------------------------------------------------------
outimg         = cfg_branch;
outimg.tag     = 'outimg';
outimg.name    = 'Specify Output File & Directory';
outimg.val     = {outdir fname};
outimg.help    = {'Specify a output filename and target directory.'};
% ---------------------------------------------------------------------
% autoimg Automatically generate Output Filename
% ---------------------------------------------------------------------
autoimg        = cfg_const;
autoimg.tag    = 'autoimg';
autoimg.name   = 'Automatically generate Output Filename';
autoimg.val    = {true};
autoimg.help   = {'Generate image filenames automatically. If the sources '...
    'are probstruct files, then their filenames will be used as basename. '...
    'Otherwise, files will be saved to the current directory using running '...
    'indices.'};
% ---------------------------------------------------------------------
% outname Output Naming Scheme
% ---------------------------------------------------------------------
outname         = cfg_choice;
outname.tag     = 'outname';
outname.name    = 'Output Naming Scheme';
outname.values  = {outimg autoimg};
outname.help    = {'Specify a output filename and target directory or use some heuristics.'};
% ---------------------------------------------------------------------
% dtype Data Type
% ---------------------------------------------------------------------
dtype         = cfg_menu;
dtype.tag     = 'dtype';
dtype.name    = 'Data Type';
dtype.help    = {'Data type of output images.'};
dtype.labels = cellstr(spm_type(spm_type));
dtype.values = num2cell(spm_type);
% ---------------------------------------------------------------------
% mapType Map Type
% ---------------------------------------------------------------------
mapType         = cfg_menu;
mapType.tag     = 'mapType';
mapType.name    = 'Map to extract';
mapType.help    = {'Specify the type of map to be extracted. Note that '...
                   'not all probStructs support all map types.'};
mapType.labels = {'ProbMap', ...
                  'Probability Index of a connection to both seed regions (PIBS)',...
                  'Probability Index of forming part of BOI (PIBI)',...
                  'PIBS - PIBI',...
                  'Probability Index of connecting fibre configuration (PICC)',...
                  'Fraction of connecting fibre configuration (FCC)'};
mapType.values = {'probMap',...
                  'sumMap',...
                  'conMap',...
                  'mergMap',...
                  'conWeight',...
                  'conFrac'};
% ---------------------------------------------------------------------
% output Output Options
% ---------------------------------------------------------------------
output         = cfg_branch;
output.tag     = 'output';
output.name    = 'Output Options';
output.val     = {outname dtype mapType};
output.help    = {'Specify a output filename and data type.'};
% ---------------------------------------------------------------------
% probstruct2nifti Convert probStruct to Nifti image(s)
% ---------------------------------------------------------------------
probstruct2nifti         = cfg_exbranch;
probstruct2nifti.tag     = 'probstruct2nifti';
probstruct2nifti.name    = 'Convert probStruct to Nifti image(s)';
probstruct2nifti.val     = {srcchoice output };
probstruct2nifti.help    = {
                        'Convert probstruct image/volume/series to NifTI compatible files'
                        'FORMAT [res, errStr]= probstruct_to_nifti(probstruct, fName)'
                        '======'
                        'This routine tries to convert probstruct data to NifTI compatible files using SPM5 functions spm_type, spm_platform, spm_write_vol. Image data will be saved in user given format (dataType). Default format is int16. If available, a NifTI compatible coordinate transformation matrix will be'
                        'computed from probstruct.edges. If this is not possible, axial orientation'
                        'is assumed and voxel sizes are read from probstruct.vox. If no voxel sizes'
                        'are given, a default of 1x1x1mm will be assumed. Image data are'
                        'rearranged in Nifti/Analyze format and slices are re-ordered depending on'
                        'the handedness of the voxel-to-world coordinate mapping.'
                        ''
                        'Input arguments'
                        'probstruct - probstruct to be converted'
                        'fName    - location and filename for created file(s). If probstruct contains series or multi-echo data, separate volumes will be created for each echo or series image/volume by appending a running series/echo number to the file name.'
                        'dataType - data type for output (''uint8'',''int16'',''int32'',''float32'','
                        '''float64'',''int8'',''uint16'',''uint32'')'
                        'Output arguments'
                        'res      - a cell array containing volume handles (see spm_vol). The array is shaped according to the series/echo dimensions of the probstruct data array.'
                        'errStr   - empty, if successful operation. Otherwise, an error message.'
                        '_______________________________________________________________________'
                        'Bjoern W. Kreher'
                        '08/05'
                        ''
                        'UNIX'
                        '_______________________________________________________________________'
                        '$Id: impexp_cfg_NiftiMrStruct.m,v 1.14 2014/11/03 15:33:05 kellnere Exp $'
}';
probstruct2nifti.prog = @(job)impexp_run_probstruct2nifti('run',job);
probstruct2nifti.vout = @(job)impexp_run_probstruct2nifti('vout',job);
% ---------------------------------------------------------------------
% srcstruct File containing probStruct
% ---------------------------------------------------------------------
srcstruct         = cfg_files;
srcstruct.tag     = 'srcstruct';
srcstruct.name    = 'File containing maskStruct';
srcstruct.help    = {'Select a .mat file containing the maskstruct.'};
srcstruct.filter = 'mat';
srcstruct.ufilter = '.*';
srcstruct.num     = [1 Inf];
% ---------------------------------------------------------------------
% srcvar Directly passed variable
% ---------------------------------------------------------------------
srcvar         = cfg_entry;
srcvar.tag     = 'srcvar';
srcvar.name    = 'maskStruct';
srcvar.help    = {'Specify a maskStruct variable.'};
srcvar.strtype = 'e';
srcvar.num     = [1 Inf];
srcvar.check   = @(job)impexp_run_maskstruct2nifti('check','ismaskstruct',job);
% ---------------------------------------------------------------------
% srcchoice Input source
% ---------------------------------------------------------------------
srcchoice         = cfg_choice;
srcchoice.tag     = 'srcchoice';
srcchoice.name    = 'Input source';
srcchoice.help    = {'Input maskStructs can be either loaded from disk or passed as a variable.'};
srcchoice.values  = {srcstruct srcvar};
% ---------------------------------------------------------------------
% outdir Output directory
% ---------------------------------------------------------------------
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output directory';
outdir.help    = {'Files produced by this function will be written into this output directory.'};
outdir.filter = 'dir';
outdir.ufilter = '.*';
outdir.num     = [1 1];
% ---------------------------------------------------------------------
% fname Output Filename
% ---------------------------------------------------------------------
fname         = cfg_entry;
fname.tag     = 'fname';
fname.name    = 'Output Filename';
fname.help    = {'The output file specified here is written to the selected output directory. '...
    'If there is more than one input maskstruct or an maskstruct contains '...
    'multiple masks, then the filename part will be extended with running indices.'};
fname.strtype = 's';
fname.num     = [1 Inf];
% ---------------------------------------------------------------------
% outimg Output File & Directory
% ---------------------------------------------------------------------
outimg         = cfg_branch;
outimg.tag     = 'outimg';
outimg.name    = 'Specify Output File & Directory';
outimg.val     = {outdir fname};
outimg.help    = {'Specify a output filename and target directory.'};
% ---------------------------------------------------------------------
% autoimg Automatically generate Output Filename
% ---------------------------------------------------------------------
autoimg        = cfg_const;
autoimg.tag    = 'autoimg';
autoimg.name   = 'Automatically generate Output Filename';
autoimg.val    = {true};
autoimg.help   = {'Generate image filenames automatically. If the sources '...
    'are maskstruct files, their filenames and mask names will be used as basename. '...
    'Otherwise, files will be saved to the current directory using running '...
    'indices.'};
% ---------------------------------------------------------------------
% outname Output Naming Scheme
% ---------------------------------------------------------------------
outname         = cfg_choice;
outname.tag     = 'outname';
outname.name    = 'Output Naming Scheme';
outname.values  = {outimg autoimg};
outname.help    = {'Specify a output filename and target directory or use some heuristics.'};
% ---------------------------------------------------------------------
% maskNo Mask Number
% ---------------------------------------------------------------------
maskNo         = cfg_entry;
maskNo.tag     = 'maskNo';
maskNo.name    = 'Mask(s) to extract';
maskNo.strtype = 'n';
maskNo.num     = [1 Inf];
maskNo.help    = {'Specify a vector of mask numbers. A value of ''Inf'' means to extract all masks in a file.'};
% ---------------------------------------------------------------------
% output Output Options
% ---------------------------------------------------------------------
output         = cfg_branch;
output.tag     = 'output';
output.name    = 'Output Options';
output.val     = {outname maskNo};
output.help    = {'Specify a output filename and data type.'};
% ---------------------------------------------------------------------
% maskstruct2nifti Convert maskstruct to Nifti image(s)
% ---------------------------------------------------------------------
maskstruct2nifti         = cfg_exbranch;
maskstruct2nifti.tag     = 'maskstruct2nifti';
maskstruct2nifti.name    = 'Convert maskStruct to Nifti image(s)';
maskstruct2nifti.val     = {srcchoice output };
maskstruct2nifti.help    = {
                        'Convert maskstruct image/volume/series to NifTI compatible files'
                        'FORMAT [res, errStr]= maskstruct_to_nifti(maskstruct, fName)'
                        '======'
                        'This routine tries to convert maskstruct data to NifTI compatible files using SPM5 functions spm_type, spm_platform, spm_write_vol. Image data will be saved in user given format (dataType). Default format is int16. If available, a NifTI compatible coordinate transformation matrix will be'
                        'computed from maskstruct.edges. If this is not possible, axial orientation'
                        'is assumed and voxel sizes are read from maskstruct.vox. If no voxel sizes'
                        'are given, a default of 1x1x1mm will be assumed. Image data are'
                        'rearranged in Nifti/Analyze format and slices are re-ordered depending on'
                        'the handedness of the voxel-to-world coordinate mapping.'
                        ''
                        'Input arguments'
                        'maskstruct - maskstruct to be converted'
                        'fName    - location and filename for created file(s). If maskstruct contains series or multi-echo data, separate volumes will be created for each echo or series image/volume by appending a running series/echo number to the file name.'
                        'dataType - data type for output (''uint8'',''int16'',''int32'',''float32'','
                        '''float64'',''int8'',''uint16'',''uint32'')'
                        'Output arguments'
                        'res      - a cell array containing volume handles (see spm_vol). The array is shaped according to the series/echo dimensions of the maskstruct data array.'
                        'errStr   - empty, if successful operation. Otherwise, an error message.'
                        '_______________________________________________________________________'
                        'Bjoern W. Kreher'
                        '08/05'
                        ''
                        'UNIX'
                        '_______________________________________________________________________'
                        '$Id: impexp_cfg_NiftiMrStruct.m,v 1.14 2014/11/03 15:33:05 kellnere Exp $'
}';
maskstruct2nifti.prog = @(job)impexp_run_maskstruct2nifti('run',job);
maskstruct2nifti.vout = @(job)impexp_run_maskstruct2nifti('vout',job);
% ---------------------------------------------------------------------
% srcimgs Source Images
% ---------------------------------------------------------------------
srcimgs         = cfg_files;
srcimgs.tag     = 'srcimgs';
srcimgs.name    = 'Source Images';
srcimgs.help    = {'Specify here files to be included in the mrStruct. The number of files must match the expected values: 1 for image/volume, any number for series2D/3D and image/volumeEchos, #series-by-#echos for series2D/3DEchos. The decision on whether 2D or 3D data are imported will be made based on the NifTI header of the images. All images need to have same orientation and voxel sizes.'};
srcimgs.filter = 'image';
srcimgs.ufilter = '.*';
srcimgs.num     = [1 Inf];
% ---------------------------------------------------------------------
% srcimgs Source Images (DTD convert
% ---------------------------------------------------------------------
srcimgsTensor         = cfg_files;
srcimgsTensor.tag     = 'srcimgsTensor';
srcimgsTensor.name    = 'Source Images';
srcimgsTensor.help    = {'Specify here tensor files (with prefix Dxx,Dxy,Dxz,Dyy,Dyz,Dzz). Additionaly you may give a seventh image containing a b0 image. All images need to have same orientation and voxel sizes.'};
srcimgsTensor.filter = 'image';
srcimgsTensor.ufilter = '.*';
srcimgsTensor.num     = [1 Inf];
% ---------------------------------------------------------------------
% srcimgs Source Images (DTD convert
% ---------------------------------------------------------------------
srcimgsDW         = cfg_files;
srcimgsDW.tag     = 'srcimgsDW';
srcimgsDW.name    = 'Source Images';
srcimgsDW.help    = {'Specify here DW image files, the image information has to be present in the .mat files. All images need to have same orientation and voxel sizes.'};
srcimgsDW.filter = 'image';
srcimgsDW.ufilter = '.*';
srcimgsDW.num     = [1 Inf];
% ---------------------------------------------------------------------
% srcimgs Source Images (DTD convert
% ---------------------------------------------------------------------
gdirTransformation         = cfg_entry;
gdirTransformation.tag     = 'gdirTransformation';
gdirTransformation.name    = 'Transformation of Gradient Directions';
gdirTransformation.help    = {'Use a transformation different from Id if there are problems with the orientation of the gradient directions'};
gdirTransformation.strtype = 'r';
gdirTransformation.num     = [3 3];
gdirTransformation.extras  = [0 1];
gdirTransformation.def     = @(x) eye(3);
% ---------------------------------------------------------------------
% FSL TOPUP choice
% ---------------------------------------------------------------------
% filename3 name of nifti
% ---------------------------------------------------------------------
% fsltopup_filename         = cfg_files;
% fsltopup_filename.tag     = 'fsltopup_filename';
% fsltopup_filename.name    = 'Select a nifti as reference for reslicing';
% fsltopup_filename.help    = {'Select the nifti (usually something like blah.nii)'};
% fsltopup_filename.filter  = 'image';
% fsltopup_filename.ufilter = '.*';
% fsltopup_filename.num     = [1 1];
fsltopup_filename         = cfg_files;
fsltopup_filename.tag     = 'fsltopup_b0file';
fsltopup_filename.name    = 'Select reversed PE b0 image';
fsltopup_filename.help    = {'Choose b0 image with reversed phase encoding direction'};
fsltopup_filename.filter  = 'image';
fsltopup_filename.ufilter = '.*';
fsltopup_filename.num     = [1 1];

fsltopup_yes            = cfg_branch;
fsltopup_yes.tag     = 'fsltopup_yes';
fsltopup_yes.name    = 'Apply fsltopup';
fsltopup_yes.val     = {fsltopup_filename};
fsltopup_yes.help    = {'Choose if b0 image with reversed PE direction was acquired.'};

fsltopup_no = cfg_const;
fsltopup_no.tag    = 'fsltopup_no';
fsltopup_no.name   = 'No. Choose if PSF Dico was already applied';
fsltopup_no.val    = {true};
fsltopup_no.help   = {'Choose if PSF Dico was already applied'};


fsltopup         = cfg_choice;
fsltopup.tag     = 'fsltopup';
fsltopup.name    = 'FSL topup Distortion Correction';
fsltopup.values  = {fsltopup_yes fsltopup_no};
fsltopup.help    = {'Choose this option if the PSF correction was not performed. A b0 image with inversed PE direction is needed.'};


% ---------------------------------------------------------------------
% thresh Mask threshold (min max)
% ---------------------------------------------------------------------
thresh         = cfg_entry;
thresh.tag     = 'thresh';
thresh.name    = 'Mask threshold (min max)';
thresh.help    = {'Minimum and maximum value to include in ROI (e.g. .8 Inf to include all voxels above 0.8)'};
thresh.strtype = 'e';
thresh.num     = [1 2];
% ---------------------------------------------------------------------
% nfval Mask value for non-finite intensities
% ---------------------------------------------------------------------
nfval         = cfg_menu;
nfval.tag     = 'nfval';
nfval.name    = 'Mask value for non-finite intensities';
nfval.labels  = {'true','false'};
nfval.values  = {true,false};
% ---------------------------------------------------------------------
% outdir Output directory
% ---------------------------------------------------------------------
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output directory';
outdir.help    = {'Files produced by this function will be written into this output directory.'};
outdir.filter = 'dir';
outdir.ufilter = '.*';
outdir.num     = [1 1];
% ---------------------------------------------------------------------
% fname Output Filename
% ---------------------------------------------------------------------
fname         = cfg_entry;
fname.tag     = 'fname';
fname.name    = 'Output Filename';
fname.help    = {'The output file specified here is written to the selected output directory. If the input mrStruct contains multiple volumes (series, echos), then the filename part will be extended with running indices.'};
fname.strtype = 's';
fname.num     = [1 Inf];
% ---------------------------------------------------------------------
% outmat Output File & Directory
% ---------------------------------------------------------------------
outmat         = cfg_branch;
outmat.tag     = 'outmat';
outmat.name    = 'Output File & Directory';
outmat.val     = {outdir fname };
outmat.help    = {'Specify a output filename and target directory.'};
% ---------------------------------------------------------------------
% outvar Output variable
% ---------------------------------------------------------------------
outvar         = cfg_const;
outvar.tag     = 'outvar';
outvar.name    = 'Output variable';
outvar.val{1}  = true;
% ---------------------------------------------------------------------
% outchoice Output destination
% ---------------------------------------------------------------------
outchoice         = cfg_choice;
outchoice.tag     = 'outchoice';
outchoice.name    = 'Output destination';
outchoice.help    = {'Output can be saved to disk or into a MATLAB variable. No check is performed whether the specified file already exists and any previous contents will be overwritten.'};
outchoice.values  = {outmat outvar};
% ---------------------------------------------------------------------
% nifti2maskstruct Convert Nifti images to MaskStruct
% ---------------------------------------------------------------------
nifti2maskstruct         = cfg_exbranch;
nifti2maskstruct.tag     = 'nifti2roistruct';
nifti2maskstruct.name    = 'Convert Nifti images to MaskStruct';
nifti2maskstruct.val     = {srcimgs thresh nfval outchoice };
nifti2maskstruct.help    = {
    'Convert NIFTI image(s) to mask structs.'
                   }';
nifti2maskstruct.prog = @(job)impexp_run_nifti2maskstruct('run',job);
nifti2maskstruct.vout = @(job)impexp_run_nifti2maskstruct('vout',job);
% ---------------------------------------------------------------------
% nifti2DTDstruct Convert Nifti images to DTDStruct
% ---------------------------------------------------------------------
nifti2DTDstruct         = cfg_exbranch;
nifti2DTDstruct.tag     = 'nifti2DTDstruct';
nifti2DTDstruct.name    = 'Convert Nifti tensor images to DTDStruct';
nifti2DTDstruct.val     = {srcimgsTensor gdirTransformation outchoice };
nifti2DTDstruct.help    = {'Convert NIFTI tensor image(s) to DTDstruct.'
                   }';
nifti2DTDstruct.prog = @(job)impexp_run_nifti2DTDstruct('run',job);
nifti2DTDstruct.vout = @(job)impexp_run_nifti2DTDstruct('vout',job);
% ---------------------------------------------------------------------
% nifti2HARDIstruct Convert Nifti images to MRStruct (HARDI)
% ---------------------------------------------------------------------
nifti2HARDIstruct         = cfg_exbranch;
nifti2HARDIstruct.tag     = 'nifti2HARDIstruct';
nifti2HARDIstruct.name    = 'Convert Nifti DW images to mrStruct (HARDI)';
nifti2HARDIstruct.val     = {srcimgsDW fsltopup gdirTransformation outchoice };
nifti2HARDIstruct.help    = {'Convert NIFTI DW images to mrStruct (HARDI).'
                   }';
nifti2HARDIstruct.prog = @(job)impexp_run_nifti2HARDIstruct('run',job);
nifti2HARDIstruct.vout = @(job)impexp_run_nifti2HARDIstruct('vout',job);
% ---------------------------------------------------------------------
% srcstruct File containing probStruct
% ---------------------------------------------------------------------
srcstruct         = cfg_files;
srcstruct.tag     = 'srcstruct';
srcstruct.name    = 'File containing maskStruct';
srcstruct.help    = {'Select a .mat file containing the maskstruct.'};
srcstruct.filter = 'mat';
srcstruct.ufilter = '.*';
srcstruct.num     = [1 1];
% ---------------------------------------------------------------------
% srcvar Directly passed variable
% ---------------------------------------------------------------------
srcvar         = cfg_entry;
srcvar.tag     = 'srcvar';
srcvar.name    = 'maskStruct';
srcvar.help    = {'Specify a maskStruct variable.'};
srcvar.strtype = 'e';
srcvar.num     = [1 1];
srcvar.check   = @(job)impexp_run_maskstruct2nifti('check','ismaskstruct',job);
% ---------------------------------------------------------------------
% srcchoice Input source
% ---------------------------------------------------------------------
srcchoice         = cfg_choice;
srcchoice.tag     = 'srcchoice';
srcchoice.name    = 'Input source';
srcchoice.help    = {'Input maskStructs can be either loaded from disk or passed as a variable.'};
srcchoice.values  = {srcstruct srcvar};
% ---------------------------------------------------------------------
% maskname Mask Name
% ---------------------------------------------------------------------
maskname         = cfg_entry;
maskname.tag     = 'masknames';
maskname.strtype = 's';
maskname.num     = [1 Inf];
maskname.name   = 'Mask Name';
% ---------------------------------------------------------------------
% masknames Mask Names
% ---------------------------------------------------------------------
masknames        = cfg_repeat;
masknames.tag    = 'masknames';
masknames.num    = [1 Inf];
masknames.name   = 'Mask Names';
masknames.values = {maskname};
masknames.help   = {'Enter one or more mask names.'};
% ---------------------------------------------------------------------
% msnames2ind Convert Mask Names to MaskStruct Indices
% ---------------------------------------------------------------------
msnames2ind         = cfg_exbranch;
msnames2ind.tag     = 'msnames2ind';
msnames2ind.name    = 'Convert Mask Names to MaskStruct Indices';
msnames2ind.val     = {srcchoice masknames };
msnames2ind.help    = {
    'Query a MaskStruct and return a list of indices for a list of names.'
                   }';
msnames2ind.prog = @(job)impexp_run_msnames2ind('run',job);
msnames2ind.vout = @(job)impexp_run_msnames2ind('vout',job);
% ---------------------------------------------------------------------
% impexp_NiftiMrStruct NiftiMrStruct
% ---------------------------------------------------------------------
impexp_NiftiMrStruct         = cfg_choice;
impexp_NiftiMrStruct.tag     = 'impexp_NiftiMrStruct';
impexp_NiftiMrStruct.name    = 'NiftiMrStruct';
impexp_NiftiMrStruct.help    = {
                                '_______________________________________________________________________'
                                ''
                                'This toolbox contains various helper functions to convert NifTI images to and from MedPhys mrStructs.'
                                ''
                                'This toolbox is free but copyright software, distributed under the terms of the GNU General Public Licence as published by the Free Software Foundation (either version 2, as given in file spm_LICENCE.man, or at your option, any later version). Further details on "copyleft" can be found at http://www.gnu.org/copyleft/.'
                                'The toolbox consists of the files listed in its Contents.m file.'
                                '_______________________________________________________________________'
                                ''
                                '@(#) $Id: impexp_cfg_NiftiMrStruct.m,v 1.14 2014/11/03 15:33:05 kellnere Exp $'
}';
impexp_NiftiMrStruct.values  = {nifti2mrstruct mrstruct2nifti bo2nifti fa2nifti trace2nifti eigVal1_2_nifti eigVal2_2_nifti eigVal3_2_nifti nifti2maskstruct ...
                   maskstruct2nifti nifti2DTDstruct nifti2HARDIstruct msnames2ind probstruct2nifti impexp_cfg_ftrstruct_Curves2dx};

