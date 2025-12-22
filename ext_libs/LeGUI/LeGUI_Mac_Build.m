function LeGUI_Mac_Build(varargin)
% -a adds files to compiled executable, -N clears matlab path except for main matlab folder, -p adds toolbox
% cfgroot is the root location of archived files (-a flag) for deployed code (mfilename('fullpath') works in a deployed environment) 

Version = '1.2';
if nargin
    Version = varargin{1};
end
disp(['Building version ',Version]);
RootDir = fileparts(mfilename('fullpath'));
BuildDir = fullfile(RootDir,'build');
MainFile = fullfile(RootDir,'LeGUI.mlapp');
DepCell = {
    'tpm'
    'nmm'
    'icons'
    'ReadMe.txt'
    'Contents.txt'
    '@file_array/private/file2mat.mexa64'
    '@file_array/private/file2mat.mexmaci64'
    '@file_array/private/file2mat.mexw32'
    '@file_array/private/file2mat.mexw64'
    '@file_array/private/init.mexa64'
    '@file_array/private/init.mexmaci64'
    '@file_array/private/init.mexw32'
    '@file_array/private/init.mexw64'
    '@file_array/private/mat2file.mexa64'
    '@file_array/private/mat2file.mexmaci64'
    '@file_array/private/mat2file.mexw32'
    '@file_array/private/mat2file.mexw64'
    'toolbox/FieldMap/pm_create_connectogram_dtj.mexa64'
    'toolbox/FieldMap/pm_create_connectogram_dtj.mexmaci64'
    'toolbox/FieldMap/pm_create_connectogram_dtj.mexw32'
    'toolbox/FieldMap/pm_create_connectogram_dtj.mexw64'
    'toolbox/FieldMap/pm_estimate_ramp.mexa64'
    'toolbox/FieldMap/pm_estimate_ramp.mexmaci64'
    'toolbox/FieldMap/pm_estimate_ramp.mexw32'
    'toolbox/FieldMap/pm_estimate_ramp.mexw64'
    'toolbox/FieldMap/pm_ff_unwrap.mexa64'
    'toolbox/FieldMap/pm_ff_unwrap.mexmaci64'
    'toolbox/FieldMap/pm_ff_unwrap.mexw32'
    'toolbox/FieldMap/pm_ff_unwrap.mexw64'
    'toolbox/FieldMap/pm_invert_phasemap_dtj.mexa64'
    'toolbox/FieldMap/pm_invert_phasemap_dtj.mexmaci64'
    'toolbox/FieldMap/pm_invert_phasemap_dtj.mexw32'
    'toolbox/FieldMap/pm_invert_phasemap_dtj.mexw64'
    'toolbox/FieldMap/pm_merge_regions.mexa64'
    'toolbox/FieldMap/pm_merge_regions.mexmaci64'
    'toolbox/FieldMap/pm_merge_regions.mexw32'
    'toolbox/FieldMap/pm_merge_regions.mexw64'
    'toolbox/FieldMap/pm_pad.mexa64'
    'toolbox/FieldMap/pm_pad.mexmaci64'
    'toolbox/FieldMap/pm_pad.mexw32'
    'toolbox/FieldMap/pm_pad.mexw64'
    'toolbox/FieldMap/pm_restore_ramp.mexa64'
    'toolbox/FieldMap/pm_restore_ramp.mexmaci64'
    'toolbox/FieldMap/pm_restore_ramp.mexw32'
    'toolbox/FieldMap/pm_restore_ramp.mexw64'
    'toolbox/FieldMap/pm_smooth_phasemap_dtj.mexa64'
    'toolbox/FieldMap/pm_smooth_phasemap_dtj.mexmaci64'
    'toolbox/FieldMap/pm_smooth_phasemap_dtj.mexw32'
    'toolbox/FieldMap/pm_smooth_phasemap_dtj.mexw64'
    'spm_add.mexa64'
    'spm_add.mexmaci64'
    'spm_add.mexw32'
    'spm_add.mexw64'
    'toolbox/OldNorm/spm_brainwarp.mexa64'
    'toolbox/OldNorm/spm_brainwarp.mexmaci64'
    'toolbox/OldNorm/spm_brainwarp.mexw32'
    'toolbox/OldNorm/spm_brainwarp.mexw64'
    'spm_bsplinc.mexa64'
    'spm_bsplinc.mexmaci64'
    'spm_bsplinc.mexw32'
    'spm_bsplinc.mexw64'
    'spm_bsplins.mexa64'
    'spm_bsplins.mexmaci64'
    'spm_bsplins.mexw32'
    'spm_bsplins.mexw64'
    'spm_bwlabel.mexa64'
    'spm_bwlabel.mexmaci64'
    'spm_bwlabel.mexw32'
    'spm_bwlabel.mexw64'
    'spm_cat.mexa64'
    'spm_cat.mexmaci64'
    'spm_cat.mexw32'
    'spm_cat.mexw64'
    'spm_conv_vol.mexa64'
    'spm_conv_vol.mexmaci64'
    'spm_conv_vol.mexw32'
    'spm_conv_vol.mexw64'
    'spm_dicom_dict.mat'
    'spm_dicom_dict.txt'
    'spm_diffeo.mexa64'
    'spm_diffeo.mexmaci64'
    'spm_diffeo.mexw32'
    'spm_diffeo.mexw64'
    'spm_dilate_erode.mexa64'
    'spm_dilate_erode.mexmaci64'
    'spm_dilate_erode.mexw32'
    'spm_dilate_erode.mexw64'
    'spm_existfile.mexa64'
    'spm_existfile.mexmaci64'
    'spm_existfile.mexw32'
    'spm_existfile.mexw64'
    'spm_field.mexa64'
    'spm_field.mexmaci64'
    'spm_field.mexw32'
    'spm_field.mexw64'
    'spm_gamrnd.mexa64'
    'spm_gamrnd.mexmaci64'
    'spm_gamrnd.mexw32'
    'spm_gamrnd.mexw64'
    'spm_get_lm.mexa64'
    'spm_get_lm.mexmaci64'
    'spm_get_lm.mexw32'
    'spm_get_lm.mexw64'
    'spm_global.mexa64'
    'spm_global.mexmaci64'
    'spm_global.mexw32'
    'spm_global.mexw64'
    'spm_hist.mexa64'
    'spm_hist.mexmaci64'
    'spm_hist.mexw32'
    'spm_hist.mexw64'
    'spm_hist2.mexa64'
    'spm_hist2.mexmaci64'
    'spm_hist2.mexw32'
    'spm_hist2.mexw64'
    'spm_krutil.mexa64'
    'spm_krutil.mexmaci64'
    'spm_krutil.mexw32'
    'spm_krutil.mexw64'
    'spm_mesh_utils.mexa64'
    'spm_mesh_utils.mexmaci64'
    'spm_mesh_utils.mexw32'
    'spm_mesh_utils.mexw64'
    'spm_mrf.mexa64'
    'spm_mrf.mexmaci64'
    'spm_mrf.mexw32'
    'spm_mrf.mexw64'
    'spm_project.mexa64'
    'spm_project.mexmaci64'
    'spm_project.mexw32'
    'spm_project.mexw64'
    'spm_render_vol.mexa64'
    'spm_render_vol.mexmaci64'
    'spm_render_vol.mexw32'
    'spm_render_vol.mexw64'
    'spm_resels_vol.mexa64'
    'spm_resels_vol.mexmaci64'
    'spm_resels_vol.mexw32'
    'spm_resels_vol.mexw64'
    'spm_sample_vol.mexa64'
    'spm_sample_vol.mexmaci64'
    'spm_sample_vol.mexw32'
    'spm_sample_vol.mexw64'
    'spm_slice_vol.mexa64'
    'spm_slice_vol.mexmaci64'
    'spm_slice_vol.mexw32'
    'spm_slice_vol.mexw64'
    'spm_unlink.mexa64'
    'spm_unlink.mexmaci64'
    'spm_unlink.mexw32'
    'spm_unlink.mexw64'
    '@xmltree/private/xml_findstr.mexa64'
    '@xmltree/private/xml_findstr.mexmaci64'
    '@xmltree/private/xml_findstr.mexw32'
    '@xmltree/private/xml_findstr.mexw64'
    '@gifti/private/zstream.mexa64'
    '@gifti/private/zstream.mexmaci64'
    '@gifti/private/zstream.mexw64'
    };

DepCellFull = cellfun(@(x,y)fullfile([' -a ',x],y),repmat({RootDir},length(DepCell),1),DepCell,'uniformoutput',false);
DepStrFull = cell2mat(DepCellFull');

MccStr = ['mcc -v -m ',MainFile,' -d ',BuildDir,' -o LeGUI_Mac',DepStrFull];
eval(MccStr);

zipname = fullfile(BuildDir,['LeGUI_Mac_v',Version,'.zip']);
exename = fullfile(BuildDir,'LeGUI_Mac.app');
atlasname = fullfile(BuildDir,'atlases');
excludedname = fullfile(BuildDir,'mccExcludedFiles.log');
readmename = fullfile(BuildDir,'readme.txt');
requiredname = fullfile(BuildDir,'requiredMCRProducts.txt');
scriptname = fullfile(BuildDir,'run_LeGUI_Mac.sh');
includedname = fullfile(BuildDir,'includedSupportPackages.txt');
unresolvedname = fullfile(BuildDir,'unresolvedSymbols.txt');

zip(zipname,{exename,atlasname,scriptname});

try
    rmdir(exename,'s')
catch
    rmdir(exename,'s')
end
delete(excludedname, readmename, requiredname, scriptname, includedname, unresolvedname);

% opts = compiler.package.InstallerOptions(...
%     'ApplicationName',['LeGUI_Mac_v',Version],...
%     'AuthorCompany','University of Utah',...
%     'AuthorName','Tyler Davis',...
%     'InstallerName',['LeGUI_Mac_Installer_v',Version],...
%     'OutputDir',BuildDir,...
%     'Version',Version,...
%     'Summary','Software for localizing intracranial electrodes');
% 
% compiler.package.installer(fullfile(BuildDir,'LeGUI_Mac.app'),fullfile(BuildDir,'requiredMCRProducts.txt'),'Options',opts)


