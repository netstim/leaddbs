function ea_patch_spm
% Patch SPM def, norm and shoot cfg file to make deformation file whose name does
% not start with 'y_*' or 'iy_*' also recognizable.

unzip(fullfile(ea_getearoot, 'ext_libs', 'spm', 'cfg_patch.zip'), fileparts(which('spm')));

% defCfgFile = [fileparts(which('spm')), filesep, 'config', filesep, 'spm_cfg_deformations.m'];
% defCfg = fileread(defCfgFile);
% if contains(defCfg, '.*y_.*\.nii$')
%     defCfg = strrep(defCfg, '.*y_.*\.nii$', '.*\.nii$');
%     fid = fopen(defCfgFile, 'w');
%     fprintf(fid, '%s', defCfg);
%     fclose(fid);
% end
% 
% normCfgFile = [fileparts(which('spm')), filesep, 'config', filesep, 'spm_cfg_norm.m'];
% normCfg = fileread(normCfgFile);
% if contains(normCfg, 'y_.*\.nii$')
%     normCfg = strrep(normCfg, 'y_.*\.nii$', '.*\.nii$');
%     fid = fopen(normCfgFile, 'w');
%     fprintf(fid, '%s', normCfg);
%     fclose(fid);
% end
% 
% shootCfgFile = [fileparts(which('spm')), filesep, 'toolbox', filesep, 'Shoot', filesep, 'tbx_cfg_shoot.m'];
% shootCfg = fileread(shootCfgFile);
% if contains(shootCfg, '^y_.*')
%     shootCfg = strrep(shootCfg, '^y_.*', '.*');
%     fid = fopen(shootCfgFile, 'w');
%     fprintf(fid, '%s', shootCfg);
%     fclose(fid);
% end
% 
% spm_jobman('initcfg');
