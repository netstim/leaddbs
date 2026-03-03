function ea_patch_spm
% Patch SPM def, norm and shoot cfg file to make deformation file whose name does
% not start with 'y_*' or 'iy_*' also recognizable.

% unzip(fullfile(ea_getearoot, 'ext_libs', 'spm', 'cfg_patch.zip'), fileparts(which('spm')));

patchedFile = {};

spmDir = fileparts(which('spm'));
defCfgFile = fullfile(spmDir, 'config', 'spm_cfg_deformations.m');
normCfgFile = fullfile(spmDir, 'config', 'spm_cfg_norm.m');
shootCfgFile = fullfile(spmDir, 'toolbox', 'Shoot', 'tbx_cfg_shoot.m');
platformCfgFile = fullfile(spmDir, 'spm_platform.m');

defCfg = fileread(defCfgFile);
if contains(defCfg, '.*y_.*\.nii$')
    patchedFile = [patchedFile, defCfgFile];
    defCfg = strrep(defCfg, '.*y_.*\.nii$', '.*\.nii$');
    fid = fopen(defCfgFile, 'w');
    fprintf(fid, '%s', defCfg);
    fclose(fid);
elseif contains(defCfg, "def          = files('Deformation Field','def','.*\.nii$',[1 1]);")
    patchedFile = [patchedFile, defCfgFile];
end

normCfg = fileread(normCfgFile);
if contains(normCfg, 'y_.*\.nii$')
    patchedFile = [patchedFile, normCfgFile];
    normCfg = strrep(normCfg, 'y_.*\.nii$', '.*\.nii$');
    fid = fopen(normCfgFile, 'w');
    fprintf(fid, '%s', normCfg);
    fclose(fid);
elseif contains(normCfg, "def.ufilter = '.*\.nii$';")
    patchedFile = [patchedFile, normCfgFile];
end

shootCfg = fileread(shootCfgFile);
if contains(shootCfg, '^y_.*')
    shootCfg = strrep(shootCfg, '^y_.*', '.*');
    fid = fopen(shootCfgFile, 'w');
    fprintf(fid, '%s', shootCfg);
    fclose(fid);
    patchedFile = [patchedFile, shootCfgFile];
elseif contains(shootCfg, "deformation.ufilter = '.*\.nii$';")
    patchedFile = [patchedFile, shootCfgFile];
end

spm_jobman('initcfg');

if startsWith(spm('version'), 'SPM12')
    platformCfg = fileread(platformCfgFile);
    platformCfgPatched = 0;
    if ~isempty(regexp(platformCfg, "'MACI64',    'unx',   0;...\s+'GLNX86'", 'once'))
        platformCfg = regexprep(platformCfg, "('MACI64',    'unx',   0;...)(\s+)('GLNX86')", "$1$2'MACA64',    'unx',   0;...$2$3");
        platformCfgPatched = 1;
    elseif contains(platformCfg, "'MACA64',    'unx',   0;...")
        platformCfgPatched = 1;
    end
    if contains(platformCfg, "'MACI64'}")
        platformCfg = strrep(platformCfg, "'MACI64'}", "'MACI64','MACA64'}");
        platformCfgPatched = 1;
    elseif contains(platformCfg, "'MACI64','MACA64'}")
        platformCfgPatched = 1;
    end
    if platformCfgPatched
        patchedFile = [patchedFile, platformCfgFile];
    end
    fid = fopen(platformCfgFile, 'w');
    fprintf(fid, '%s', platformCfg);
    fclose(fid);
end

if ~isempty(patchedFile)
    savejson('pathedFile', patchedFile, fullfile(ea_prefsdir, 'SPMPatched.json'));
end
