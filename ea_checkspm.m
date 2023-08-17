function ea_checkspm

if ~isdeployed
    try
        ver = spm('version');
    catch
        ea_error('SPM12 seems not installed. Please install SPM12 and add it to the MATLAB path before using Lead-DBS.');
    end

    try
        spm_check_installation('basic'); % have SPM check its path and binaries and give proper warnings and hints to the user if something is not good
    catch ME
        if strcmp(computer('arch'), 'maca64') && contains(ME.message, 'MEX files')
            % Add maca64 MEX files from upstream
            ea_cprintf('CmdWinWarnings', 'Adding missing maca64 MEX files for SPM12...\n');
            unzip(fullfile(ea_getearoot, 'ext_libs', 'spm', 'spm_mexmaca64.zip'), fileparts(which('spm')));
        else
            rethrow(ME);
        end
    end

    % Patch SPM cfg files
    cfg = dir(fullfile(fileparts(which('spm')), 'toolbox', 'Shoot', 'tbx_cfg_shoot.m'));
    if cfg.bytes ~= 35562
        ea_patch_spm;
        ea_cprintf('CmdWinWarnings', 'Patched SPM cfg files for use in LeadDBS.\n')
    end

    if ~any(ismember(ver(8:11),'.')) % old version format
        if str2double(ver(8:11))<6906
            msgbox('Some functions (such as SPM SHOOT and DARTEL) may not be available using your SPM version. Please upgrade SPM12 to at least revision 6906 (or simply update to newest release).');
        end
    end
end