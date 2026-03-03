function ea_checkspm

if ~isdeployed
    try
        ver = spm('version');
    catch
        ea_error('SPM seems not installed. Please install SPM and add it to the MATLAB path before using Lead-DBS.');
    end

    if startsWith(ver, 'SPM12')
        verNumber = regexp(ver, '(?<=SPM12 \()(\d+)(?=\))', 'match', 'once');
        if ~isempty(verNumber) % Will be empty when using dev version of SPM
            if str2double(verNumber)<6906 % Old version of SPM detected
                msgbox('Some functions (such as SPM SHOOT and DARTEL) may not be available using your SPM version. Please upgrade SPM12 to at least revision 6906 (or simply update to newest release).');
            end
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
    end

    % Patch SPM cfg files
    SPMPatched = 0;
    spmDir = fileparts(which('spm'));
    targets = {
        fullfile(spmDir, 'config', 'spm_cfg_deformations.m')
        fullfile(spmDir, 'config', 'spm_cfg_norm.m')
        fullfile(spmDir, 'toolbox', 'Shoot', 'tbx_cfg_shoot.m')
    };
    if startsWith(ver, 'SPM12')
        targets{end+1} = fullfile(spmDir, 'spm_platform.m');
    end
    flagFile = fullfile(ea_prefsdir, 'SPMPatched.json');
    if isfile(flagFile)
        flag = loadjson(flagFile);
        if all(contains(flag.pathedFile, targets))
            SPMPatched = 1;
        end
    end

    if ~SPMPatched
        ea_patch_spm;
        ea_cprintf('CmdWinWarnings', 'Patched SPM cfg files for use in LeadDBS.\n')
    end
end