function ea_switchspace(~,~,spacename,mute)
if strcmp(spacename(1:2),'->')
    return
end

if ~exist('mute','var')
    answ=questdlg('Please be aware that switching the default template space is a critical and more or less complicated step. Not all functions may perfectly work if you switch to a different default template as opposed to the ICBM2009b Nonlinear Asymmetric series. Are you sure you wish to switch to a different space?','Switch anatomical space','Sure','Cancel','Cancel');
else
    answ='Sure';
end
if strcmp(answ,'Sure')
    ea_storemachineprefs('space',spacename);

    % Get space definitions
    spacedef = ea_getspacedef;

    % Check if there's space specific togglestates setting
    if ~isfield(spacedef, 'togglestates') % Fallback
        mnispace = [ea_getearoot,'templates',filesep,'space',filesep,'MNI_ICBM_2009b_NLIN_ASYM'];
        load([mnispace, filesep, 'ea_space_def.mat'], 'spacedef');
    end
    ea_storemachineprefs('togglestates', spacedef.togglestates);

    % Check if there's space specific view setting
    if ~isfield(spacedef, 'view') % Fallback
        mnispace = [ea_getearoot,'templates',filesep,'space',filesep,'MNI_ICBM_2009b_NLIN_ASYM'];
        load([mnispace, filesep, 'ea_space_def.mat'], 'spacedef');
    end
    ea_storemachineprefs('view', spacedef.view);

    % Check if there's space specific default atlas setting
    if ~isfield(spacedef, 'defaultatlas') % Fallback
        mnispace = [ea_getearoot,'templates',filesep,'space',filesep,'MNI_ICBM_2009b_NLIN_ASYM'];
        load([mnispace, filesep, 'ea_space_def.mat'], 'spacedef');
    end
    ea_storemachineprefs('defaultatlas', spacedef.defaultatlas);

    if ~exist('mute','var')
        disp('Restarting Lead Neuroimaging Suite...');
        close all force
        lead;
    end
end
