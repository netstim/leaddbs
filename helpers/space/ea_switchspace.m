function ea_switchspace(~,~,spacename,mute)
if startsWith(spacename,'->')
    return
end

if ~exist('mute','var')
    answ=questdlg('Please be aware that switching the default template space is a critical and more or less complicated step. Not all functions may perfectly work if you switch to a different default template as opposed to the ICBM2009b Nonlinear Asymmetric series. Are you sure you wish to switch to a different space?','Switch anatomical space','Sure','Cancel','Cancel');
else
    answ='Sure';
end
if strcmp(answ,'Sure')
    ea_setprefs('space',spacename);

    % Get space definitions
    spacedef = ea_getspacedef;

    % Check if there's space specific togglestates setting
    if ~isfield(spacedef, 'togglestates') % Fallback
        mnispace = [ea_getearoot,'templates',filesep,'space',filesep,'MNI152NLin2009bAsym'];
        load([mnispace, filesep, 'spacedef.mat'], 'spacedef');
    end
    ea_setprefs('togglestates', spacedef.togglestates);

    % Check if there's space specific view setting
    if ~isfield(spacedef, 'view') % Fallback
        mnispace = [ea_getearoot,'templates',filesep,'space',filesep,'MNI152NLin2009bAsym'];
        load([mnispace, filesep, 'spacedef.mat'], 'spacedef');
    end
    ea_setprefs('view', spacedef.view);

    % Check if there's space specific default atlas setting
    if ~isfield(spacedef, 'defaultatlas') % Fallback
        mnispace = [ea_getearoot,'templates',filesep,'space',filesep,'MNI152NLin2009bAsym'];
        load([mnispace, filesep, 'spacedef.mat'], 'spacedef');
    end
    ea_setprefs('defaultatlas', spacedef.defaultatlas);

    % Set tensorFileName here
    if strcmp(spacename, 'MNI152NLin2009bAsym')
        ea_setprefs('vatsettings.butenko_tensorFileName', 'IITMeanTensor.nii.gz');
    elseif strcmp(spacename, 'WaxholmSpaceSDRat')
        ea_setprefs('vatsettings.butenko_tensorFileName', 'JohnsonWS.nii.gz');
    else
        ea_setprefs('vatsettings.butenko_tensorFileName', '')
    end

    if ~exist('mute','var')
        disp('Restarting Lead Neuroimaging Suite...');
        close all force
        lead;
    end
end
