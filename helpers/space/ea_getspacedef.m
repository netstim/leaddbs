function spacedef = ea_getspacedef

try
    load([ea_space, 'spacedef.mat'], 'spacedef');
catch ME
    if ~isempty(which('lead_dbs.mlapp')) && contains(ea_getspace, '_')
        % LeadDBS v3.0+ installed but prefs from classic version found.
        ea_restoreprefs(1, 1);
        ea_cprintf('CmdWinWarnings', 'Incompatible prefs found! Restored to default.\n');
        load([ea_space, 'spacedef.mat'], 'spacedef');
    else
        rethrow(ME);
    end
end
