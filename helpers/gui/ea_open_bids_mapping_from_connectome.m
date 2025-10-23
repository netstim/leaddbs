function ea_open_bids_mapping_from_connectome(BIDSRoot, subjId)
% Open main Lead-DBS GUI and switch to the BIDS Import (mapping) tab.
% Tries to preselect dataset and subject when possible.

if ~exist('BIDSRoot','var') || isempty(BIDSRoot) || ~isfolder(BIDSRoot)
    lead; % fallback
    return;
end

try
    % Open main GUI and get figure handle
    hLead = lead; % returns figure handle in OutputFcn
    drawnow;

    % Try to get handles struct
    h = guihandles(hLead);

    % Set dataset root if available
    if isfield(h, 'datasetselect') && isgraphics(h.datasetselect)
        try
            h.datasetselect.String = BIDSRoot;
            h.datasetselect.TooltipString = BIDSRoot;
        end
    end

    % Switch to Import tab if present
    if isfield(h, 'processtabgroup') && isfield(h, 'importtab')
        try
            h.processtabgroup.SelectedTab = h.importtab;
        end
    end

    % Preselect subject in patient list if provided
    if exist('subjId','var') && ~isempty(subjId) && isfield(h,'patientlist') && isgraphics(h.patientlist)
        try
            if ischar(subjId), subjId = {subjId}; end
            % Build/refresh patient list using existing utilities if available
            try
                options.prefs = ea_prefs;
                ea_getpatients(options,h);
            catch
            end
            % Set selection
            if isprop(h.patientlist,'Data')
                ids = h.patientlist.Data.subjId;
                sel = find(ismember(ids, subjId{1}));
                if ~isempty(sel)
                    h.patientlist.Selection = sel(1);
                end
            end
        end
    end
catch
    % Silently ignore - mapper GUI is still opened
end


