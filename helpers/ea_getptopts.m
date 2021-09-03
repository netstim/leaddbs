function options = ea_getptopts(directory,options)
% Generate minimal options struct from a patient directory.

if isempty(directory)
    directory = pwd;
end

directory = regexp(directory, ['.*derivatives\' filesep 'leaddbs\' filesep 'sub-[^' filesep ']*'], 'match', 'once'); % remove after subid if present

options.earoot = ea_getearoot;
options.native = 0;

if contains(directory, ['derivatives', filesep, 'leaddbs'])
    % Construct BIDS class, get subjId
    BIDSRoot = regexp(directory, ['^.*(?=\', filesep, 'derivatives)'], 'match', 'once');
    bids = BIDSFetcher(BIDSRoot);
    subjId = regexp(directory, ['(?<=leaddbs\', filesep, 'sub-).*'], 'match', 'once');

    % Set prefs (not really needed, leave it for now)
    options.prefs = ea_prefs(['sub-', subjId]);

    % Set modality field
    uiPrefsFile = bids.getPrefs(subjId, 'uiprefs', 'mat');
    if isfile(uiPrefsFile)
        uiprefs = load(uiPrefsFile);
        options.modality = uiprefs.modality;
        options.sides = uiprefs.sides;
    else
        options.modality = bids.settings.preferMRCT;
        options.sides = [1,2];
    end

    % Set subj BIDS struct
    options.subj = bids.getSubj(subjId, options.modality);

    % Set elmodel and elspec
    recon = bids.getRecon(subjId);
    if isfile(recon.recon)
        load(recon.recon);
        try
            options.elmodel = ea_get_first_notempty_elmodel(reco.props);
        catch
            options.elmodel = 'Medtronic 3389';
        end
        options = ea_resolve_elspec(options);
    end
else
    error('Not a BIDS dataset!')
end
