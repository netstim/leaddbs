function record = ea_lc_coreg_record(options)
% Read the current B0 registration contract. Never infer it from file dates.
directory = fullfile(options.root, options.patientname);
protocol = fullfile(directory, 'coregistration', 'log', 'ea_lc_b0_coreg.mat');
if ~isfile(protocol)
    error('LeadDBS:MissingB0CoregProtocol', ...
        'No current B0 coregistration protocol. Rerun B0 coregistration in Check Coregistration once.');
end
data = load(protocol, 'record');
record = data.record;
fields = {'moving', 'fixed', 'anchor', 'output', 'forward', 'inverse'};
for i = 1:numel(fields)
    path = fullfile(directory, record.(fields{i}));
    if ~isfile(path)
        error('LeadDBS:MissingB0CoregResult', ...
            'The B0 coregistration protocol requires %s. Rerun B0 coregistration.', path);
    end
    record.(fields{i}) = path;
end
end
