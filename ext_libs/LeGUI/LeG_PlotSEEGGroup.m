function LeG_PlotSEEGGroup(rootDir, patients, renderMode)
% Plot multiple SEEG reconstructions in one standard Lead-DBS MNI figure.
%
% rootDir should point to the BIDS dataset root.
%
% USAGE:
%   LeG_PlotSEEGGroup('/path/to/dataset', [], 'sphere')
%   LeG_PlotSEEGGroup('/path/to/dataset', [], 'electrode')
%
%   patients = ["sub-72","sub-84","sub-90"];
%   LeG_PlotSEEGGroup('/path/to/dataset', patients, 'sphere')


if nargin < 3 || isempty(renderMode)
    renderMode = 'sphere';
end

if nargin < 2
    patients = [];
end

if ~ismember(lower(renderMode), {'sphere','electrode'})
    error('renderMode must be ''sphere'' or ''electrode''.');
end


% Choose dataset folder if needed

if nargin < 1 || isempty(rootDir)

    rootDir = uigetdir(pwd, ...
        'Select BIDS dataset folder');

    if isequal(rootDir,0)
        return
    end
end


% Lead-DBS derivatives folder

leaddbsDir = fullfile(rootDir, 'derivatives', 'leaddbs');

if ~isfolder(leaddbsDir)
    error('Could not find derivatives/leaddbs inside:\n%s', rootDir);
end


% Find reconstruction files

if isempty(patients)

    recoFiles = dir(fullfile( ...
        leaddbsDir, '**', '*_desc-reconstruction.mat'));

else

    patients = string(patients);
    recoFiles = [];

    for p = 1:numel(patients)

        patientName = char(patients(p));

        patientDir = fullfile(leaddbsDir, patientName);

        r = dir(fullfile( ...
            patientDir, ...
            'reconstruction', ...
            '*_desc-reconstruction.mat'));

        if isempty(r)
            warning('No reconstruction found for %s.', patientName);
            continue
        end

        recoFiles = [recoFiles; r(:)];

    end
end

if isempty(recoFiles)
    error('No reconstruction files found.');
end

fprintf('\nFound %d reconstruction(s).\n', numel(recoFiles));


% Lead-DBS MNI figure

resultfig = ea_mnifigure;

figure(resultfig);
hold on;


% Rendering options

options = struct;

options.prefs = ea_prefs('');

options.d3.elrendering = 1;
options.d3.hlactivecontacts = 0;

switch lower(renderMode)

    case 'sphere'
        options.prefs.machine.d2.seegRenderMode = 'Spheres';

    case 'electrode'
        options.prefs.machine.d2.seegRenderMode = 'Electrode model';

end


% One unique color per subject/reconstruction

nSubjects = numel(recoFiles);

subjectHues = linspace(0, 1, nSubjects + 1)';
subjectHues(end) = [];

% Offset so the first subject is not always red
subjectHues = mod(subjectHues + 0.58, 1);

subjectColors = hsv2rgb([ ...
    subjectHues, ...
    0.70*ones(nSubjects,1), ...
    0.85*ones(nSubjects,1)]);


% Loop through every reconstruction

for i = 1:numel(recoFiles)

    recoFile = fullfile( ...
        recoFiles(i).folder, ...
        recoFiles(i).name);

    S = load(recoFile, 'reco');

    if ~isfield(S,'reco') || ...
       ~isfield(S.reco,'mni') || ...
       ~isfield(S.reco.mni,'coords_mm')

        fprintf('Skipping %s\n', recoFiles(i).name);
        continue
    end


    % Build the SEEG elstruct

    elstruct = struct;

    elstruct.coords_mm = S.reco.mni.coords_mm;


    % Electrode model information
    if isfield(S.reco,'props')
        elstruct.props = S.reco.props;
    end


    % Any stored electrode model field
    if isfield(S.reco,'elmodel')
        elstruct.elmodel = S.reco.elmodel;
    end


    figure(resultfig);
    hold on;

    fprintf('\nRendering %s\n', recoFiles(i).name);

    options.subjectColor = subjectColors(i,:);

    ea_trajectory_seeg(elstruct, options);

    fprintf('Plotted %s\n', recoFiles(i).name);

end

figure(resultfig);

axis equal;
axis vis3d;

fprintf('\nDone — plotted all SEEG reconstructions in one MNI figure.\n');

end