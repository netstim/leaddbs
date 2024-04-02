function settings = ea_switch2VATgrid(options,settings,side)

% change connectome fibers to regular grid aligned with the electrode
% for parameters, see ea_till_creategridforelectrode

reco = load(options.subj.recon.recon);
reco = reco.reco;
% we actually need to ensure that this 
if options.native
    ea_till_creategridforelectrode(reco,side+1,'scrf',options)
else
    ea_till_creategridforelectrode(reco,side+1,'mni',options)
end

settings.axonLength = [15;15;20];

% re-create data
% Get paths of tracts
connName = 'OSSDBSgrid';
connFolder = [options.subj.subjDir,filesep,'connectomes',filesep,'dMRI_MultiTract',filesep,connName];
tracts = ea_regexpdir(connFolder, '\.mat$', 0);
data1 = struct;
data2 = struct;
fibersFound = zeros(numel(tracts),2);
for t=1:numel(tracts)
    tract = tracts{t};
    [~, tractName] = fileparts(tract);
    settings.connectomeTractNames{t} = tractName;
    fprintf('Loading connectome: %s, Tract: %s ...\n', connName, tractName);
    conn = load(tract);

    % no need to convert
    % if options.native
    %     originalFib = conn;
    %     % Convert connectome fibers from MNI space to anchor space
    %     fprintf('Convert connectome into native space...\n\n');
    %     fibersMNIVox = ea_mm2vox(conn.fibers(:,1:3), [ea_space, options.primarytemplate, '.nii'])';
    %     conn.fibers(:,1:3)  = ea_map_coords(fibersMNIVox, ...
    %         [ea_space, options.primarytemplate, '.nii'], ...
    %         [options.subj.subjDir, filesep, 'forwardTransform'], ...
    %         preopAnchor)';
    % end

    % Filter fibers based on the spherical ROI
    if options.native
        fiberFiltered = ea_filterfiber_stim(conn, coords_mm, stimProtocol, 'kuncel', 1, preopAnchor);
    else
        fiberFiltered = ea_filterfiber_stim(conn, coords_mm, stimProtocol, 'kuncel', 1, [ea_space, options.primarytemplate, '.nii']);
    end

    % Filter fibers based on the minimal length
    fiberFiltered = ea_filterfiber_len(fiberFiltered, settings.axonLength(t));

    % Move original fiber id to the 5th column, the 4th column will be 1:N
    for i=1:length(fiberFiltered)
        if ~isempty(fiberFiltered{i}.fibers)
            fibers = zeros(size(fiberFiltered{i}.fibers,1),5);
            fibers(:,[1,2,3,5]) = fiberFiltered{i}.fibers;
            fibers(:,4) = repelem(1:length(fiberFiltered{i}.idx), fiberFiltered{i}.idx)';
            fiberFiltered{i}.fibers = fibers;

            % store the original number of fibers
            % to compute percent activation
            fiberFiltered{i}.origNum = size(conn.idx,1);

            fibersFound(t,i) = 1;
        end
    end
    data1.(tractName) = fiberFiltered{1};
    data2.(tractName) = fiberFiltered{2};
end

% Create output folder
settings.connectomePath = [outputDir, filesep, connName];
ea_mkdir(settings.connectomePath);
settings.connectome = ['Multi-Tract: ', connName];

% Save filtered fibers
save([settings.connectomePath, filesep, 'data1.mat'], '-struct', 'data1', '-v7.3');
save([settings.connectomePath, filesep, 'data2.mat'], '-struct', 'data2', '-v7.3');


%if settings.axonLength > 1
%    settings.axonLength = repmat(options.prefs.machine.vatsettings.butenko_axonLength(1),3,1);
%endI 
if settings.fiberDiameter > 1
    settings.fiberDiameter = repmat(options.prefs.machine.vatsettings.butenko_fiberDiameter(1),3,1);
end
