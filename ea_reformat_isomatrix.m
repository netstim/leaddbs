function isom=ea_reformat_isomatrix(isom,M,options)
% this function tries to put the isomatrix provided into the correct format
% to be used by LEAD-DBS. The isomatrix is a matrix of values that can be
% displayed either at each active contact, each contact or in the spacings
% between contacts (e.g. for power-analyses).
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

stimmat = cell(1, length(options.sides));
if ~iscell(isom) % check if isomatrix is a cell ({[right_matrix]},{[left_matrix]}), if not convert to one.
    if size(isom,1) == length(M.patient.list) && size(isom,2) == 1
        % single vector (1 value for each patient)
        for iside=1:length(options.sides)
            side=options.sides(iside);
            stimmat{side}=init_isoMatrixMask(M.elstruct, side);
            stimmat{side}=bsxfun(@times,stimmat{side},isom);
        end
    elseif size(isom,1) == length(M.patient.list) && size(isom,2) == length(options.sides)
        % nx2 matrix (1 column for each hemisphere)
        for iside=1:length(options.sides)
            side=options.sides(iside);
            stimmat{side}=init_isoMatrixMask(M.elstruct, side);
            stimmat{side}=bsxfun(@times,stimmat{side},isom(:,side));
        end
    elseif size(isom,1) == length(M.patient.list) && ...
            (size(isom,2) == (get_maxNumContacts(M.elstruct)-1)*2 || size(isom,2) == get_maxNumContacts(M.elstruct)*2)
        % (1 value for each contact pair) or (1 value for each contact) * patientlist
        stimmat{1}=isom(:,1:size(isom,2)/2);
        stimmat{2}=isom(:,(size(isom,2)/2)+1:end);
    else
        ea_error('Wrong Isomatrix provided.');
    end
end

if options.normregressor>1 % apply normalization to regressor data
    for iside=1:length(options.sides)
        side=options.sides(iside);
        if options.normregressor==2 % apply z-score
            stimmat{side}=reshape(nanzscore(stimmat{side}(:)),size(stimmat{side},1),size(stimmat{side},2));
        elseif options.normregressor==3 % apply normal method from van albada 2008
            nanidx=isnan(stimmat{side});

            stimmat{side}=reshape(ea_normal(stimmat{side}(:)),size(stimmat{side},1),size(stimmat{side},2));
            stimmat{side}(nanidx)=nan;
        end
    end
end

isom=stimmat;


function z=nanzscore(data)

datawonan = data(~isnan(data));
datamean = mean(datawonan);
datasd = std(datawonan);
z = (data-datamean)/datasd;


function maxNumContacts = get_maxNumContacts(elstruct, side)
if ~exist('side', 'var') % Check both sides
    coords = {elstruct.coords_mm};
    coords = horzcat(coords{:})';
    maxNumContacts = max(cellfun(@(x) size(x,1), coords));
    return;
end

coords = {elstruct.coords_mm}';

% Get max number of contacts for the group of patients
maxNumContacts = 0;
for p = 1:length(coords)
    for s = side
        if size(coords{p}{side},1) > maxNumContacts
            maxNumContacts = size(coords{p}{side},1);
        end
    end
end


function isoMatrixMask = init_isoMatrixMask(elstruct, side)

coords = {elstruct.coords_mm}';

% Get max number of contacts for the group of patients
maxNumContacts = get_maxNumContacts(elstruct, side);

isoMatrixMask = nan(length(coords), maxNumContacts);
for p = 1:length(coords)
    isoMatrixMask(p, 1:size(coords{p}{side},1)) = 1;
end
