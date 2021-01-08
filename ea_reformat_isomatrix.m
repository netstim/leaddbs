function isom=ea_reformat_isomatrix(isom,M,options)
% this function tries to put the isomatrix provided into the correct format
% to be used by LEAD-DBS. The isomatrix is a matrix of values that can be
% displayed either at each active contact, each contact or in the spacings
% between contacts (e.g. for power-analyses).
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ~iscell(isom) % check if isomatrix is a cell ({[right_matrix]},{[left_matrix]}), if not convert to one.
    if min(size(isom))==1 && length(size(isom))==2 % single vector (1 value for each patient)
        for iside=1:length(options.sides)
            side=options.sides(iside);
            try
                stimmat{side}=cat(1,M.stimparams(:,1).U);
            catch
                stimmat{side}=ones(length(M.patient.list),4);
            end
            stimmat{side}=bsxfun(@times,stimmat{side}>0,isom);
        end
    elseif min(size(isom))==2 && length(size(isom))==2 % 2xn matrix (1 value for each hemisphere)
        for iside=1:length(options.sides)
            side=options.sides(iside);
            try
                stimmat{side}=cat(1,M.stimparams(:,1).U);
            catch
                stimmat{side}=ones(length(M.patient.list),4);
            end
            stimmat{side}=bsxfun(@times,stimmat{side}>0,isom(:,side));
        end
    elseif (min(size(isom))==6 || min(size(isom))==8) && max(size(isom))==length(M.patient.list) % 6 (1 value for each contact pair) or 8 (1 value for each contact) * patientlist
        if size(isom,2)==length(M.patient.list)
            isom=isom';
        end

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
