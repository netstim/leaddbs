function [fibidx,peakvals,meanvals,sumvals,peak5vals,binvals] = ea_explorer_efoverlapfiber(vol,img,fibers)
fibers = single(fibers);
if ~isequal(sort(fibers(:,4)),fibers(:,4))
    keyboard
end
%% convert fibers to vx-coordinates
fibersorg = fibers;
fibers = vertcat(fibers(:,1:3)',ones(1,size(fibers,1)));
fibers = round(vol.mat \ fibers);
fibers = fibers';
fibers(:,4) = fibersorg(:,4);
clear fibersorg
%% Remove values outside the template
for dim = 1:3
    fibers(fibers(:,dim)<1,:)=[];
    fibers(fibers(:,dim)>vol.dim(dim),:)=[];
end
%% get voxelvals for each fiber coordinate by first getting linear indexes
linidx = sub2ind(size(img),fibers(:,1),fibers(:,2),fibers(:,3));
fibvals = img(linidx);
%% sort fibers and their corresponding values into cell array with one cell per fiber
idxnew = diff(find([true,diff(fibers(:,4)')~=0,true]))';
fibercell=mat2cell(fibers(:,4),idxnew);
valscell= mat2cell(fibvals,idxnew);
%% calculate maximum value for each cell
peakwrapper = @(x) max(x,[],'omitnan');
peakvals = cellfun(peakwrapper,valscell);
meanwrapper = @(x) mean(x,'all','omitnan');
meanvals = cellfun(meanwrapper,valscell);
sumwrapper = @(x) sum(x,'all','omitnan');
sumvals = cellfun(sumwrapper,valscell);
peak5vals=cellfun(@(x) mean(maxk(x,ceil(0.05*numel(x)))), valscell);
fibidx = unique(fibers(:,4));
binvals = ones(size(fibidx));