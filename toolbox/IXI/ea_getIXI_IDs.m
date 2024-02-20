function IDs=ea_getIXI_IDs(n,age,sd)
% This function returns n age-matched IDs from the IXI dataset.

load([ea_getearoot,'toolbox',filesep,'IXI',filesep,'IXI_demographics.mat']); % adds vars IXI.AGE and IXI.ID

% % excludeIDs with wrong warpings:
% [~,idx]=ismember(IXI.excl,IXI.uID);
% IXI.uID(idx)=[];
% IXI.uAGE(idx)=[];

if ~exist('n','var')
   IDs=IXI.uID;
   IDs=ea_convID2cell(IDs);
   return
end

if ~exist('age','var') % give n random IDs back
    numids=randperm(length(IXI.uID));
    IDs=IXI.uID(numids(1:n));
    IDs=ea_convID2cell(IDs);
    return
end


if ~exist('sd','var') % five n closest IDs to specified age back 
    agediff=abs(IXI.uAGE-age);
    [agediff,idx]=sort(agediff,'ascend');
    IDs=IXI.uID(idx(1:n));
    IDs=ea_convID2cell(IDs);
    return
end

ea_error('Standard-Deviation Matching not yet supported.');


function ID_cell=ea_convID2cell(IDs)
prefs=ea_prefs('');
try
basedir=prefs.ixi.dir;
catch
   ea_error('IXI database path not set. Please specify the path to the IXI database by using the variable prefs.ixi.dir in your prefs');
end

% check if basedir is mounted:
if ~exist(basedir,'file')
   ea_error('IXI database path is specified but the folder cannot be accessed. Please make sure that the absolute path in your prefs (prefs.ixi.dir) points to the correct folder containing the IXI database.');
end

for i=1:length(IDs)
ID_cell{i}=[basedir,'IXI',sprintf('%03.0f',IDs(i))];
end
