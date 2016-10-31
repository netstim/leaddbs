% eventRelatedMultiple.m
%
%      usage: view = eventRelatedMultiple(view,params)
%         by: farshad moradi (based on code by justin gardner)
%       date: 02/05/07
%    purpose: eventrelated for multiple files
%
function [view d] = eventRelatedMultiple(view,params)

% check arguments
if ~any(nargin == [1 2])
  help(mfilename)
  return
end
d = [];

% First get parameters
if ieNotDefined('params')
  % Initialize analysis parameters with default values
  params = eventRelatedGUI('groupName',viewGet(view,'groupName'));
end

% Reconcile params with current status of group and ensure that it has
% the required fields. 
params = mrParamsReconcile([],params);

% Abort if params empty
if ieNotDefined('params')
  mrMsgBox('eventRelated cancelled');
  return
end

if ~isfield(params, 'includeScans')
    params.includeScans = params.scanNum;
end;

if length(params.includeScans)==1
    % pass
elseif ~params.scanParams{1}.sameForAll
    mrMsgBox('eventRelated with inplace concat requires same parameters for all scans');
    return
end

params.hdrlen = params.scanParams{1}.hdrlen;


% step 1: estimate hrf using all included scans
ehdr = 0;
meanintensity = 0;
[fullscm, scms, volumes, nhdr] = getFullDesign(view, params);
if nhdr<=0
    mrMsgBox('incompatible number of conditions');
    return
end

set(viewGet(view,'figNum'),'Pointer','watch');drawnow;
for i=params.includeScans
    % load the scan
    d = loadScan(view,i);
    d.scms = scms;
    d.volumes = volumes{i};
    % do any called for preprocessing
    d = eventRelatedPreProcess(d,params.scanParams{i}.preprocess);    
    % compute the hemodynamic responses
    d.scanNum = i;
    d = calcleastsqestimate(d);
    ehdr = ehdr+d.ehdr;
    meanintensity = meanintensity+d.meanintensity;
end

meanintensity = meanintensity/length(params.includeScans);

% step 2: calculate r2
unexplainedVariance = 0;
totalVariance = 0;

for i=params.includeScans
    d = loadScan(view,i);
    d.ehdr = ehdr;
    d.scm = scms{i};
    d.volumes = volumes{i};
    % do any called for preprocessing
    d = eventRelatedPreProcess(d,params.scanParams{i}.preprocess);    
    % compute the residual error    
    d = calcVariances(d);
    unexplainedVariance = unexplainedVariance+d.unexplainedVariance;
    totalVariance = totalVariance+d.totalVariance;
end

d = rmfield(d,'data');
d.nhdr = nhdr;
d.hdrlen = ceil(params.hdrlen/d.tr);

ehdr = reshape(ehdr, [d.dim(1:3), d.hdrlen, nhdr]);
ehdr = permute(ehdr, [1,2,3,5,4]);

d.r2 = 1-unexplainedVariance./totalVariance;
d.meanintensity = meanintensity;
d.ehdr = ehdr./repmat(meanintensity, [1,1,1,nhdr, d.hdrlen])*100;
S2 = unexplainedVariance/(size(fullscm,1)-size(fullscm,2));
% now distribute that error to each one of the points
% in the hemodynamic response according to the inverse
% of the covariance of the stimulus convolution matrix.

ehdrste = sqrt(S2(:)*diag(pinv(fullscm'*fullscm))');
ehdrste = reshape(ehdrste, [d.dim(1:3), d.hdrlen, nhdr]);
ehdrste = permute(ehdrste, [1,2,3,5,4]);
d.ehdrste = ehdrste./repmat(meanintensity, [1,1,1,nhdr, d.hdrlen])*100;

% create the r2 overlay
dateString = datestr(now);
r2.name = 'r2';
r2.groupName = params.groupName;
r2.function = 'eventRelated';
r2.reconcileFunction = 'mrParamsReconcile';
r2.data = cell(1,viewGet(view,'nScans'));
for i=params.includeScans(1)
    r2.data{i}=d.r2;
end;
r2.date = dateString;
% r2.params = params;
r2.range = [0 1];
r2.clip = [0 1];
r2.params = cell(1,viewGet(view,'nScans'));
for i=params.includeScans
    r2.params{i} = params.scanParams{i};
    r2.params{i}.scanNum = i;
end;
% colormap is made with a little bit less on the dark end
r2.colormap = hot(312);
r2.colormap = r2.colormap(end-255:end,:);
r2.alpha = 1;
r2.colormapType = 'setRangeToMax';
r2.interrogator = 'eventRelatedPlot';


% install analysis
erAnal.name = 'erAnal';  % This can be reset by editAnalysisGUI
erAnal.type = 'erAnal';
erAnal.groupName = params.groupName;
erAnal.function = 'eventRelated';
erAnal.reconcileFunction = 'eventRelatedReconcileParams';
erAnal.guiFunction = 'eventRelatedGUI';
erAnal.params = params;
erAnal.overlays = r2;
erAnal.curOverlay = 1;
erAnal.date = dateString;
erAnal.ehdr = d.ehdr;
erAnal.nhdr = nhdr;
erAnal.d = {d};

view = viewSet(view,'newAnalysis',erAnal);

% Save it
if ~isfield(params, 'suppress_save')
    saveAnalysis(view,erAnal.name);
end
set(viewGet(view,'figNum'),'Pointer','arrow');drawnow


%%%%%%%%%%%%%%%%%%%
% getFullDesign
%%%%%%%%%%%%%%%%%%%
function [fullscm, scms, volumes, nhdr] = getFullDesign(view, params)
scms=[];
volumes = [];
nhdr = [];
accumLength = 1;

d.stimfile = cell(1, max(params.includeScans));
d.concatInfo.runTransition = zeros(max(params.includeScans), 2);
d.tr = [];

for i=params.includeScans
    stimfile = viewGet(view,'stimfile',i);
    if length(stimfile) == 1
        d.stimfile{i} = stimfile{1};
    end;
    nFrames = viewGet(view,'nFrames',i);
    d.concatInfo.runTransition(i,:) = accumLength+[0 nFrames-1];
    accumLength = accumLength+nFrames;
    if isempty(d.tr)
        d.tr = viewGet(view,'framePeriod',i);
    elseif d.tr ~= viewGet(view,'framePeriod',i)
        disp('Incompatible TRs');
        retrun;
    end;
    d.dim = viewGet(view,'scandims',i);
    volumes{i}=1:nFrames;
end;
% get the stim volumes, if empty then abort
d = getStimvol(d);
if isempty(d.stimvol),mrWarnDlg('No stim volumes found');return,end
nhdr=length(d.stimvol);
d.dim(4) = accumLength;
% make a stimulation convolution matrix
d = makescm(d,ceil(params.hdrlen/d.tr));
fullscm = d.scm;
for i=params.includeScans
    if d.concatInfo.runTransition(i,1)>0
        scms{i}=d.scm(d.concatInfo.runTransition(i,1):d.concatInfo.runTransition(i,2),:);
    end;
end;


%%%%%%%%%%%%%%%%%%%
% calcVariances
%
%      usage: d = calcVariances(d)
%         by: farshad moradi
%       date: 02/05/07
%    purpose: calculates unexplained and total variance
%
%%%%%%%%%%%%%%%%%%%
function d = calcVariances(d)

% no roi, so just do whole volume
slices = 1:d.dim(3);
yvals = 1:d.dim(2);
yvaln = length(yvals);
% initialize values
uv = [];
tv = [];
% preallocate memory
d.unexplainedVariance = zeros(d.dim(1),d.dim(2),d.dim(3));
d.totalVariance = zeros(d.dim(1),d.dim(2),d.dim(3));
% display string
disppercent(-inf,'Calculating goodness of fit');
% cycle through images calculating the estimated hdr and r^s of the 
% estimate.
%
onesmatrix = ones(length(d.volumes),1);
for j = yvals
    disppercent(max((j-min(yvals))/yvaln,0.1));
    for k = slices
        ehdr = squeeze(d.ehdr(:,j,k,:))';
        % get the time series we are working on
        timeseries = squeeze(d.data(:,j,k,d.volumes))';
        % subtract off column means
        colmeans = mean(timeseries,1);
        timeseries = timeseries - onesmatrix*colmeans;
        % calculate variance accounted for by the estimated hdr
        uv{j,k} = sum((timeseries-d.scm*ehdr).^2);
        tv{j,k} = sum(timeseries.^2);
    end
end
disppercent(inf);
% reshape matrix. 
for j = yvals
  for k = slices
    % now reshape into a matrix
    d.unexplainedVariance(:,j,k) = uv{j,k};
    d.totalVariance(:,j,k) = tv{j,k};
  end
end

%%%%%%%%%%%%%%%%%%%
% calcleastsqestimate
%
%      usage: d = calcleastsqestimate(d)
%         by: farshad moradi
%       date: 02/05/07
%    purpose: estimates hrf
%
%%%%%%%%%%%%%%%%%%%
function d = calcleastsqestimate(d)

% init some variables
ehdr=[];

% precalculate the normal equation (this dramatically speeds up things)
% Use the full design
n=zeros(1, 1+length(d.scms));
scm = [];
for i=1:length(d.scms),
  n(i+1)=size(d.scms{i}, 1);
  scm = [scm; d.scms{i}];
end;
n = cumsum(n);
precalcmatrix = pinv(scm);
precalcmatrix = precalcmatrix(:, n(d.scanNum)+1:n(d.scanNum+1));
% no roi, so just do whole volume
slices = 1:d.dim(3);
xvals = 1:d.dim(1);xvaln = length(xvals);
yvals = 1:d.dim(2);yvaln = length(yvals);
  
% preallocate memory
d.ehdr = zeros(d.dim(1),d.dim(2),d.dim(3),size(precalcmatrix,1));
d.meanintensity = zeros(d.dim(1),d.dim(2),d.dim(3));

% display string
disppercent(-inf,'Calculating hdr');
% cycle through images calculating the estimated hdr and r^s of the 
% estimate.
onesmatrix = ones(length(d.volumes),1);
for j = yvals
  disppercent(max((j-min(yvals))/yvaln,0.1));
  for k = slices
    % get the time series we are working on
    % this includes all the rows of one column from one slice
    % and all data points for each of these
    % thus the time series is a nxm matrix where each of the m columns
    % contains the n time points recording for that voxel
    timeseries = squeeze(d.data(:,j,k,d.volumes))';
    % subtract off column means
    colmeans = mean(timeseries,1);
    timeseries = timeseries - onesmatrix*colmeans;
    % get hdr for the each voxel
    ehdr{j,k} = precalcmatrix*timeseries;
    d.meanintensity(:,j,k)=colmeans(:);
  end
end
disppercent(inf);
% reshape matrix. this also seems the fastest way to do things. we
% could have made a matrix in the above code and then reshaped here
% but the reallocs needed to continually add space to the matrix
% seems to be slower than the loops needed here to reconstruct
% the matrix from the {} arrays.
disppercent(-inf,'Reshaping matrices');
for i = xvals
    disppercent((i-min(xvals))/xvaln);
    for j = yvals
        for k = slices
            % now reshape into a matrix
            d.ehdr(i,j,k,:) = ehdr{j,k}(:,i);
        end
    end
end

% display time took
disppercent(inf);

