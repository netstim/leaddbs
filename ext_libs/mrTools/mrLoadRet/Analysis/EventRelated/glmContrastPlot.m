function glmContrastPlot(view,overlayNum,scan,x,y,s,roi)
% glmContrastPlot.m
%
%      usage: glmContrastPlot()
%         by: farshad moradi
%       date: 09/14/07
%    purpose: 
%

% check arguments
if ~any(nargin == [1:7])
  help glmContrastPlot
  return
end

% get the analysis structure
analysis = viewGet(view,'analysis');
d = analysis.d{scan};
if isempty(d)
  disp(sprintf('(glmContrastPlot) Event related not for scan %i',scan));
  return
end
d.r2 = analysis.overlays(1).data{scan};
% select the window to plot into
selectGraphWin;

% check to see if there is a regular event related analysis
for anum = 1:viewGet(view,'nAnalyses')
  if strcmp(viewGet(view,'analysisType',anum),'erAnal')
    % get the event related headers
    erAnalysis = viewGet(view,'analysis',anum);
    der = erAnalysis.d{scan};
    clear erAnalysis
  end
end

global MLR;
fignum = MLR.graphFigure;

% turn off menu/title etc.
set(fignum,'NumberTitle','off');
set(fignum,'Name','glmContrastPlot');

% set roi coords
for roinum = 1:length(roi)
  % get scan coordinates
  roi{end}.scanCoords = getROICoordinates(view,roi{roinum},scan);
end

% get cutoff value
cutoffr2 = viewGet(view,'overlayMin');

if isempty(d)
  disp('No analysis');
  return
end


% reconstruct what the glm produces
if ~isempty(analysis.params.scanParams{scan}) && isfield(analysis.params.scanParams{scan},'segmentBegin')
  % first want to get stimvols for just the beginning of the 
  % stimulus so that we can run an event-related analysis
  % on the glm modeled timecourse
  %
  % so, construct an appropriate params variable
  params = analysis.params.scanParams{scan};
  params.segmentNum = params.segmentBegin;
  params = rmfield(params,'segmentEnd');
  % get the stimvols 
  dglm = loadScan(view,scan,viewGet(view,'curGroup'),0); 
  dglm = getStimvol(dglm,params);
  % make stimulus convolution matrix
  dglm = makescm(dglm);
  % check to make sure dimensions match
  if size(d.scm,2) ~= size(d.ehdr,4)
    % no match, give up
    [ehdr time ehdrste] = gethdr(d,x,y,s);
    ehdrTitle = 'Scaled HRF';
  else
    % now we want to create a time course that is the glm design
    % multiplied by the beta weights (this is the estimated time
    % course)
    estTimecourse = d.scm*squeeze(d.ehdr(x,y,s,:));
    % now deconvolve to get the stimulus triggered averages
    % of the estimated timecourse
    estResponses = pinv(dglm.scm)*(estTimecourse-mean(estTimecourse));
    % and reshape
    ehdr = reshape(estResponses,dglm.hdrlen,dglm.nhdr)';
    time = dglm.tr/2:dglm.tr:(dglm.hdrlen*dglm.tr);
    ehdrste = [];
    ehdrTitle = 'GLM modeled response';
  end
else
  % get the estimated hdr from the beta weight times
  % the hdr function
  [ehdr time ehdrste] = gethdr(d,x,y,s);
  ehdrTitle = 'Scaled HRF';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the hemodynamic response for voxel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,2,1);
plotEcontrast(shiftdim(d.ehdr(x,y,s,:,:), 3),shiftdim(d.ehdrste(x,y,s,:,:), 3));

subplot(3,2,3);

% display ehdr with out lines if we have a fit
% since we also need to plot fit
if isfield(d,'peak') & isfield(d.peak,'fit') & ~any(isnan(d.peak.amp(x,y,s,:)))
  plotEhdr(time,ehdr,ehdrste,'');
  for r = 1:d.nhdr
    d.peak.fit{x,y,s,r}.smoothX = 1:.1:d.hdrlen;
    fitTime = d.tr*(d.peak.fit{x,y,s,r}.smoothX-0.5);
    plot(fitTime+d.tr/2,d.peak.fit{x,y,s,r}.smoothFit,getcolor(r,'-'));
  end
else
  % if there is deconvolution data, display that too
  if exist('der','var') && (length(der.stimvol) == length(d.stimvol))

    [deconvEhdr deconvTime deconvEhdrste] = gethdr(der,x,y,s);
    % note that we need to subtract the mean of the glm ehdr
    % to match the event related and the glm data. This is
    % because the glm has the mean subtracted from the columns
    % i.e. the glm gives an estimate of the response *before*
    % mean subtraction
    if size(ehdr,1) == size(deconvEhdr,1)
      deconvEhdr = deconvEhdr-repmat(mean(ehdr,2),1,size(deconvEhdr,2));
    else
      deconvEhdr = deconvEhdr-repmat(mean(ehdr(:)),size(deconvEhdr,1),size(deconvEhdr,2));
    end      

    % chop down to the same size as the glm responses
    displayLen = min(size(deconvEhdr,2),size(ehdr,2));
    ehdr = ehdr(:,1:displayLen);
    time = time(1:displayLen);
    % now display
    plotEhdr(time,ehdr,ehdrste,'-',1);
    subplot(3,2,5);
    plotEhdr(deconvTime,deconvEhdr,deconvEhdrste,'-');
    title('Deconvolution');
    legend(der.stimNames);
    subplot(3,2,3);
  else
    plotEhdr(time,ehdr,ehdrste);
  end
end
title(sprintf('%s: (%i,%i,%i): r2=%0.3f',ehdrTitle,x,y,s,analysis.overlays(1).data{scan}(x,y,s)));
xaxis(0,max(time));
% add peaks if they exist to the legend
if isfield(d,'stimNames')
  stimNames = d.stimNames;
  if isfield(d,'peak')
    for i = 1:d.nhdr
      stimNames{i} = sprintf('%s: %s=%0.2f',stimNames{i},d.peak.params.method,d.peak.amp(x,y,s,i));
    end
  end
  legend(stimNames);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if there is an roi at this voxel
% then plot mean response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for roinum = 1:length(roi)
  subplot(3,2,4);
  ehdr = [];
  econt = [];
  roin = 0;
  for voxnum = 1:size(roi{roinum}.scanCoords,2)
    x = roi{roinum}.scanCoords(1,voxnum);
    y = roi{roinum}.scanCoords(2,voxnum);
    s = roi{roinum}.scanCoords(3,voxnum);
    if d.r2(x,y,s) >= cutoffr2
      roin = roin+1;
      [ehdr(roin,:,:) time] = gethdr(d,x,y,s);
      econt(roin,:,:)  = shiftdim(d.ehdr(x,y,s,:,:), 3);
      % if there is a peak field, calculate average peak
      if isfield(d,'peak')
	for i = 1:d.nhdr
	  amp(i,roin) = d.peak.amp(x,y,s,i);
	end
      end
    end
  end
  % plot the average of the ehdrs that beat the r2 cutoff
  if roin
    subplot(3,2,4);
    plotEhdr(time,shiftdim(mean(ehdr),1),shiftdim(std(ehdr),1)/sqrt(size(roi{roinum}.scanCoords,2)));
    subplot(3,2,2);
    plotEcontrast(shiftdim(mean(econt),1),shiftdim(std(econt),1)/sqrt(size(roi{roinum}.scanCoords,2)));
  end
  title(sprintf('%s (n=%i/%i)',roi{roinum}.name,roin,size(roi{roinum}.scanCoords,2)),'Interpreter','none');
  % create a legend (only if peaks exist) to display mean amplitudes
  if isfield(d,'peak')
    for i = 1:d.nhdr
      % get the stimulus name
      if isfield(d,'stimNames')
	stimNames{i} = d.stimNames{i};
      else
	stimNames{i} = '';
      end
      % and now append the peak info
      stimNames{i} = sprintf('%s: median=%0.2f',stimNames{i},median(amp(i,:)));
    end
    legend(stimNames);
  end
end

drawnow;


%%%%%%%%%%%%%%%%%%%%%%%%%
% function to plot ehdr
%%%%%%%%%%%%%%%%%%%%%%%%%
function plotEhdr(time,ehdr,ehdrste,lineSymbol,drawSymbols)

% whether to plot the line inbetween points or not
if ~exist('lineSymbol','var'),lineSymbol = '-';,end
if ~exist('drawSymbols','var'),drawSymbols = 1;,end

% and display ehdr
for i = 1:size(ehdr,1)
  if (nargin == 2) || isempty(ehdrste)
    % draw with symbols, no error bars
    if drawSymbols
      h=plot(time,ehdr(i,:),getcolor(i,getsymbol(i,lineSymbol)),'MarkerSize',8);
    % draw without symbols, no error bars
    else
      h=plot(time,ehdr(i,:),getcolor(i,lineSymbol));
    end
  else
    % draw with symbols, and error bars
    if drawSymbols
      h=errorbar(time,ehdr(i,:),ehdrste(i,:),ehdrste(i,:),getcolor(i,getsymbol(i,lineSymbol)),'MarkerSize',8);
    % draw without symbols, and error bars
    else
      h=errorbar(time,ehdr(i,:),ehdrste(i,:),ehdrste(i,:),getcolor(i,lineSymbol));
    end
  end
  set(h,'MarkerFaceColor',getcolor(i));
  hold on
end
xlabel('Time (sec)');
ylabel('% Signal change');

function plotEcontrast(econt,econtste)

% display econt
cla;
h = bar(econt');
if size(econt,2)==1
    econt = econt';
    econtste = econtste';
    xaxisLabel = 'contrast';
else
    xaxisLabel = 'hrf component';
end
if nargin==2
    hold on;
    for i=1:length(h)
        % location of the bar
        x = get(get(h(i),'Children'),'XData');
        % find the center of the bar
        x = (x(2,:)+x(3,:))/2;
        herr = errorbar(x, econt(i,:), econtste(i,:), 'k');
        temp = get(herr, 'Children');
        set(temp(1), 'visible', 'off');
    end
end
xaxis(0.5,length(econt)+0.5);
xlabel(xaxisLabel);
ylabel('beta');


