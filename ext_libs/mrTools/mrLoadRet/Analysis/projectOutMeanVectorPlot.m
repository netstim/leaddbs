% projectOutMeanVectorPlot.m
%
%        $Id$
%      usage: projectOutMeanVectorPlot()
%         by: justin gardner
%       date: 05/02/08
%    purpose: 
%
function projectOutMeanVectorPlot(view,overlayNum,scan,x,y,s,roi)

% check arguments
if ~any(nargin == [1:7])
  help projectOutMeanVectorPlot
  return
end

% get the analysis structure
analysis = viewGet(view,'analysis');
d = analysis.d{scan};
if isempty(d)
  disp(sprintf('(projectOutMeanVectorPlot) projectionAnal not for scan %i',scan));
  return
end
frameNums = [];
% get overlay num and name
overlayNum = viewGet(view,'curOverlay');
overlayName = viewGet(view,'overlayName');

% figure out which d it corresponds to.
if strcmp(overlayName,'r')
  dnums = 1;
  frameNums{1} = [];
  name = '';
elseif ~isempty(regexp(overlayName,'[0-9]+'))
  numloc = regexp(overlayName,'[0-9]+');
  dnums = str2num(overlayName(numloc:end));
  concatInfo = viewGet(view,'concatInfo',scan);
  frameNums{1} = concatInfo.runTransition(dnums(1),:);
  name = sprintf('(Scan %i)',dnums(1));
elseif strcmp(overlayName,'mean_r');
  concatInfo = viewGet(view,'concatInfo',scan);
  dnums = 1:length(d);
  for dnum = 1:length(dnums)
    frameNums{dnum} = concatInfo.runTransition(dnums(dnum),:);
  end
  name = '(mean)';
else
  disp(sprintf('(projectOutMeanVectorPlot) Unknown overlay name %s',overlayName));
  return
end

% check to make sure overlay exists
if isempty(analysis.overlays(overlayNum).data{scan})
  disp(sprintf('(projectOutMeanVectorPlot) Overlay does not exist for scan %i',scan));
  return
end

% look for this voxel in projection info
linearCoord = sub2ind(viewGet(view,'scanDims'),x,y,s);
thisLinearCoord = find(d{dnums(1)}.linearCoords == linearCoord);

% get projection magnitude
normMag = analysis.overlays(overlayNum).data{scan}(x,y,s);

% get the tSeries, reconMag and the sourceMeanVector
for dnum = 1:length(dnums)
  % get tSeries
  tSeries(dnum,:) = squeeze(loadTSeries(view,scan,s,frameNums{dnum},x,y))';

  % de mean normalize the tseries
  tSeries(dnum,:) = tSeries(dnum,:)-mean(tSeries(dnum,:));
  tSeries(dnum,:) = tSeries(dnum,:)/norm(tSeries(dnum,:));

  % get recon magnitude and source mean vector
  if ~isempty(thisLinearCoord)
    reconMag(dnum) = d{dnums(dnum)}.reconProjectionMagnitude(thisLinearCoord);
  end
  sourceMeanVector(dnum,:) = d{dnums(dnum)}.sourceMeanVector;
end

% take average if we are showing mean projection
if length(dnums) > 1
  % keep full tSeries for possible event-related processing below
  full.tSeries = reshape(tSeries',1,prod(size(tSeries)));
  full.sourceMeanVector = reshape(sourceMeanVector',1,prod(size(sourceMeanVector)));
  if ~isempty(thisLinearCoord)
    full.original = full.tSeries + reshape(sourceMeanVector'*diag(reconMag),1,prod(size(sourceMeanVector)));
    full.original = full.original/norm(full.original);
  else
    full.original = nan(size(full.tSeries));
  end
  % get mean of tSeries/reconMag and sourceMeanVector
  tSeries = mean(tSeries);
  if ~isempty(thisLinearCoord)
    reconMag = mean(reconMag);
  end
  sourceMeanVector = mean(sourceMeanVector);
end

% this voxel didn't have anything projected out
if isempty(thisLinearCoord)
  reconMag = 0;
end

% compute original
original = tSeries+reconMag*sourceMeanVector;
original = original/norm(original);

% compute what is left of tSeries along projection vector
leftoverMag = tSeries*sourceMeanVector';

% now check for erAnal
analysisTypes = viewGet(view,'analysisTypes');
erAnalNums = find(strcmp('erAnal',analysisTypes)||strcmp('deconvAnal',analysisTypes));
if isempty(erAnalNums) || (length(dnums)==1)
  nrows = 2;
  ncols = 3;
  doErAnal = 0;
else
  doErAnal = 1;
  nrows = 3;
  ncols = 3;
end
% select the window to plot into
selectGraphWin;

% plot the tSeries
subplot(nrows,ncols,1:ncols);
plot(tSeries,'k.-');hold on
if reconMag == 0
  plot(sourceMeanVector,'r.-');
  legend(sprintf('tSeries %s',name),sprintf('Projection vector %s',name));
else
  plot(sourceMeanVector*reconMag,'r.-');
  plot(original,'g.-');
  legend(sprintf('tSeries %s',name),sprintf('Projection vector %s',name),sprintf('tSeries before projection %s',name));
end

xlabel('Time (TR)');
ylabel('Normalized magnitude');
title(sprintf('%s Voxel: [%i %i %i] r: %0.3f Projection ROI: %s\nMagnitude of component left in projection direction: %f',name,x,y,s,normMag,d{dnums(1)}.sourceName,leftoverMag),'Interpreter','none');
axis tight

subplot(nrows,ncols,(1:ncols)+ncols);
plot(abs(fft(tSeries)),'k.-');hold on
plot(abs(fft(sourceMeanVector)),'r.-');
if reconMag ~= 0
  plot(abs(fft(original)),'g.-');
end
xaxis(0,size(tSeries,2)/2);
xlabel('Frequency (Fourier component number)')
ylabel('Magnitude');

drawnow;

if doErAnal
  % get the event related analysis
  % take first one, if there are multiples.
  erAnal = viewGet(view,'analysis',erAnalNums(1));
  % now get the parameters for event related processing
  erAnalD = erAnal.d{scan};
  if ~isempty(erAnalD)
    scm = erAnalD.scm;
    nhdr = erAnalD.nhdr;
    hdrlen = erAnalD.hdrlen;
    tr = erAnalD.tr;
  end

  % compute evented related responses
  full.er = getr2timecourse([full.tSeries;full.original;full.sourceMeanVector]+1,nhdr,hdrlen,scm,tr);

  % and display ehdr
  for tnum = 1:3
    subplot(nrows,ncols,(nrows-1)*ncols+tnum);
    for rnum = 1:nhdr
      h=errorbar(full.er.time,squeeze(full.er.ehdr(tnum,rnum,:)),squeeze(full.er.ehdrste(tnum,rnum,:)),squeeze(full.er.ehdrste(tnum,rnum,:)),getcolor(rnum,getsymbol(rnum,'-')),'MarkerSize',8);
      set(h,'MarkerFaceColor',getcolor(rnum));
      hold on
    end
    a = axis;
  end

  % set some titles
  subplot(nrows,ncols,(nrows-1)*ncols+1);
  ylabel(sprintf('Normalized Signal change.\n*NOT* percent signal change'));
  xaxis(0,max(full.er.time)+tr/2);
  title('tSeries');
  subplot(nrows,ncols,(nrows-1)*ncols+2);
  xaxis(0,max(full.er.time)+tr/2);
  xlabel('Time (sec)');
  title('Original');
  subplot(nrows,ncols,(nrows-1)*ncols+3);
  xaxis(0,max(full.er.time)+tr/2);
  title('Projection tSeries');

end