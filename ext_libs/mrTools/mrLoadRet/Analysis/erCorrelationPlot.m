% erCorrelationPlot.m
%
%      usage: erCorrelationPlot()
%         by: justin gardner
%       date: 05/12/08
%    purpose: interrogator function that displays the two er's. i.e. the
%             one from which the correlation was made and the current voxel
%
function retval = erCorrelationPlot(v,overlayNum,scan,x,y,s,roi,varargin)


% get the analysis structure
analysis = viewGet(v,'analysis');
d = analysis.d{scan};
if isempty(d)
  disp(sprintf('(eventRelatedPlot) Event related not for scan %i',scan));
  return
end
d.r2 = analysis.overlays(1).data{scan};
rOverlay = viewGet(v,'overlay',viewGet(v,'overlayNum','r'));
d.r = rOverlay.data{scan};

% select the window to plot into
fignum = selectGraphWin;

% plot this one
subplot(1,3,2);
[ehdr time ehdrste] = gethdr(d,x,y,s);
plotEhdr(time,ehdr,ehdrste,'-');
title(sprintf('Original\n[%i %i %i] r2=%0.3f r=%0.3f',x,y,s,d.r2(x,y,s),d.r(x,y,s)));

% plot this one
subplot(1,3,3);
ehdr = reshape(d.projectedOutEhdr(x,y,s,:),d.nhdr,d.hdrlen);
plotEhdr(time,ehdr,[],'-');
title(sprintf('Projected out\n[%i %i %i]',x,y,s));

% plot the one this correlation is based on
subplot(1,3,1);
ehdr = d.projectionEhdr;
plotEhdr(time,ehdr,[],'-');
ylabel('Normalized response');
if isempty(rOverlay.params.roiName)
  title(sprintf('Seed voxel\n[%i %i %i] r2=%0.3f',rOverlay.params.x,rOverlay.params.y,rOverlay.params.s,d.r2(rOverlay.params.x,rOverlay.params.y,rOverlay.params.s)));
else
  title(sprintf('Seed roi: %s (n=%i)',rOverlay.params.roiName,rOverlay.params.roiN),'interpreter','none');
end

makeEqualYaxis(1,3);

%%%%%%%%%%%%%%%%%%%%%%%%%
% function to plot ehdr
%%%%%%%%%%%%%%%%%%%%%%%%%
function plotEhdr(time,ehdr,ehdrste,lineSymbol)

% whether to plot the line inbetween points or not
if ~exist('lineSymbol','var'),lineSymbol = '-';,end

% and display ehdr
for i = 1:size(ehdr,1)
  if ieNotDefined('ehdrste')
    h=plot(time,ehdr(i,:),getcolor(i,getsymbol(i,lineSymbol)),'MarkerSize',8);
  else
    h=errorbar(time,ehdr(i,:),ehdrste(i,:),ehdrste(i,:),getcolor(i,getsymbol(i,lineSymbol)),'MarkerSize',8);
  end
  set(h,'MarkerFaceColor',getcolor(i));
  hold on
end
xlabel('Time (sec)');
ylabel('% Signal change');

