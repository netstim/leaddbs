% plotStims.m
%
%      usage: plotStims(stimOnsets, stimDurations, framePeriod, colors, axesHandle)
%         by: julien besle 
%       date: 02/12/2010
%    purpose: Plots stimulus matrix (in cell vector or matrix form). Onsets and durations must be integers and expressed in frames;
%              $Id$

function [hStims,hTransitions] = plotStims(stimOnsets, stimDurations, framePeriod, colors, axesHandle, runTransitions)


hTransitions = [];

if ieNotDefined('axesHandle')
  axesHandle = gca;
end
if ieNotDefined('framePeriod')
  framePeriod = 1;
end
if ieNotDefined('runTransitions')
  runTransitions = [];
end

%onsets and durations have to be the same size, otherwise durations will be ignored
if ieNotDefined('stimDurations') || any(size(stimDurations)~=size(stimOnsets))
  stimDurations=[];
end

if iscell(stimOnsets) %If the times and durations are in array of cell form, transform to matrix form
  stimOnsets = stimCell2Mat(stimOnsets, stimDurations,runTransitions);
end

if ieNotDefined('colors')
  colors = randomColors(size(stimOnsets,2));
end

if verLessThan('matlab','8.4')
  set(axesHandle,'drawMode','fast');
else
  set(axesHandle,'SortMethod','childorder');
end
hold(axesHandle,'on');
Ylims = get(axesHandle,'Ylim');
bottom = Ylims(1) - diff(Ylims);
top = Ylims(2) + diff(Ylims);
Ylims(1) = Ylims(1) - .2*diff(Ylims);
Ylims(2) = Ylims(2) + .2*diff(Ylims);
% bottom = Ylims(1);
% top = Ylims(2);

stimOnsets = logical(stimOnsets);
[patterns, dump, whichPatterns] = unique(stimOnsets,'rows'); 
% % % %remove empty pattern
% % % emptyPattern = find(all(patterns==0,2));
% % % if ~isempty(emptyPattern)
% % %   patterns(emptyPattern,:)=[];
% % %   whichPatterns(whichPatterns==emptyPattern)=0;
% % %   whichPatterns(whichPatterns>emptyPattern)=whichPatterns(whichPatterns>emptyPattern)-1;
% % % end  

%for each pattern
for iPattern = 1:size(patterns,1)
  if ~all(patterns(iPattern,:)==0)
    %find onsets and offsets for each block of this pattern
    blocks = reshape(find(diff([0;whichPatterns;0]==iPattern)),2,[])-1;
    %and draw each block
    for iBlock = 1:size(blocks,2)
      hStims(patterns(iPattern,:)') = plotFrame(axesHandle,framePeriod*blocks(1,iBlock), framePeriod*blocks(2,iBlock), bottom, top, colors(patterns(iPattern,:)',:));
    end
  end
end

if ~isempty(runTransitions)
  for iRun = 1:size(runTransitions,1)-1
    hTransitions = plot(axesHandle,[runTransitions(iRun,2) runTransitions(iRun,2)]*framePeriod,[bottom top],':k');%,'Color',[.3 .3 .3]);
  end
end
xlabel(axesHandle,'Time (s)');
set(axesHandle,'Xlim',[0 size(stimOnsets,1)*framePeriod]); 
set(axesHandle,'Ylim',Ylims);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% plotFrame %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = plotFrame(axesHandle,left,right,bottom,top,colors)

if size(colors,1)==1
   h = patch('vertices',[left top; right top; right bottom; left bottom], 'faces',[1 2 3 4]...
      ,'edgecolor','none','facecolor',colors,'parent',axesHandle... %,'faceAlpha',.5
      ); 
else
  numberColors = size(colors,1);
  h = zeros(1,numberColors);
  stripeWidth = (right-left)/numberColors;
  stripeHeight = (top-bottom)/numberColors;
  
  for iColor = 1:numberColors
     h(iColor) = patch('vertices',...
       [left+stripeWidth*(iColor-1) top; left+stripeWidth*iColor top; left top-stripeHeight*(iColor); left top-stripeHeight*(iColor-1)]...
       ,'faces',[1 2 3 4],'edgecolor','none','facecolor',colors(iColor,:),'parent',axesHandle... %,'faceAlpha',.5...
       );
     patch('vertices',...
       [left+stripeWidth*(iColor-1) bottom; left+stripeWidth*iColor bottom; right top-stripeHeight*(iColor); right top-stripeHeight*(iColor-1)]...
       ,'faces',[1 2 3 4],'edgecolor','none','facecolor',colors(iColor,:),'parent',axesHandle... %,'faceAlpha',.5...
       );
  end
end

% to remove

% % % if ~iscell(stimOnsets) %If the times and durations are in matrix form, transform to array of cell form
% % %   maxTime = size(stimOnsets,1)*framePeriod;
% % %   stimOnsetsMatrix = stimOnsets;
% % %   stimOnsets = cell(1,size(stimOnsetsMatrix,2));
% % %   for iStim  =1:size(stimOnsets,2)
% % %     stimOnsets{iStim} = find(stimOnsetsMatrix(:,iStim));
% % % %     stimOnsets{iStim} = diff([0 stimOnsetsMatrix(:,iStim)]); %could find onsets and duration to limit the 
% % % %     stimOnsets{iStim} = stimOnsets{iStim}(1:2:end);          %number of graphic object, but not worth the trouble
% % %   end
% % % else
% % %   maxTime=0;
% % % end
% % % if ieNotDefined('stimDurations')
% % %   for iStim  =1:size(stimOnsets,2)
% % %     stimDurations{iStim} = ones(size(stimOnsets{iStim}));
% % %   end
% % % end
% % % 
% % % if ieNotDefined('colors')
% % %   colors = randomColors(length(stimOnsets));
% % % end
% % % 
% % % 
% % % axes(axesHandle);
% % % hold on;
% % % Ylims = get(axesHandle,'Ylim');
% % % bottom = Ylims(1) - 10*diff(Ylims);
% % % top = Ylims(2) + 10*diff(Ylims);
% % % for iEventType = 1:length(stimOnsets)
% % %   if ~isempty(stimOnsets{iEventType})
% % %     
% % %     for iEvent = 1:length(stimOnsets{iEventType})
% % %       left = (stimOnsets{iEventType}(iEvent)-1)*framePeriod;
% % %       right = (stimOnsets{iEventType}(iEvent)-1+stimDurations{iEventType}(iEvent))*framePeriod;
% % %       h(iEventType) = patch('vertices',[left bottom; left top; right top; right bottom], 'faces',[1 2 3 4],...
% % %         'faceAlpha',.5,'edgecolor','none','facecolor',colors(iEventType,:));
% % %     end
% % %     maxTime = max(maxTime,max(right));
% % %   end
% % % end
% % % xlabel('Time (s)');
% % % set(axesHandle,'Xlim',[0 maxTime]);
% % % set(axesHandle,'Ylim',Ylims);

