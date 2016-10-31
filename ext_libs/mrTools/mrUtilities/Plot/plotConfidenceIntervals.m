% plotConfidenceIntervals - plot confidence intervals computed by roisConfidenceIntervals
%
%        $Id$
%      usage: [  ] = roisConfidenceInterval(actualValue,bootstrapMean,lowerErrorBar,upperErrorBar,axesValues,axesUnits,axesLabels,<figNum>,<plotMode>)
%         by: julien besle
%       date: 2010-10-23
%     inputs: 
%    outputs: 
%
%    taken out of roisConfidenceIntervals so that I can use it outside mrLoadRet
%


function plotConfidenceIntervals(actualValue,bootstrapMean,lowerErrorBar,upperErrorBar,statisticalThreshold,axesNames,axesValues,axesUnits,axesLabels,figNum,plotMode)

  if ~ismember(nargin,[9 10 11])
    help plotConfidenceIntervals
  end

  if ieNotDefined('figNum')
    figNum = figure;
  end
  if ieNotDefined('plotMode')
    figureData.currentPlotMode = 'CI Curves';
  else
    figureData.currentPlotMode = plotMode;
  end
  
  % turn off menu/title etc.
  set(figNum,'NumberTitle','off');
  set(figNum,'DefaultAxesDrawMode','fast');
  set(figNum,'DefaultAxesXlimMode','manual');
  set(figNum,'DefaultAxesYlimMode','manual');

  figureData.axesLabels = axesLabels;
  figureData.axesUnits = axesUnits;
  figureData.axesValues = axesValues;
  figureData.axesNames = axesNames;
  figureData.minValue = min(bootstrapMean(:) - lowerErrorBar(:));
  figureData.maxValue = max(bootstrapMean(:) + upperErrorBar(:));
  figureData.actualValue = actualValue;
  figureData.bootstrapMean = bootstrapMean;
  figureData.lowerErrorBar = lowerErrorBar;
  figureData.upperErrorBar = upperErrorBar;
  figureData.statisticalThreshold = statisticalThreshold;
  figureData.handles = struct;
  
   
  figureData.dimensionsOrder = [1 2 3 4]; %order of dimensions we initialize the figure with
  uicontrol('Parent',figNum, 'Unit','normalized', 'Style','Text', 'Position',[.91 .95 .08 .04], 'String','Figure Y Axis');
  figureData.chooseYsubPlotDimHandle = uicontrol('Parent',figNum, 'Unit','normalized', 'Style','popupmenu', 'Position',[.91 .91 .08 .04],...
      'String',axesNames,'Value',4);
  uicontrol('Parent',figNum, 'Unit','normalized', 'Style','Text', 'Position',[.91 .87 .08 .04], 'String','Figure X Axis');
  figureData.chooseXsubPlotDimHandle = uicontrol('Parent',figNum, 'Unit','normalized', 'Style','popupmenu', 'Position',[.91 .83 .08 .04],...
      'String',axesNames(1:3),'Value',3);
  uicontrol('Parent',figNum, 'Unit','normalized', 'Style','Text', 'Position',[.91 .79 .08 .04], 'String','Plots X Axis');
  figureData.chooseXaxisDimHandle = uicontrol('Parent',figNum, 'Unit','normalized', 'Style','popupmenu', 'Position',[.91 .75 .08 .04],...
   'String',axesNames(1:2),'Value',1);
  set(figNum,'userdata',figureData);
 
  set(figureData.chooseYsubPlotDimHandle,'Callback',{@chooseYsubPlotDim,figNum});
  set(figureData.chooseXsubPlotDimHandle,'Callback',{@chooseXsubPlotDim,figNum});
  set(figureData.chooseXaxisDimHandle,'Callback',{@chooseXaxisDim,figNum});
  uicontrol('Parent',figNum, 'Unit','normalized', 'Style','Pushbutton', 'Position',[.91 .7 .08 .04],...
    'String','Redraw Plots','CreateFcn',{@drawPlots,figNum},'Callback',{@drawPlots,figNum});
  
  
  uicontrol('Parent',figNum, 'Unit','normalized', 'Style','Text', 'Position',[.91 .6 .08 .04], 'String','Plot Type');
  uicontrol('Parent',figNum, 'Unit','normalized', 'Style','Text', 'Position',[.91 .35 .08 .04], 'String','Max Y scale');
  uicontrol('Parent',figNum, 'Unit','normalized', 'Style','Text', 'Position',[.91 .2 .08 .04], 'String','Min Y scale');
  uicontrol('Parent',figNum, 'Unit','normalized', 'Style','Pushbutton', 'Position',[.91 .05 .08 .04],...
    'String','Save Data','Callback',{@saveData,figNum});

end

function chooseYsubPlotDim(handle,eventdata,figNum)

  figureData = get(figNum,'userdata');
  ySubPlotValue = get(handle,'value');

  figureData.dimensionsOrder(4) = ySubPlotValue;
  figureData.dimensionsOrder(1:3) = setdiff([1 2 3 4],ySubPlotValue);
  set(figNum,'userdata',figureData);  
  
  ySubPlotDims = get(handle,'string');
  set(figureData.chooseXsubPlotDimHandle,'String',ySubPlotDims(setdiff([1 2 3 4],ySubPlotValue)));
  chooseXsubPlotDim(figureData.chooseXsubPlotDimHandle,[],figNum)

end

function chooseXsubPlotDim(handle,eventdata,figNum)

  figureData = get(figNum,'userdata');
  xSubPlotValue = get(handle,'value');

  remainingDimensions = sort(figureData.dimensionsOrder(1:3));
  figureData.dimensionsOrder(3) = remainingDimensions(xSubPlotValue);
  figureData.dimensionsOrder(1:2) = setdiff(remainingDimensions,remainingDimensions(xSubPlotValue));
  set(figNum,'userdata',figureData);  

  xSubPlotDims = get(handle,'string');
  set(figureData.chooseXaxisDimHandle,'String',xSubPlotDims(setdiff([1 2 3],xSubPlotValue)));
  chooseXaxisDim(figureData.chooseXaxisDimHandle,[],figNum)

end

function chooseXaxisDim(handle,eventdata,figNum)

  figureData = get(figNum,'userdata');
  xAxisValue = get(handle,'value');
  remainingDimensions = sort(figureData.dimensionsOrder(1:2));
  figureData.dimensionsOrder(1) = remainingDimensions(xAxisValue);
  figureData.dimensionsOrder(2) = setdiff(remainingDimensions,remainingDimensions(xAxisValue));
  set(figNum,'userdata',figureData);  

end

function drawPlots(handle,eventdata,figNum)

  set(figNum,'Pointer','watch');drawnow;
  figureData = get(figNum,'userdata');
  plotModes = {'CI Histograms','CI Curves','Estimate Histograms', 'Estimate Curves'};
  plotModeValue = find(ismember(plotModes,figureData.currentPlotMode));

  %delete objects to replace
  handleFields = fieldnames(figureData.handles);
  for iHandles = 1:length(handleFields)
    delete(figureData.handles.(handleFields{iHandles}));
    figureData.handles = rmfield(figureData.handles,handleFields{iHandles}); %this is for when the program crashes during debugging
  end
  set(figNum,'userdata',figureData);                                         %this too

  %Let's reorganize the data according to the order of dimensions
  %Plot order: (X,Z,Xplot,Yplot), 
  %where Yplot and Xplot are the number of subplots in the figure
  % X is the X dimension in each subplot
  % Z is the dimension of the overlapping plots in each subplot
  dimensionsOrder = figureData.dimensionsOrder;
  yPlotLabels = figureData.axesLabels{dimensionsOrder(4)};
  xPlotLabels = figureData.axesLabels{dimensionsOrder(3)};

  zLabels = figureData.axesLabels{dimensionsOrder(2)};
  %make sure we have strings in the cells
  for iLabel = 1:length(zLabels)
    if iscell(zLabels{iLabel})
      zLabels{iLabel} = char(zLabels{iLabel});
    end
  end
  zValues = figureData.axesValues{dimensionsOrder(2)};

  xUnits = figureData.axesUnits{dimensionsOrder(1)};
  xValues = figureData.axesValues{dimensionsOrder(1)};

  actualValue = permute(figureData.actualValue,dimensionsOrder);
  bootstrapMean = permute(figureData.bootstrapMean,dimensionsOrder);
  lowerErrorBar = permute(figureData.lowerErrorBar,dimensionsOrder);
  upperErrorBar = permute(figureData.upperErrorBar,dimensionsOrder);
  intervalSize = upperErrorBar + lowerErrorBar;


  ySubPlots = length(yPlotLabels);
  xSubPlots = length(xPlotLabels);
  zScale = [0 floor(zValues(end)+1)];
  xScale = [0 floor(xValues(end)+1)];
  handles.hPlots = zeros(ySubPlots,xSubPlots);
  hBaselines = zeros(ySubPlots,xSubPlots);

  estimateHistData.hActual = zeros(length(zValues),ySubPlots,xSubPlots);
  estimateHistData.hIntervals = zeros(length(xValues),length(zValues),ySubPlots,xSubPlots);
  intervalHistData.hActual = zeros(length(zValues),ySubPlots,xSubPlots);
  intervalHistData.hIntervals = [];

  estimateCurveData.hActual = zeros(length(zValues),ySubPlots,xSubPlots);
  estimateCurveData.hIntervals = zeros(length(zValues),ySubPlots,xSubPlots);
  intervalCurveData.hActual = zeros(length(zValues),ySubPlots,xSubPlots);
  intervalCurveData.hIntervals = [];

  if ~isempty(figureData.statisticalThreshold)
    statisticalThreshold = permute(figureData.statisticalThreshold,dimensionsOrder);
    %only for estimate histograms
    estimateHistData.hThreshold = zeros(length(xValues),ySubPlots,xSubPlots);
  else
    estimateHistData.hThreshold = [];
    statisticalThreshold = [];
  end
  estimateCurveData.hThreshold = [];
  intervalHistData.hThreshold = [];
  intervalCurveData.hThreshold = [];

  %set the colors
  colors = color2RGB;
  colors = colors([7 5 6 8 4 3 2 1 10]); %remove white and re-order
  for i_color = 1:length(colors)
     colorOrder(i_color,:) = color2RGB(colors{i_color});
  end
  if length(zValues)>size(colorOrder,1)
     colorOrder = repmat(colorOrder,ceil(length(zValues)/size(colorOrder,1)),1);
  end
  colorOrder = colorOrder(1:length(zValues),:);
  %set color order for plots
  set(figNum,'DefaultAxesColorOrder',colorOrder);
  %for bars, need to set the colormap
  set(figNum,'colormap',colorOrder);

  %plot legends
  if length(zLabels)>1
    legendWidth = .08;
    %Histograms legend
    subplot('position',[0 .95-2*length(zLabels)*.03 legendWidth length(zLabels)*.03]);
    h=bar([1 2],zeros(2,length(zLabels)),'grouped','edgecolor','none');
    for i=1:length(h)
      set(get(h(i),'baseline'),'visible','off');
    end
    axis off
    handles.hHistLegendBox = legend(zLabels,'position',[0 .95-2*length(zLabels)*.03 legendWidth length(zLabels)*.03]);
    set(handles.hHistLegendBox,'Interpreter','none','box','off');
    figureData.hHistogramsLegend = get(handles.hHistLegendBox,'children');
    set(figureData.hHistogramsLegend,'visible','off');

    %Curves legend
    subplot('position',[0 .95-length(zLabels)*.03 legendWidth length(zLabels)*.03]);
    plot(zeros(length(zLabels)),zeros(length(zLabels)));
    axis off
    handles.hCurvesLegendBox = legend(zLabels,'position',[0 .95-length(zLabels)*.03 legendWidth length(zLabels)*.03]);
    set(handles.hCurvesLegendBox,'Interpreter','none','box','off');
    figureData.hCurvesLegend = get(handles.hCurvesLegendBox,'children');
    set(figureData.hCurvesLegend,'visible','off');
  else
    legendWidth = 0;
    figureData.hHistogramsLegend = [];
    figureData.hCurvesLegend = [];
  end

  %subplot dimensions
  outerMargins = [legendWidth+.05 .1 .1 .03]; %(left bottom right top)
  innerMargins = [.01 .01]; %(horizontal vertical)
  subPlotDims = [1-sum(outerMargins([1 3])) 1-sum(outerMargins([2 4]))];%(horizontal vertical)
  subPlotPosition(3) = (subPlotDims(1) - (xSubPlots-1)*innerMargins(1))/xSubPlots;
  subPlotPosition(4) = (subPlotDims(2) - (ySubPlots-1)*innerMargins(2))/ySubPlots;

  for ySubPlot = 1:ySubPlots
    for xSubPlot = 1:xSubPlots
      subPlotPosition(1) = outerMargins(1)+ (xSubPlot-1)*(innerMargins(1)+ subPlotPosition(3));
      subPlotPosition(2) = outerMargins(2)+ (ySubPlots-ySubPlot)*(innerMargins(2)+ subPlotPosition(4));
      %handles.hPlots(ySubPlot,xSubPlot) = subplot(ySubPlots, xSubPlots,(ySubPlot-1)*xSubPlots+xSubPlot);
      handles.hPlots(ySubPlot,xSubPlot) = subplot('position',subPlotPosition);
      axis([zScale figureData.minValue figureData.maxValue])
      hold on;

      %plot estimate histograms
      if size(actualValue,1)==1
        set(gca,'nextPlot','add');
        normalisedZvalues = xValues + (zValues-mean(zValues))/(zValues(end) - zValues(1));
        for z = 1:size(actualValue,2)
          estimateHistData.hActual(z,ySubPlot,xSubPlot) = bar(...
            normalisedZvalues(z),...
            actualValue(1,z,xSubPlot,ySubPlot)',...
            'faceColor',colorOrder(z,:),'edgecolor','none',...
            'barwidth',.8/length(zValues),'visible','off')';
        end
      else
        estimateHistData.hActual(:,ySubPlot,xSubPlot) = bar(...
          xValues,actualValue(:,:,xSubPlot,ySubPlot),...
          'grouped','edgecolor','none','visible','off')';
      end
      hBaselinesToDelete = get(estimateHistData.hActual(:,ySubPlot,xSubPlot),'baseline');
      if iscell(hBaselinesToDelete)
        hBaselinesToDelete = cell2mat(hBaselinesToDelete);
      end
      delete(hBaselinesToDelete);

      %plot interval size histograms
      if size(actualValue,1)==1
        set(gca,'nextPlot','add');
        for z = 1:size(actualValue,2)
          intervalHistData.hActual(z,ySubPlot,xSubPlot) = bar(...
            normalisedZvalues(z),...
            intervalSize(1,z,xSubPlot,ySubPlot)',...
            'faceColor',colorOrder(z,:),'edgecolor','none',...
            'barwidth',.8/length(zValues),'visible','off')';
        end
      else
        intervalHistData.hActual(:,ySubPlot,xSubPlot) = bar(...
          xValues,intervalSize(:,:,xSubPlot,ySubPlot),...
          'grouped','edgecolor','none','visible','off')';
      end
      hBaselinesToDelete = get(intervalHistData.hActual(:,ySubPlot,xSubPlot),'baseline');
      if iscell(hBaselinesToDelete)
        hBaselinesToDelete = cell2mat(hBaselinesToDelete);
      end
      delete(hBaselinesToDelete);
      %plot my own baseline
      hBaselines(ySubPlot,xSubPlot) = plot(xScale, [0 0],'--k','visible','off'); 

      if length(xValues)==1
        %plot estimate curves
        estimateCurveData.hActual(:,ySubPlot,xSubPlot) = plot(...
          normalisedZvalues', actualValue(:,:,xSubPlot,ySubPlot)',...
          'k','visible','off');
        %plot interval size curves
        intervalCurveData.hActual(:,ySubPlot,xSubPlot) = plot(...
          normalisedZvalues', intervalSize(:,:,xSubPlot,ySubPlot)',...
          'k','visible','off');
      else
        %handle NaNs
        curvesXvalues = repmat((xValues)',1,size(actualValue,2));
        curvesXvalues(isnan(actualValue(:,:,xSubPlot,ySubPlot)))=nan;
        %plot estimate curves
        estimateCurveData.hActual(:,ySubPlot,xSubPlot) = plot(...
          curvesXvalues,...
          actualValue(:,:,xSubPlot,ySubPlot),...
          'visible','off');
        %plot interval size curves
        intervalCurveData.hActual(:,ySubPlot,xSubPlot) = plot(...
          curvesXvalues,...
          intervalSize(:,:,xSubPlot,ySubPlot),...
          'visible','off');
      end

      %plot estimate histogram error bars (=confidence intervals)
      for z = 1:size(actualValue,2)
        set(gca,'nextPlot','add');
        barPositions = get(get(estimateHistData.hActual(z,ySubPlot,xSubPlot),'Children'),'XData');
        estimateHistData.hIntervals(z,ySubPlot,xSubPlot) = errorbar(...
          mean(barPositions),...
          bootstrapMean(:,z,xSubPlot,ySubPlot),...
          lowerErrorBar(:,z,xSubPlot,ySubPlot),...
          upperErrorBar(:,z,xSubPlot,ySubPlot),...
          'xk','visible','off');

        if ~isempty(statisticalThreshold)
          estimateHistData.hThreshold(:,z,ySubPlot,xSubPlot) = plot(...
            barPositions(2:3,:),...
            [statisticalThreshold(:,z,xSubPlot,ySubPlot) statisticalThreshold(:,z,xSubPlot,ySubPlot)]',...
            '--k','visible','off');
        end
      end

      %plot estimate curve error bars (=confidence intervals)
      if length(xValues)==1
        set(gca,'nextPlot','add');
        for z = 1:length(normalisedZvalues);
        estimateCurveData.hIntervals(z,ySubPlot,xSubPlot) = errorbar(...
          normalisedZvalues(z),...
          bootstrapMean(:,z,xSubPlot,ySubPlot),...
          lowerErrorBar(:,z,xSubPlot,ySubPlot),...
          upperErrorBar(:,z,xSubPlot,ySubPlot),...
          'x','color',colorOrder(z,:),'visible','off')';
        end
      else
        estimateCurveData.hIntervals(:,ySubPlot,xSubPlot) = errorbar(...
          curvesXvalues,...
          bootstrapMean(:,:,xSubPlot,ySubPlot),...
          lowerErrorBar(:,:,xSubPlot,ySubPlot),...
          upperErrorBar(:,:,xSubPlot,ySubPlot),...
          'x','visible','off')';
      end

      %set the titles, scale, units and ticks
      set(handles.hPlots(ySubPlot,xSubPlot),'Xlim',xScale);
      set(gca,'xTick',xValues);
      if ySubPlot == 1
        title(xPlotLabels{xSubPlot});
      end
      if xSubPlot == 1
        ylabel(yPlotLabels{ySubPlot},'interpreter','none');
      else
        set(gca,'yTickLabel',[]);
      end
      if ySubPlot == size(handles.hPlots,1)
        xlabel(xUnits);
      else
        set(gca,'xTickLabel',[]);
      end

    end
  end


  if ~isempty(statisticalThreshold)
    handles.showThreshold=uicontrol('Parent',figNum, 'Unit','normalized', 'Style','checkBox', 'Position',[.91 .45 .08 .04],...
    'String','Show Statistical Threshold','Value',1);
  else
    handles.showThreshold=[];
  end

  handles.showBaseline = uicontrol('Parent',figNum, 'Unit','normalized', 'Style','checkBox', 'Position',[.91 .4 .08 .04],...
        'String','Show Baseline','Value',1,'CreateFcn',{@showObject,hBaselines},'Callback',{@showObject,hBaselines});
  handles.maxYEdit = uicontrol('Parent',figNum, 'Unit','normalized', 'Style','Edit', 'Position',[.91 .3 .08 .04],...
       'CreateFcn',{@setYScale,handles.hPlots,'max'},'Callback',{@setYScale,handles.hPlots,'max'}, 'String',num2str(figureData.maxValue));
  handles.minYEdit = uicontrol('Parent',figNum, 'Unit','normalized', 'Style','Edit', 'Position',[.91 .25 .08 .04],...
     'CreateFcn',{@setYScale,handles.hPlots,'min'},'Callback',{@setYScale,handles.hPlots,'min'}, 'String',num2str(figureData.minValue));
  figureData.handles = handles;
  figureData.estimateHistData = estimateHistData;
  figureData.estimateCurveData = estimateCurveData;
  figureData.intervalHistData = intervalHistData;
  figureData.intervalCurveData = intervalCurveData;

  set(figNum,'userdata',figureData);

  figureData.handles.hPlotMode = uicontrol('Parent',figNum, 'Unit','normalized', 'Style','popupmenu', 'Position',[.91 .55 .08 .04],...
        'String',plotModes,'Value',plotModeValue,...
        'CreateFcn',{@switchHistogramCurves,figNum},...
        'Callback',{@switchHistogramCurves,figNum});
  set(figNum,'userdata',figureData);
  set(figNum,'Pointer','arrow');drawnow;

end

function setYScale(handle,eventdata,hPlots, minmax)

  scale = get(hPlots(1,1),'Ylim');
  new_value = str2num(get(handle,'string'));
  if strcmp(minmax,'max')
     scale(2) = new_value;
  else
     scale(1) = new_value;
  end
  for i_plot = 1:numel(hPlots)
     set(hPlots(i_plot),'Ylim',scale);
  end

end

function showObject(handle,eventdata,hObject)

  if get(handle,'value')
    set(hObject,'visible','on');
  else
    set(hObject,'visible','off');
  end

end

function switchHistogramCurves(handle,eventdata,figNum)

figureData = get(figNum,'userdata');

  switch figureData.currentPlotMode
    case 'Estimate Histograms'
      invisibleObjects = figureData.estimateHistData;
      invisibleLegend = figureData.hHistogramsLegend;
    case 'Estimate Curves'
      invisibleObjects = figureData.estimateCurveData;
      invisibleLegend = figureData.hCurvesLegend;
    case 'CI Histograms'
      invisibleObjects = figureData.intervalHistData;
      invisibleLegend = figureData.hHistogramsLegend;
    case 'CI Curves'
      invisibleObjects = figureData.intervalCurveData;
      invisibleLegend = figureData.hCurvesLegend;
  end

  plotModes = get(handle,'string');
  plotModeValue = get(handle,'value');
  switch plotModes{plotModeValue}
    case 'Estimate Histograms'
      visibleObjects = figureData.estimateHistData;
      visibleLegend = figureData.hHistogramsLegend;
    case 'Estimate Curves'
      visibleObjects = figureData.estimateCurveData;
      visibleLegend = figureData.hCurvesLegend;
    case 'CI Histograms'
      visibleObjects = figureData.intervalHistData;
      visibleLegend = figureData.hHistogramsLegend;
    case 'CI Curves'
      visibleObjects = figureData.intervalCurveData;
      visibleLegend = figureData.hCurvesLegend;
  end
  figureData.currentPlotMode = plotModes{plotModeValue};

  set(invisibleObjects.hActual,'visible','off');
  set(invisibleObjects.hIntervals,'visible','off');
  set(invisibleLegend,'visible','off');
  set(invisibleObjects.hThreshold,'visible','off');
  set(figureData.handles.showThreshold,'CreateFcn',{@showObject,visibleObjects.hThreshold},'Callback',{@showObject,visibleObjects.hThreshold});
  set(visibleObjects.hActual,'visible','on');
  set(visibleObjects.hIntervals,'visible','on');
  set(visibleLegend,'visible','on');

  set(figNum,'userdata',figureData);

end

function saveData(handle,eventdata,figNum)

  [filename, pathname] = uiputfile([pwd '/*.mat'],'Mat file name');

  if ~isnumeric(filename)
    figureData = get(figNum,'userdata');
    save([pathname filename],'-struct','figureData',...
      'actualValue','upperErrorBar','lowerErrorBar','bootstrapMean',...
      'axesNames','axesLabels','axesUnits','axesValues');
    if isfield(figureData,'statisticalThreshold') && ~isempty(figureData.statisticalThreshold)
      save([pathname filename],'-struct','figureData','statisticalThreshold','-append');
    end
  end

end
