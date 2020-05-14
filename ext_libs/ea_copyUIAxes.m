function h = copyUIAxes(varargin)
%  The use of Matlab's copyobj() is not supported with UI axes. 
% use COPYUIAXES to copy the content of a UI axis to a new figure.
% COPYUIAXES receives the handle to a UI axes and copies all of its
% children and most of its properties to a new axis. If the UI axis
% has a legend, the legend is copied too (but see tech notes below). 
% The handle to a colobar may also be included to copy the colorbar.
% Requires Matlab r2016a or later. 
%
% COPYUIAXES(uiax) creates a new figure and axis in default position
% and copies the content of the UI axis onto the new axis. 'uiax' is
% required to be a UIAxis handle of class matlab.ui.control.UIAxes.
% If the axis handle contains a non-empty legend property, the legend
% is copied too (but see tech notes below).
%
% COPYUIAXES(uiax,destination) copies the content of the UI axis to 
% destination object which can either be a figure or axes handle. If
% the handle is a figure the axes will be created in default position.
%
% h = COPYUIAXES(___) returns a scalar structure containing all graphics
% handles that were produced.  
%
% COPYUIAXES(___, 'legend', h) specifies the legend handle (h) to 
% copy to the new destination. If the lengend handle is provided 
% here, it overrides any legend detected from within the axis handle.
%
% COPYUIAXES(___, 'colorbar', h) specifies the colorbar handle (h)  
% to copy to the new destination (but see tech notes below). 
%
% COPYUIAXES(___, 'listIgnoredProps', TF) when TF is true, a table
% of graphic objects appears in the command window listing the properties
% of each object that were ignored. Some properties are not editable
% while others are intentionally ignored.  
%
% Examples (for r2016b or later)
%     fh = uifigure();
%     uiax = uiaxes(fh);
%     hold(uiax,'on')
%     grid(uiax,'on')
%     x = linspace(0,3*pi,200);
%     y = cos(x) + rand(1,200);
%     ph = scatter(uiax,x,y,25,linspace(1,10,200),'filled');
%     lh = plot(uiax, x, smooth(x,y), 'k-'); 
%     title(uiax,'copyUIaxes Demo','FontSize',18)
%     xlabel(uiax,'x axis')
%     ylabel(uiax,'y axis')
%     zlabel(uiax,'z axis')
%     cb = colorbar(uiax); 
%     cb.Ticks = -2:2:12;
%     caxis(uiax,[-2,12])
%     ylabel(cb,'cb ylabel')
%     lh = legend(uiax,[ph,lh],{'Raw','lowess'},'Location','NorthEast'); 
%     drawnow()% all graphics to catch up before copying
% 
%     % Ex 1: Copy axis to another figure
%     h = copyUIAxes(uiax);
%     
%     % Ex 2: specify ax axis as destination
%     fNew = figure();
%     axh = subplot(2,2,1,'Parent',fNew);
%     h = copyUIAxes(uiax,axh);
% 
%     % Ex 3: Provide colorbar and legend handles
%     h = copyUIAxes(uiax, 'colorbar', cb, 'legend', lh);
%     
%     % Ex 4: See which properties were not copied
%     h = copyUIAxes(uiax, 'colorbar', cb, 'listIgnoredProps', true);
%    
% To report bugs, feature requests, or high fives, the author's 
% email address is within the file. 

% Adam Danz Oct 2019.
% Copyright (c) 2019, Adam Danz
% All rights reserved

% Please contact me if you experience any bugs. 
% Email: (Run the line below)
% char([97 100 97 109 46 100 97 110 122 64 103 109 97 105 108 46 99 111 109]) 

% To follow discussion on matlab central see note [1].

% Version history
% 191000  vs 1.0.0 first upload to FEX
% 200209  vs 1.1.0 return fig handle; improved error handling; improved listIgnoredProps text.
% 200319  vs 1.2.0 Now works with categoricalHistograms; if set() results in warnings, a final 
%       warning msg appears; added some commented-out troubleshooting code; new newFieldOrder();
%       Added class to addRequired(uiax) error msg; drawnow() updates figs before copy; axes
%       made visible at the end only if UIAxes are visible;

%% Tech notes
% Mimicking copyobj() for UIaxes is not trivial.  Many of the axes properties
% and sub-handles are hidden or obscured so some properties are not copied. 
% Copying the legend is not currently supported so we have to create a new one.
% This is also not trivial since we do not have a 1:1 match between line 
% object handles and their DisplayName values from within the legend handle. 
% So we must re-create the legend and search for the DisplayName values in 
% the axis children.  This may result in inconsistencies with the original 
% legend.  For example, if there is more than 1 object with the same DisplayName
% value then the legend will pair with the first one which may not be the same 
% object in the original legend.  In these cases it is better to delete the new
% legend and just recreate it.  Similar problems exist with copying the 
% colorbar.  The colorbar cannot be copied so we must recreated it and do 
% our best to copy the properties. 
%
% I've chosen not to copy the axis position property.  Typically UIaxes are 
% smaller because they are embedded in an app.  If you want the axes to 
% match in size you can make those adjustments after the recreation.  
%
% See footnotes for additional details regarding certain lines of the code. 

%% input parser
p = inputParser();
p.FunctionName = mfilename;
addRequired(p, 'uiax', @(h)assert(isa(h,'matlab.ui.control.UIAxes'),...
    sprintf(['copyUIAxes() only copies from UIAxes. You''re attemping to copy a %s ',...
    'object. Use copyobj() to copy objects from other axis types.'], class(h))))
addOptional(p, 'destination', [], @(h)assert(isgraphics(h,'axes')||isgraphics(h,'figure'),... % see notes [9,10]
    'Destination must be an axis handle or figure handle.'));  
addParameter(p, 'legend', [], @(h)isgraphics(h,'legend'));
addParameter(p, 'colorbar', [], @(h)isgraphics(h,'colorbar')); 
addParameter(p, 'listIgnoredProps', false, @islogical)
parse(p,varargin{:})

%% Produce figure and axes if needed
if isempty(p.Results.destination)
    h.figure = figure('Visible','off'); 
    h.axes = axes(h.figure); 
    figureCreatedInternal = true; 
elseif isgraphics(p.Results.destination, 'figure')
    h.figure = p.Results.destination; 
    h.axes = axes(h.figure); 
    figureCreatedInternal = true; 
else % we assume it's an axis handle; see notes [9,10]
    h.figure = p.Results.destination.Parent; 
    h.axes = p.Results.destination; 
    figureCreatedInternal = false; 
end

drawnow() % To avoid lagging graphics problems

%% Adjust axes with CategoricalRulers [16]
% As far as I'm aware, this only affects histogram() plots with 
% categorical data along the x or y axis. I've tested pie() and
% scatter plots with categorical data and there are no problems. 
% histogram2() and other hist plots do not allow categorical inputs (as of r2019b). 
% This section must come before assigning data to the axis. 

% Get handle(s) to categorical histogram in UIAxes.
catHistHandle = findall(p.Results.uiax.Children, 'Type', 'categoricalhistogram'); 

% Detect if X and Y axis are CategoricalRuler
if ~isempty(catHistHandle) && strcmpi(class(p.Results.uiax.XAxis), 'matlab.graphics.axis.decorator.CategoricalRuler')
    matlab.graphics.internal.configureAxes(h.axes, catHistHandle(1).Data, 1); % set xaxis to categorical [16]
end
if ~isempty(catHistHandle) && strcmpi(class(p.Results.uiax.YAxis), 'matlab.graphics.axis.decorator.CategoricalRuler')
    matlab.graphics.internal.configureAxes(h.axes, 1, catHistHandle(1).Data); % set yaxis to categorical [16]
end

%% Copy axis children and (most) properties
% Copy all children from UIAxes to new axis
copyobj(p.Results.uiax.Children, h.axes)

% Anonymous func used to move selected fields to the end of a structure 
% INPUTS: 
%   S: input structure
%   propList: nx1 or 1xn cell array of chars listing fields of S to put at the end. If field is missing, ignored.
% OUTPUT: The structure S with the reordered fields. 
% EXAMPLE; Snew = newFieldOrder(uiaxGoodParams, {'XLim','YLim','ZLim'})
newFieldOrder = @(S, propList)orderfields(S,[find(~ismember(fields(S),propList)); find(ismember(fields(S),propList))]); 

% Choose selected properties to copy
copyLaterList = {'Title';'XLabel';'YLabel';'ZLabel'}; % see note [13]
excludeList = [{'Parent';'Children';'XAxis';'YAxis';'ZAxis';'Position';'OuterPosition';'Units';'Toolbar'};copyLaterList]; 
[uiaxGoodParams, uiaxbadProps] = getGoodParams(p.Results.uiax, h.axes, excludeList, 'axis');  % see note [2]
uiaxbadProps(ismember(uiaxbadProps.axis,copyLaterList),:) = []; %rm the fields that we'll copy later.

% move axis limit fields last or they may not be set properly
uiaxGoodParams = newFieldOrder(uiaxGoodParams, {'XLim','YLim','ZLim'}); 

% clear the last warning in order to detect new warnings.
[lastWarnMsg, lastWarnID] = lastwarn(); 
lastwarn('') % Clear last warning so we can detect if one appears.  

% set properties 
set(h.axes, uiaxGoodParams)

% For trouble shooting, or to loop through each property rather then setting them all at once,
% replace the line of the code above with this for-loop below which indicates the property
% being edited in the command window.
% To remove fields for testing purposes, uiaxGoodParams = rmfield(uiaxGoodParams, {'XTick', 'XLim'}); 
%     allf = fieldnames(uiaxGoodParams);
%     for i = 1:length(allf)
%         fprintf('Property #%d: %s\n', i, allf{i});
%         if i == 48
%             pause(0.1)
%         end
%         % get(p.Results.uiax, allf{i})    % show current value for UIAxis
%         % get(h.axes, allf{i})            % show current value for normal axis
%         set(h.axes, allf{i}, uiaxGoodParams.(allf{i}))
%     end

% Copy title and axis labels; see note [13]
h.axesTitle = title(h.axes, p.Results.uiax.Title.String);
[axTtlGoodParams, axTtlbadProps] = getGoodParams(p.Results.uiax.Title, h.axesTitle, {'Parent'; 'Position'}, 'axesTitle');
set(h.axesTitle, axTtlGoodParams)

h.axesXLabel = xlabel(h.axes, p.Results.uiax.XLabel.String);
[axXLGoodParams, axXLbadProps] = getGoodParams(p.Results.uiax.XLabel, h.axesXLabel, {'Parent'; 'Position'}, 'axesXLabel');
set(h.axesXLabel, axXLGoodParams)

h.axesYLabel = ylabel(h.axes, p.Results.uiax.YLabel.String);
[axYLGoodParams, axYLbadProps] = getGoodParams(p.Results.uiax.YLabel, h.axesYLabel, {'Parent'; 'Position'}, 'axesYLabel');
set(h.axesYLabel, axYLGoodParams)

h.axesZLabel = zlabel(h.axes, p.Results.uiax.ZLabel.String);
[axZLGoodParams, axZLbadProps] = getGoodParams(p.Results.uiax.ZLabel, h.axesZLabel, {'Parent'; 'Position'}, 'axesZLabel');
set(h.axesZLabel, axZLGoodParams)

%% Detect legend and copy if one exists (see note [14])
if (any(strcmpi(properties(p.Results.uiax),'Legend')) && ~isempty(p.Results.uiax.Legend)) ... % see notes [6,8]
        || ~isempty(p.Results.legend) 
    % if Legend was provided, use that handle, otherwise use the detected one. 
    if ~isempty(p.Results.legend)
        legHand = p.Results.legend;
    else
        legHand = p.Results.uiax.Legend; 
    end
    % Search for objects in new axes that have matching displayNames values as legend strings (see note [11])
    newChildren = h.axes.Children; 
    dispNames = get(newChildren,'DisplayName'); 
    [~,legIdx] = ismember(legHand.String, dispNames); 
    legObjHands = newChildren(legIdx); 
    
    % Create new legend and copy selected properties
    h.legend = legend(h.axes,legObjHands,legHand.String); % see note [7]
    excludeList = {'String';'Parent'; 'Children'; 'Position'; 'Units'; 'UIContextMenu'}; 
    [legGoodParams, legbadProps] = getGoodParams(legHand, h.legend, excludeList, 'legend');  % see note [3]
    set(h.legend, legGoodParams)
else
    legbadProps = table(); 
end

%% Detect colorbar and copy if one exists (see note [14])
% Note, as of r2019b, there is no way I know of to detect colorbar or get its handle from ui axes (email me if you know how).
if  ~isempty(p.Results.colorbar) 
    %Copy colorbar & selected properties
    h.colorbar = colorbar(h.axes); % see note [3]
    excludeList = {'Parent'; 'Children'; 'Position'; 'Units'; 'UIContextMenu'};
    [cbGoodParams, cbbadProps] = getGoodParams(p.Results.colorbar, h.colorbar, excludeList, 'colorbar'); 
    set(h.colorbar, cbGoodParams)
       
    % Copy title
    h.colorbarTitle = title(h.colorbar, p.Results.colorbar.Title.String);
    [cbTtlGoodParams, cbTtlbadProps] = getGoodParams(p.Results.colorbar.Title, h.colorbarTitle, {'Parent'; 'Children';'Position'}, 'colorbarTitle'); 
    set(h.colorbarTitle, cbTtlGoodParams)
    
    % Copy ylabels
    h.colorbarYlabel = ylabel(h.colorbar, p.Results.colorbar.YLabel.String);
    [cbYLabGoodParams, cbYLabbadProps] = getGoodParams(p.Results.colorbar.YLabel, h.colorbarYlabel, {'Parent'; 'Children';'Position'},'colorbarYlabel'); 
    set(h.colorbarYlabel, cbYLabGoodParams)
else
    cbbadProps = table(); 
    cbTtlbadProps = table();
    cbYLabbadProps = table(); 
end

%% If we just created the figure|axes, turn on its visibility
% Exception: when the UIAxis visibility is off.
if figureCreatedInternal && strcmpi(p.Results.uiax.Visible, 'on')
    h.figure.Visible = 'on'; 
    h.axes.Visible = 'on'; 
end

%% Show summary of properties that were not copied
if p.Results.listIgnoredProps
    % Put all badProps arrays into table
    badProps = {uiaxbadProps, axTtlbadProps, axXLbadProps, axYLbadProps, axZLbadProps, ...
        legbadProps, cbbadProps, cbTtlbadProps, cbYLabbadProps};
    numProps = cellfun(@numel,badProps);
    badProps(numProps==0) = []; 
    badParams = cellfun(@(c)[c;repmat({' '},max(numProps)-numel(c),1)],badProps,'UniformOutput',false);
    badParams = [badParams{:}];
    % Display table
    if ~isempty(badParams)
        fprintf(['\n%s\nThe following properties (rows) were not copied from the following objects (columns)\n'...
            'either because they are not editable or because they were intentionally ignored.\n'...
            'See %s for details.\n'], char(42*ones(1,70)), sprintf('<a href="matlab: help(''%s'') ">%s.m</a>', which(mfilename), mfilename))
        disp(badParams)
        fprintf('%s\n\n',char(42*ones(1,70)))
    else
        fprintf('\n %s\n All properties were copied.\n See %s for details.\n%s\n\n ',char(42*ones(1,70)), ...
            sprintf('<a href="matlab: help(''%s'') ">%s.m</a>', which(mfilename), mfilename), char(42*ones(1,70)))
    end
end

%% Check if a warning was thrown
pause(0.05) % set() needs time to complete prior to determining whether it caused an error.
if ~isempty(lastwarn())
    % A warning was thrown; tell user to inspect copied figure    
    msg = [newline, 'Based on the warning messages shown above, some properties may not have '...
        'been copied properly.  This can happen when a property change triggers changes '...
        'to other properties or due to the order in which the properties were set. Inspect ' ...
        'the copied axis closely and set neglected properties after running copyUIAxes().', newline];
    defaultMode = warning('query', 'backtrace');  %store default
    warning off backtrace                         %turn off backtrace
    fprintf('\n')
    warning(msg)
    warning(defaultMode.state, 'backtrace')       %turn back on default
    
else
    % A warning was not thrown; return the lastwarning.
    lastwarn(lastWarnMsg, lastWarnID)
   
end

%% Local functions
function [goodParams, badParams] = getGoodParams(hObjSource, hObjDest, badList, objName)
% The goal is to create a structure of parameter values to be changed in the destination object.  
% INPUTS
% hObjSource: handle to object being copied (ie, the old legend)
% hObjDest: handle to the new, existing copied object (ie, the new legend)
% badList: a nx1 cell array of strings that match property names (case sensitive).  These
%   properties will not be copied and will be removed from the structure.  
% objName: a char array identifying the object in hObhDest (for badParams output headers)
% OUTPUTS
% goodParams: a strcture of parameter names (fields) and their values that will be used
%   to update the parameters in hObjDest. To apply the parameter values, set(hObjDest, goodParams).
% badParams: a mx1 table of strings listing the properties that will not be copied either because
%   they were listed in badList or because they are not editable. Headers are defined by objName.

% List all parameters of the source object
params = get(hObjSource);
paramNames = fieldnames(hObjSource);

% Get list of editable params in destination object
editableParams = fieldnames(set(hObjDest));

% Remove the params that aren't editable or are unwanted in the destination obj
badParams = paramNames(~ismember(paramNames, editableParams));
badParams = unique([badParams; badList(ismember(badList,paramNames))]); % see note [4]
goodParams = rmfield(params,unique(badParams)); % see note [5]
badParams = table(sort(badParams),'VariableNames', {objName}); % see note [12]

% Change the order of properties (for troubleshooting purposes only)
% Properties in col1, if they exist in goodParams, should be set
% prior to their paired property in col2.  
%     AbeforeB = {
%         'ColorOrderIndex',           'ColorOrder'           % [15]
%         'LineStyleOrderIndex',       'LineStyleOrder'       % [15]
%         }; 
% 
%     params = fields(goodParams);
%     for i = 1:size(AbeforeB,1)
%         [~, paramIdx, ABidx] = intersect(params, AbeforeB(i,:));
%         if ABidx(1) > ABidx(2)
%             % Switch positions
%             paramOrder = (1:numel(params)).';
%             paramOrder(paramIdx) = flipud(paramIdx);
%             goodParams = orderfields(goodParams,paramOrder);
%         end
%     end


%% Notes
% [1] >580 views in past 30-days as of Oct 2019; >660 as of Feb 2019.
%   https://www.mathworks.com/matlabcentral/answers/281318-how-can-i-save-a-figure-within-app-designer
% [2] Some axes properties are read-only so we can't copy those over to the new axes.  Also, there are 
%   some axes properties that you don't want to copy over such as parent, children, X/Y/ZAxis, which will 
%   cause errors.  In addition to those, we will not copy the "Position" and "OuterPosition" properties
%   since the new axis position is pre-defined. 
% [3] Similar to [2], some legend and colorbar properties should not be copied.  I've chosen to not copy 
%   position due to the potential size differences between the new and old axes and we don't know what units 
%   the legend/cb contains so it could end up way off the plot.  User can adjust position afterwards.
% [4] ZAxis (and maybe other properties) didn't exist in earlier releases (r2016a for example).
% [5] unique() needed because in earlier maltab release some of the properties in my additional bad field 
%   list are already included and duplicates cause error in rmfield().
% [6] Legend is currently a property of the UI axis. Previously it was not.  
%   Legends were not supported in app designer until r2016b:
%   https://www.mathworks.com/help/releases/R2016b/matlab/creating_guis/graphics-support-in-app-designer.html
% [7] In 2019a you can just call legend(h) with no strings and the legend obj will be created.  Testing in 
%   r2016b reveals that this method just returns an empty graphics obeject and a string is required.  This
%   is why we're copying the original legend strings over as soon as we create the legend. 
% [8] Unfortunately we can't use isprop(p.Results.uiax,'Legend') because it returns true in r2016a but then
%   when you try to access that property you get error: You cannot get the 'Legend' property of UIAxes. 
% [9] If I ever allow non-UIAxes, this test would have be a lot more flexible. The axis might be a polaraxis
%   or geoaxis etc but 'axes' will only confirm cartesian axes.  Here's a sloppy alternative:
%   @(h)~isempty(regexp(class(h),'^matlab.graphics.axis.','once'))
% [10] List of supported axes and plots for app designer: 
% 	https://www.mathworks.com/help/matlab/creating_guis/graphics-support-in-app-designer.html
% [11] Legends must be copied with the uiaxes in copyobj() but since we can't use copyobj()
%   with UI axes, we cannot copy the legend. We created a new one but merely copying the 
%   legend properties does not preserve the 1:1 relationship between object handles and
%   thier legend strings.  For example, if the source plot has 2 objects but only the 2nd
%   one is represented in the legend, when the new legend is created it will detect that 
%   only 1 object is represented but it will show the displayname for the 1st object. 
%   Getting around this is really messy and the solution here may not always work.  If 
%   that's the case, it may be better to rebuild the legend. 
% [12] Several of the remaining fields contain sub-handles to graphics objects.  I'm not sure if copying
%   their values will results in any problems. To see which fields contain handles, 
%   goodParamIsHand = structfun(@ishghandle,goodParams,'UniformOutput',false);
% [13] The title, x/y/zLabel properties would copy fine from the source-axis but they have a 'parent' property 
%   that, when changed, would remove those objects from the source-axis rather than copying them.  There are also
%   some read-only properties of these objects that can't be copied.  Copying the position property can results in 
%   the label/title going off the figure border.  
% [14] Below are links that show the rollout states of UIAxis support
%   * 2016a - 2017a : https://www.mathworks.com/help/releases/R2017a/matlab/creating_guis/graphics-support-in-app-designer.html
%   * 2017b - 2018a : https://www.mathworks.com/help/releases/R2018a/matlab/creating_guis/graphics-support-in-app-designer.html
%   * 2018b         : https://www.mathworks.com/help/releases/R2018b/matlab/creating_guis/graphics-support-in-app-designer.html
%   * 2019a         : https://www.mathworks.com/help/releases/R2019a/matlab/creating_guis/graphics-support-in-app-designer.html
%   * current       : https://www.mathworks.com/help/matlab/creating_guis/graphics-support-in-app-designer.html
% [15] See the link below regarding ColorOrderIndex & LineStyleOrderIndex.
%   https://www.mathworks.com/help/matlab/ref/matlab.graphics.axis.axes-properties.html#mw_23890c9d-6e04-4aab-92b5-f873dc229766
% [16] Prior to vs 1.2 user reported an error copying a categorical histogram (link A below). Categorical histograms use
%   a CategoricalRuler along the x or y axes [link B] while the default is a NumericRuler when an axis is created.  Categorical
%   data can be copied to a numericRuler axes but an error occurs when you copy the axis properties from a CategoricalRuler to
%   a NumericRuler.  Specifically, the XTick and XLim (or y-) values are problematic. Other categorical data plot (pie, scatter,
%   see link C) do not have this problem, only histogram() to my knowledge (as of vs 1.2.0).  To fix this, NumericRuler must be
%   converted to CategoricalRuler **prior to** copying the data to the axes.  The line that converts the axes was pulled from 
%   histogram (line 150) > categorical/categoricalHistogram (line 108) in r2019b, update 2. 
%   [A] https://www.mathworks.com/matlabcentral/answers/510081-how-to-export-xaxistick-labels-which-are-cell-arrays-to-figure-or-powerpoint-in-app-designer
%   [B] https://www.mathworks.com/help/matlab/ref/matlab.graphics.axis.decorator.categoricalruler-properties.html
%   [C] https://www.mathworks.com/help/matlab/matlab_prog/plot-categorical-data.html
% [17] Copying the XAxis and ZAxis come with all sorts of problems.  Many of the XAxis and YAxis properties are copied
%   when we copy the main axis handles properties, anyway.  For some reason, the YAxis is read-only and can't be copied
%   (see https://www.mathworks.com/matlabcentral/answers/438687-how-do-i-set-yaxis-of-axes-from-a-numeric-ruler-to-a-datetime-ruler).
%   The code below attempts to copy the X/ZAxes but I'm leaving it out for now. 
%     Try to copy the XAxis & ZAxis
%     YAxis is read-only for some reason [17].  Sometimes copying the X/ZAxis works but not always.
%     xzcopied = [false, false]; 
%     xzcopyProp = {'XAxis','ZAxis'};  % Must be exact prop names, case sensitive.
%     try
%         h.axes.XAxis = copy(p.Results.uiax.XAxis);
%         xzcopied(1) = true; 
%     catch
%         % No plan-B for now. Some XAxis properties will be copied later (see uiaxGoodParams).
%     end
%     try
%         h.axes.ZAxis = copy(p.Results.uiax.ZAxis);
%         xzcopied(2) = true; 
%     catch
%         % No plan-B for now. Some ZAxis properties will be copied later (see uiaxGoodParams). 
%     end
% 
%     % THIS LINE BELOW ALSO NEEDS TO BE ADDED LATER IN THE CODE
%     uiaxbadProps(ismember(uiaxbadProps.axis, xzcopyProp(xzcopied)),:) = []; %rm X|ZAxis if they were copied earlier



