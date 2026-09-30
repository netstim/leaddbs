function cwin = ea_seegcortexcontrol(resultfig,baseColor,alpha)
% Controls for SEEG-style cortex visualization.

seegH = getappdata(resultfig,'seegcortex');

if isempty(seegH)
    cwin = [];
    return
end


% If already open, bring existing window forward

oldwin = getappdata(resultfig,'seegcortexcontrol');

if ~isempty(oldwin) && isgraphics(oldwin)
    figure(oldwin);
    cwin = oldwin;
    return
end


% Create window

cwin = figure( ...
    'Name','SEEG Cortex Settings', ...
    'NumberTitle','off', ...
    'MenuBar','none', ...
    'ToolBar','none', ...
    'Resize','off', ...
    'Position',[100 100 310 265], ...
    'Color',get(0,'defaultUicontrolBackgroundColor'));

setappdata(cwin,'resultfig',resultfig);


% Hemisphere visibility

uicontrol(cwin, ...
    'Style','text', ...
    'String','Hemisphere', ...
    'FontWeight','bold', ...
    'HorizontalAlignment','left', ...
    'Position',[20 220 100 22]);

rightCheck = uicontrol(cwin, ...
    'Style','checkbox', ...
    'String','Right', ...
    'Value',1, ...
    'Position',[110 220 75 22], ...
    'Callback',@(src,evt)setHemisphereVisible( ...
        resultfig,1,src.Value));

leftCheck = uicontrol(cwin, ...
    'Style','checkbox', ...
    'String','Left', ...
    'Value',1, ...
    'Position',[195 220 75 22], ...
    'Callback',@(src,evt)setHemisphereVisible( ...
        resultfig,2,src.Value));


% Color

uicontrol(cwin, ...
    'Style','text', ...
    'String','Color', ...
    'HorizontalAlignment','left', ...
    'Position',[20 165 80 22]);

colorPreview = uicontrol(cwin, ...
    'Style','pushbutton', ...
    'String','', ...
    'BackgroundColor',baseColor, ...
    'Position',[110 165 35 25], ...
    'Callback',@(src,evt)chooseCortexColor( ...
        resultfig,src));

uicontrol(cwin, ...
    'Style','pushbutton', ...
    'String','Choose Color...', ...
    'Position',[155 163 120 28], ...
    'Callback',@(src,evt)chooseCortexColor( ...
        resultfig,colorPreview));


% Alpha

uicontrol(cwin, ...
    'Style','text', ...
    'String','Alpha', ...
    'HorizontalAlignment','left', ...
    'Position',[20 110 80 22]);

alphaSlider = uicontrol(cwin, ...
    'Style','slider', ...
    'Min',0, ...
    'Max',1, ...
    'Value',alpha, ...
    'SliderStep',[0.01 0.10], ...
    'Position',[110 115 125 20]);

alphaEdit = uicontrol(cwin, ...
    'Style','edit', ...
    'String',sprintf('%.2f',alpha), ...
    'Position',[245 109 45 28]);

alphaSlider.Callback = ...
    @(src,evt)alphaSliderChanged( ...
        resultfig,src,alphaEdit);

alphaEdit.Callback = ...
    @(src,evt)alphaEditChanged( ...
        resultfig,src,alphaSlider);


% Restore defaults

uicontrol(cwin, ...
    'Style','pushbutton', ...
    'String','Restore Defaults', ...
    'Position',[90 40 130 32], ...
    'Callback',@(src,evt)restoreDefaults( ...
        resultfig, ...
        rightCheck, ...
        leftCheck, ...
        colorPreview, ...
        alphaSlider, ...
        alphaEdit));

end



% HEMISPHERE VISIBILITY

function setHemisphereVisible(resultfig,side,value)

seegH = getappdata(resultfig,'seegcortex');

if side > numel(seegH)
    return
end

if isempty(seegH{side}) || ~isgraphics(seegH{side})
    return
end

if value
    set(seegH{side},'Visible','on');
else
    set(seegH{side},'Visible','off');
end

end



% CHOOSE COLOR

function chooseCortexColor(resultfig,previewH)

color = ea_uisetcolor;

if numel(color) ~= 3
    return
end

set(previewH,'BackgroundColor',color);

applyCortexColor(resultfig,color);

end



% APPLY COLOR

function applyCortexColor(resultfig,color)

seegH = getappdata(resultfig,'seegcortex');

for s = 1:numel(seegH)

    if isempty(seegH{s}) || ~isgraphics(seegH{s})
        continue
    end

    ud = get(seegH{s},'UserData');

    if ~isstruct(ud) || ~isfield(ud,'GlowWeight')
        continue
    end

    w = ud.GlowWeight;

    % Preserve original dark-cortex → white-rim appearance
    C = color + w .* (1-color);

    set(seegH{s}, ...
        'FaceVertexCData',C, ...
        'FaceColor','interp');

    ud.BaseColor = color;
    set(seegH{s},'UserData',ud);
end

drawnow

end



% ALPHA SLIDER

function alphaSliderChanged(resultfig,sliderH,editH)

alpha = sliderH.Value;

editH.String = sprintf('%.2f',alpha);

setAlpha(resultfig,alpha);

end



% ALPHA TEXT BOX

function alphaEditChanged(resultfig,editH,sliderH)

alpha = str2double(editH.String);

if isnan(alpha) || alpha < 0 || alpha > 1
    editH.String = sprintf('%.2f',sliderH.Value);
    return
end

sliderH.Value = alpha;

setAlpha(resultfig,alpha);

end



% APPLY ALPHA

function setAlpha(resultfig,alpha)

seegH = getappdata(resultfig,'seegcortex');

for s = 1:numel(seegH)

    if ~isempty(seegH{s}) && isgraphics(seegH{s})
        set(seegH{s},'FaceAlpha',alpha);
    end
end

drawnow

end



% RESTORE DEFAULTS

function restoreDefaults( ...
    resultfig, ...
    rightCheck, ...
    leftCheck, ...
    colorPreview, ...
    alphaSlider, ...
    alphaEdit)

defaultColor = [0 0 0];
defaultAlpha = 0.4;

% Both hemispheres on
rightCheck.Value = 1;
leftCheck.Value = 1;

setHemisphereVisible(resultfig,1,1);
setHemisphereVisible(resultfig,2,1);

% Black base
set(colorPreview,'BackgroundColor',defaultColor);
applyCortexColor(resultfig,defaultColor);

% Alpha 0.4
alphaSlider.Value = defaultAlpha;
alphaEdit.String = sprintf('%.2f',defaultAlpha);

setAlpha(resultfig,defaultAlpha);

end