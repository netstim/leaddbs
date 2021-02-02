function ea_axon_viz(axons, resultfig)

[~, fname] = fileparts(axons);
axons = load(axons, 'fibers');

set(0, 'CurrentFigure', resultfig);

% Check and creat toolbox if needed
PL = getappdata(resultfig,'PL');
if isempty(PL) || ~isfield(PL, 'ht') || ~isvalid(PL.ht)
    PL.ht = uitoolbar(resultfig);
end

toolbar = PL.ht;

if isfield(PL, 'axon')
    axonInd = length(PL.axon);
else
    PL.axon = cell(3,1);
    axonInd = 0;
end

prefs = ea_prefs;

% Activated fibers
Ind = axons.fibers(:,5)==1;
fibers = axons.fibers(Ind,1:3);
if ~isempty(fibers)
    [~,~,idx] = unique(axons.fibers(Ind,4));
    idx = accumarray(idx,1);
    PL.axon{axonInd+1} = showAxons(fibers, idx, prefs.d3.axon_activated_color, [fname, '_Activated'], toolbar);
else
    fprintf('\n')
    warning('off', 'backtrace');
    warning('No activated fiber found!\n');
    warning('on', 'backtrace');
end

% Non-activated fibers
Ind = axons.fibers(:,5)==0;
fibers = axons.fibers(Ind,1:3);
if ~isempty(fibers)
    [~,~,idx] = unique(axons.fibers(Ind,4));
    idx = accumarray(idx,1);
    PL.axon{axonInd+2} = showAxons(fibers, idx, prefs.d3.axon_nonactivated_color, [fname, '_Nonactivated'], toolbar);
else
    fprintf('\n')
    warning('off', 'backtrace');
    warning('No non-activated fiber found!\n');
    warning('on', 'backtrace');
end

% Damaged fibers
Ind = axons.fibers(:,5)==-1;
fibers = axons.fibers(Ind,1:3);
if ~isempty(fibers)
    [~,~,idx] = unique(axons.fibers(Ind,4));
    idx = accumarray(idx,1);
    PL.axon{axonInd+3} = showAxons(fibers, idx, prefs.d3.axon_damaged_color, [fname, '_Damaged'], toolbar);
else
    fprintf('\n')
    warning('off', 'backtrace');
    warning('No damaged fiber found!\n');
    warning('on', 'backtrace');
end

setappdata(resultfig, 'PL', PL);


function objs = showAxons(fibers, idx, color, name, toolbar)
objs = ea_showfiber(fibers, idx, color);
axis fill;

uitoggletool(toolbar,'CData',ea_get_icn('fibers'),...
    'TooltipString',name,...
    'OnCallback',{@ea_atlasvisible,objs},...
    'OffCallback',{@ea_atlasinvisible,objs},'State','on');
