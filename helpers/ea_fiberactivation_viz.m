function ea_fiberactivation_viz(fiberActivation, resultfig)

[~, fname] = fileparts(fiberActivation);
fiberActivation = load(fiberActivation, 'fibers');

set(0, 'CurrentFigure', resultfig);

% Check and creat toolbox if needed
PL = getappdata(resultfig,'PL');
if isempty(PL) || ~isfield(PL, 'ht') || ~isvalid(PL.ht)
    PL.ht = uitoolbar(resultfig);
end

toolbar = PL.ht;

if isfield(PL, 'fiberActivation')
    fiberInd = length(PL.fiberActivation);
else
    PL.fiberActivation = cell(3,1);
    fiberInd = 0;
end

prefs = ea_prefs;

% Activated fibers
Ind = fiberActivation.fibers(:,5)==1;
fibers = fiberActivation.fibers(Ind,1:3);
if ~isempty(fibers)
    [~,~,idx] = unique(fiberActivation.fibers(Ind,4));
    idx = accumarray(idx,1);
    PL.fiberActivation{fiberInd+1} = showFibers(fibers, idx, prefs.d3.fiber_activated_color, [fname, '_Activated'], toolbar);
else
    fprintf('\n')
    warning('off', 'backtrace');
    warning('No activated fiber found!');
    warning('on', 'backtrace');
end

% Non-activated fibers
Ind = fiberActivation.fibers(:,5)==0;
fibers = fiberActivation.fibers(Ind,1:3);
if ~isempty(fibers)
    [~,~,idx] = unique(fiberActivation.fibers(Ind,4));
    idx = accumarray(idx,1);
    PL.fiberActivation{fiberInd+2} = showFibers(fibers, idx, prefs.d3.fiber_nonactivated_color, [fname, '_Nonactivated'], toolbar);
else
    fprintf('\n')
    warning('off', 'backtrace');
    warning('No non-activated fiber found!');
    warning('on', 'backtrace');
end

% Damaged fibers (encapsulation)
Ind = fiberActivation.fibers(:,5)==-1;
fibers = fiberActivation.fibers(Ind,1:3);
if ~isempty(fibers)
    [~,~,idx] = unique(fiberActivation.fibers(Ind,4));
    idx = accumarray(idx,1);
    PL.fiberActivation{fiberInd+3} = showFibers(fibers, idx, prefs.d3.fiber_damaged_color, [fname, '_Damaged'], toolbar);
else
    fprintf('\n')
    warning('off', 'backtrace');
    warning('No damaged fiber found!');
    warning('on', 'backtrace');
end


% neuron passes through CSF
Ind = fiberActivation.fibers(:,5)==-2;
fibers = fiberActivation.fibers(Ind,1:3);
if ~isempty(fibers)
    [~,~,idx] = unique(fiberActivation.fibers(Ind,4));
    idx = accumarray(idx,1);
    PL.fiberActivation{fiberInd+4} = showFibers(fibers, idx, prefs.d3.fiber_csf_color, [fname, '_CSF'], toolbar);
else
    fprintf('\n')
    warning('off', 'backtrace');
    warning('No fibers in CSF found!');
    warning('on', 'backtrace');
end

% Outside of the domain (or intersect with the electrode avoiding encapsulation, happens due to sparse sampling)
Ind = fiberActivation.fibers(:,5)==-3;
fibers = fiberActivation.fibers(Ind,1:3);
if ~isempty(fibers)
    [~,~,idx] = unique(fiberActivation.fibers(Ind,4));
    idx = accumarray(idx,1);
    PL.fiberActivation{fiberInd+5} = showFibers(fibers, idx, prefs.d3.fiber_outside_color, [fname, '_outside'], toolbar);
else
    fprintf('\n')
    warning('off', 'backtrace');
    warning('No fibers outside of the domain found!');
    warning('on', 'backtrace');
end

setappdata(resultfig, 'PL', PL);


function objs = showFibers(fibers, idx, color, name, toolbar)
objs = ea_showfiber(fibers, idx, color);
axis fill;

uitoggletool(toolbar,'CData',ea_get_icn('fibers'),...
    'TooltipString',name,...
    'OnCallback',{@(src, evt) ea_atlasvisible(objs)},...
    'OffCallback',{@(src, evt) ea_atlasinvisible(objs)},...
    'State','on',...
    'UserData','fiberactivation');
