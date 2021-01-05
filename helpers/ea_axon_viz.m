function ea_axon_viz(axons, resultfig)

[~, fname] = fileparts(axons);
axons = load(axons, 'fibers');

% Check and creat toolbox if needed
PL = getappdata(resultfig,'PL');
if isempty(PL) || ~isfield(PL, 'ht') || ~isvalid(PL.ht)
    PL.ht = uitoolbar(resultfig);
    setappdata(resultfig, 'PL', PL);
end

toolbar = PL.ht;

prefs = ea_prefs;

% Activated fibers
Ind = axons.fibers(:,5)==1;
fibers = axons.fibers(Ind,1:3);
[~,~,idx] = unique(axons.fibers(Ind,4));
idx = accumarray(idx,1);
showAxons(fibers, idx, prefs.d3.axon_activated_color, [fname, '_Activated'], toolbar);

% Non-activated fibers
Ind = axons.fibers(:,5)==0;
fibers = axons.fibers(Ind,1:3);
[~,~,idx] = unique(axons.fibers(Ind,4));
idx = accumarray(idx,1);
showAxons(fibers, idx, prefs.d3.axon_nonactivated_color, [fname, '_Nonactivated'], toolbar);

% Damaged fibers
Ind = axons.fibers(:,5)==-1;
fibers = axons.fibers(Ind,1:3);
[~,~,idx] = unique(axons.fibers(Ind,4));
idx = accumarray(idx,1);
showAxons(fibers, idx, prefs.d3.axon_damaged_color, [fname, '_Damaged'], toolbar);


function showAxons(fibers, idx, color, name, toolbar)
objs = ea_showfiber(fibers, idx, color);
axis fill;

uitoggletool(toolbar,'CData',ea_get_icn('fibers'),...
    'TooltipString',name,...
    'OnCallback',{@ea_atlasvisible,objs},...
    'OffCallback',{@ea_atlasinvisible,objs},'State','on');
