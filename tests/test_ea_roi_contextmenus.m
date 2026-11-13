function tests = test_ea_roi_contextmenus
% Batch Java toolbar binding must preserve distinct ROIs with duplicate tags.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repo = fileparts(fileparts(mfilename('fullpath')));
oldPath = path;
addpath(repo,fullfile(repo,'classes','ea_roi'),fullfile(repo,'helpers'), ...
    fullfile(repo,'ext_libs','isMatlabVer'),fullfile(repo,'ext_libs','findjobj'));
testCase.addTeardown(@() path(oldPath));
end

function testBatchBindsEachDuplicateTagOnceAndRestoresTooltips(testCase)
assumeTrue(testCase,usejava('jvm') && usejava('awt'), ...
    'Java toolbar context menus require JVM and AWT.');
assumeTrue(testCase,isMatlabVer('<',[25,1]), ...
    'Java toolbar context menus apply before MATLAB R2025a.');

fig = figure('Visible','on','WindowStyle','normal','MenuBar','none', ...
    'ToolBar','none','NumberTitle','off','Name','Lead-DBS context-menu test', ...
    'Position',[20,20,220,160]);
testCase.addTeardown(@() delete(fig));
axes('Parent',fig);
toolbar = uitoolbar(fig);
rois = {cachedRoi(fig,toolbar,'palhy.nii.gz',[0.2,0.4,0.6]), ...
    cachedRoi(fig,toolbar,'palhy.mat',[0.6,0.4,0.2])};
for r = 1:numel(rois)
    testCase.addTeardown(@() delete(rois{r}));
    rois{r}.update_roi('initialize');
    rois{r}.breathelife(false);
end
expectedTooltips = cellfun(@(roi) roi.toggleH.TooltipString,rois, ...
    'UniformOutput',false);
verifyEqual(testCase,rois{1}.Tag,rois{2}.Tag);
verifyEqual(testCase,expectedTooltips{1},expectedTooltips{2});

% Measure only the binding call, preserving any pre-existing profiler data.
prior = profile('info');
previousStatus = profile('status');
profileCleanup = onCleanup(@() restoreProfiler(previousStatus));
profile resume
try
    ea_roi.bindcontextmenus(rois,fig);
catch exception
    profile off
    rethrow(exception)
end
profile off
measured = profile('info');
verifyEqual(testCase,callCount(measured,'findjobj')-callCount(prior,'findjobj'),1, ...
    'Batch binding must not fall back to a Java lookup for each ROI.');
clear profileCleanup

verifyEqual(testCase,cellfun(@(roi) roi.toggleH.TooltipString,rois, ...
    'UniformOutput',false),expectedTooltips);

% Inspect the resulting Java callbacks outside the profiled interval.
jtoggles = findjobj(fig,'class','ToggleButton','nomenu');
matched = zeros(size(rois));
for j = 1:numel(jtoggles)
    callback = get(jtoggles(j),'MouseReleasedCallback');
    if ~iscell(callback) || numel(callback) ~= 2 || ...
            ~isa(callback{1},'function_handle') || ~isa(callback{2},'ea_roi')
        continue
    end
    for r = 1:numel(rois)
        if callback{2} == rois{r}
            matched(r) = matched(r)+1;
            verifyEqual(testCase,char(jtoggles(j).getToolTipText),expectedTooltips{r});
        end
    end
end
verifyEqual(testCase,matched,[1,1], ...
    'Each distinct ROI must own exactly one toolbar right-click callback.');
end

function roi = cachedRoi(fig,toolbar,filename,color)
roi = ea_roi;
roi.name = 'palhy';
roi.Tag = 'palhy_right';
roi.niftiFilename = filename;
roi.color = color;
roi.binary = true;
roi.smooth = 0;
roi.hullsimplify = 0;
roi.threshold = 0.5;
roi.fv.vertices = [0,0,0;1,0,0;0,1,0;0,0,1];
roi.fv.faces = [1,2,3;1,2,4;1,3,4;2,3,4];
roi.sfv = roi.fv;
roi.plotFigureH = fig;
roi.htH = toolbar;
end

function count = callCount(info,name)
entries = info.FunctionTable;
if isempty(entries)
    count = 0;
else
    matching = strcmp({entries.FunctionName},name);
    count = sum([entries(matching).NumCalls]);
end
end

function restoreProfiler(status)
if strcmp(status.ProfilerStatus,'on')
    profile resume
else
    profile off
end
end
