function tests=test_ea_atlaslabels
% Presets and missing labels must not break scene-indexed label visibility.
tests=functiontests(localfunctions);
end

function setupOnce(testCase)
repo=fileparts(fileparts(mfilename('fullpath')));
previousPath=path;
addpath(fullfile(repo,'helpers','gui'));
testCase.addTeardown(@() path(previousPath));
end

function testInitialHideIgnoresNonTextEntries(testCase)
[~,~,labels]=fixture(testCase);
labels(1).Visible='on';
labels(3).Visible='on';
labels(4).Visible='on';
ea_atlaslabelvisible([],[],labels,'off');
verifyEqual(testCase,labelStates(labels),{'off';'off';'off'});
end

function testOnUsesSceneIndicesAndExplicitButtons(testCase)
[labelbutton,buttons,labels]=fixture(testCase);
buttons(1).State='off';
buttons(2).State='on'; % This scene row intentionally has no Text label.
buttons(3).State='on';
buttons(4).State='off';
ea_atlaslabelvisible(labelbutton,[],labels,'on');
verifyEqual(testCase,labelStates(labels),{'off';'on';'off'});

% Switching a preset changes the exact scene toggles, including labels that
% were absent from the initial selection. All labels remain addressable.
buttons(1).State='on';
buttons(3).State='off';
buttons(4).State='on';
ea_atlaslabelvisible(labelbutton,[],labels,'on');
verifyEqual(testCase,labelStates(labels),{'on';'off';'on'});
end

function testOffHidesEveryLabelWithoutChangingToggles(testCase)
[labelbutton,buttons,labels]=fixture(testCase);
set(buttons,'State','on');
ea_atlaslabelvisible(labelbutton,[],labels,'on');
ea_atlaslabelvisible(labelbutton,[],labels,'off');
verifyEqual(testCase,labelStates(labels),{'off';'off';'off'});
verifyTrue(testCase,all(arrayfun(@(button) strcmp(button.State,'on'),buttons)));
end

function testMissingToggleLeavesItsLabelHidden(testCase)
[labelbutton,buttons,labels,fig]=fixture(testCase);
setappdata(fig,'colorbuttons',buttons(1:3));
labels(4).Visible='on';
ea_atlaslabelvisible(labelbutton,[],labels,'on');
verifyEqual(testCase,char(labels(4).Visible),'off');
end

function [labelbutton,buttons,labels,fig]=fixture(testCase)
fig=figure('Visible','off');
testCase.addTeardown(@() delete(fig));
ax=axes('Parent',fig);
labels=gobjects(1,4);
labels(1)=text(ax,0,0,0,'First','Visible','on');
labels(3)=text(ax,0,0,0,'Third','Visible','on');
labels(4)=text(ax,0,0,0,'Fourth','Visible','on');

% Labels need not share a toolbar with the atlas toggles, and unrelated
% tools can appear in either toolbar without changing the scene mapping.
labeltoolbar=uitoolbar(fig);
uitoggletool(labeltoolbar,'State','on');
labelbutton=uitoggletool(labeltoolbar,'State','off');
uipushtool(labeltoolbar);
atlastoolbar=uitoolbar(fig);
buttons=gobjects(1,4);
for index=1:4
    buttons(index)=uitoggletool(atlastoolbar,'State','on');
end
uipushtool(atlastoolbar);
uitoggletool(atlastoolbar,'State','off');
setappdata(fig,'colorbuttons',buttons);
end

function states=labelStates(labels)
isLabel=arrayfun(@(label) isa(label,'matlab.graphics.primitive.Text'),labels);
states=arrayfun(@(label) char(label.Visible),labels(isLabel), ...
    'UniformOutput',false)';
end
