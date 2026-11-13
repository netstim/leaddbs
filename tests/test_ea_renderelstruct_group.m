function tests = test_ea_renderelstruct_group
% Group rendering uses saved coordinates even without individual BIDS data.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repo = fileparts(fileparts(mfilename('fullpath')));
oldPath = path;
addpath(repo, fullfile(repo,'helpers'), fullfile(repo,'helpers','gui'), ...
    fullfile(repo,'classes','ea_trajectory'), ...
    fullfile(repo,'templates','electrode_models'), ...
    fullfile(repo,'ext_libs'), fullfile(repo,'ext_libs','nan'));
testCase.addTeardown(@() path(oldPath));
end

function setup(testCase)
folder = tempname;
mkdir(folder);
testCase.addTeardown(@() rmdir(folder,'s'));
fig = figure('Visible','off','MenuBar','none','ToolBar','none');
testCase.addTeardown(@() delete(fig));
axes('Parent',fig);
hold on;
testCase.TestData.folder = folder;
testCase.TestData.fig = fig;
end

function testMissingBidsUsesSavedSidesAndCoordinates(testCase)
[options,elstruct] = fixture(testCase,2);
patient = options.patient_list{1};
verifyFalse(testCase,isfolder(patient));
[rendered,~,sides] = ea_renderelstruct(options,testCase.TestData.fig,elstruct,1);
verifyEqual(testCase,sides,2);
verifyEqual(testCase,numel(rendered),1);
verifyEqual(testCase,rendered.side,2);
verifyEqual(testCase,rendered.options.patientname,'sub-offline');
verifyEqual(testCase,rendered.elstruct.coords_mm,elstruct.coords_mm);
verifyPoint(testCase,rendered.elpatch,elstruct.coords_mm{2}(2,:));
% Rendering a standalone group must not create a reconstruction directory.
verifyFalse(testCase,isfolder(patient));
end

function testMirroredUnilateralUsesItsSavedContralateralSide(testCase)
[options,original] = fixture(testCase,1);
mirrored = original;
mirrored.name = [original.name '_mirrored'];
mirrored.coords_mm = {[], original.coords_mm{1}.*[-1 1 1]};
mirrored.trajectory = {[], original.trajectory{1}.*[-1 1 1]};
mirrored.markers(2) = original.markers(1);
for field = {'head','tail','x','y'}
    mirrored.markers(2).(field{1})(1) = -mirrored.markers(2).(field{1})(1);
    mirrored.markers(1).(field{1}) = [];
end
options.sides = 1; % The on-disk/original hemisphere must not hide its mirror.
options.patient_list = [options.patient_list; options.patient_list];
options.d3.mirrorsides = 1;
[rendered,~,sides] = ea_renderelstruct(options,testCase.TestData.fig, ...
    [original mirrored],2);
verifyEqual(testCase,sides,2);
verifyEqual(testCase,numel(rendered),1);
verifyEqual(testCase,rendered.side,2);
verifyPoint(testCase,rendered.elpatch,mirrored.coords_mm{2}(2,:));
end

function testDiskReconstructionCannotOverrideGroupSidesOrBeRewritten(testCase)
[options,elstruct] = fixture(testCase,2);
patient = options.patient_list{1};
mkdir(fullfile(patient,'reconstruction'));
filename = fullfile(patient,'reconstruction','sub-offline_desc-reconstruction.mat');
reco.native.markers = struct('head',[20 -10 -5]);
% The disk reconstruction deliberately disagrees with the saved group: its
% native data only has a right electrode, while the group has only a left.
save(filename,'reco','-v7');
before = readBytes(filename);
[rendered,~,sides] = ea_renderelstruct(options,testCase.TestData.fig,elstruct,1);
verifyEqual(testCase,sides,2);
verifyEqual(testCase,rendered.side,2);
verifyPoint(testCase,rendered.elpatch,elstruct.coords_mm{2}(2,:));
verifyEqual(testCase,readBytes(filename),before, ...
    'Displaying cached group data must not rewrite individual reconstructions.');
end

function testColumnCellCoordinatesRenderBothSides(testCase)
[options,right] = fixture(testCase,1);
[~,left] = fixture(testCase,2);
elstruct = right;
elstruct.coords_mm = {right.coords_mm{1}; left.coords_mm{2}};
elstruct.trajectory = {right.trajectory{1}; left.trajectory{2}};
elstruct.activecontacts = {right.activecontacts{1}; left.activecontacts{2}};
elstruct.markers(2) = left.markers(2);
[rendered,~,sides] = ea_renderelstruct(options,testCase.TestData.fig,elstruct,1);
verifyEqual(testCase,sides,[1 2]);
verifyEqual(testCase,[rendered.side],[1 2]);
verifyPoint(testCase,rendered(1).elpatch,elstruct.coords_mm{1}(2,:));
verifyPoint(testCase,rendered(2).elpatch,elstruct.coords_mm{2}(2,:));
end

function [options,elstruct] = fixture(testCase,side)
options.leadprod = 'group';
options.native = 0;
options.patient_list = {fullfile(testCase.TestData.folder,'sub-offline')};
options.root = [testCase.TestData.folder filesep];
options.patientname = 'group';
options.sides = [1 2];
options.elmodel = 'Medtronic 3389';
options.d3.elrendering = 3;
options.d3.prolong_electrode = 0;
options.d3.mirrorsides = 0;
options.d3.hlactivecontacts = 1;
options.d3.showactivecontacts = 1;
options.d3.showpassivecontacts = 0;
options.d3.colorpointcloud = 0;
options.d3.isomatrix = {};
options.prefs.d3.pointcloudstyle = 'plain';
options = ea_resolve_elspec(options);
elstruct.name = 'sub-offline';
elstruct.elmodel = options.elmodel;
elstruct.coords_mm = {[],[]};
elstruct.trajectory = {[],[]};
elstruct.activecontacts = {zeros(4,1),zeros(4,1)};
elstruct.activecontacts{side}(2) = 1;
coords = [(3-2*side)*20*ones(4,1), -10*ones(4,1), (-5:2:1)'];
elstruct.coords_mm{side} = coords;
elstruct.trajectory{side} = [coords(1,:);coords(end,:)+[0 0 10]];
blank = struct('head',[],'tail',[],'x',[],'y',[]);
elstruct.markers = repmat(blank,1,2);
elstruct.markers(side).head = coords(1,:);
elstruct.markers(side).tail = coords(end,:);
elstruct.markers(side).x = coords(1,:)+[0.65 0 0];
elstruct.markers(side).y = coords(1,:)+[0 0.65 0];
end

function verifyPoint(testCase,points,expected)
assertEqual(testCase,numel(points),1);
verifyEqual(testCase,[points.XData points.YData points.ZData],expected);
end

function bytes = readBytes(filename)
fid = fopen(filename,'r');
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
bytes = fread(fid,Inf,'*uint8');
end
