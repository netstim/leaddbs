function tests = test_ea_roi_rendering
tests = functiontests(localfunctions);


function setupOnce(testCase)
repo = fileparts(fileparts(mfilename('fullpath')));
oldPath = path;
addpath(repo,fullfile(repo,'classes','ea_roi'),fullfile(repo,'helpers'),...
    fullfile(repo,'ext_libs','isMatlabVer'),fullfile(repo,'ext_libs','nan'));
testCase.addTeardown(@() path(oldPath));


function testInitializeUsesStoredGeometryWithoutReadingVolume(testCase)
[roi,~] = fixture(testCase);
expected = roi.sfv;
roi.nii = []; % A cache hit must not need the source volume.
[~,geometryRebuilt] = roi.update_roi('initialize');

verifyFalse(testCase,geometryRebuilt);
verifyEqual(testCase,roi.fv,expected);
verifyEqual(testCase,roi.sfv,expected);
verifyEqual(testCase,roi.patchH.Faces,expected.faces);
verifyEqual(testCase,roi.patchH.Vertices,expected.vertices);
verifyEqual(testCase,roi.toggleH.Tag,roi.Tag);


function testAppearanceChangesDoNotUploadGeometry(testCase)
[roi,fig] = fixture(testCase);
roi.update_roi('initialize');
roi.breathelife(false);
listenForMeshChanges(testCase,roi,fig);

roi.Visible = 'off';
roi.alpha = 0.25;
roi.color = [0.8,0.2,0.1];
roi.edgecolor = [0.1,0.1,0.1];
roi.SpecularStrength = 0.6;
roi.SpecularExponent = 7;
roi.SpecularColorReflectance = 0.2;
roi.DiffuseStrength = 0.3;
roi.AmbientStrength = 0.4;

verifyEqual(testCase,getappdata(fig,'meshUpdates'),0);
verifyEqual(testCase,char(roi.patchH.Visible),'off');
verifyEqual(testCase,char(roi.toggleH.State),'off');
verifyEqual(testCase,roi.patchH.FaceAlpha,0.25);
verifyEqual(testCase,roi.patchH.FaceColor,roi.color);
verifyEqual(testCase,roi.patchH.EdgeColor,roi.edgecolor);
verifyEqual(testCase,roi.patchH.SpecularStrength,0.6);
verifyEqual(testCase,roi.patchH.SpecularExponent,7);
verifyEqual(testCase,roi.patchH.SpecularColorReflectance,0.2);
verifyEqual(testCase,roi.patchH.DiffuseStrength,0.3);
verifyEqual(testCase,roi.patchH.AmbientStrength,0.4);


function testBreathelifeInstallsOneListenerPerProperty(testCase)
[roi,fig] = fixture(testCase);
roi.update_roi('initialize');
roi.breathelife(false);
roi.breathelife(false);
setappdata(fig,'visibleUpdates',0);
listener = addlistener(roi.patchH,'Visible','PostSet',...
    @(~,~) incrementCount(fig,'visibleUpdates'));
testCase.addTeardown(@() delete(listener));

roi.Visible = 'off';

verifyEqual(testCase,getappdata(fig,'visibleUpdates'),1);


function testMissingMeshRebuildsAndThresholdStillUpdates(testCase)
[roi,~] = fixture(testCase);
roi.fv = [];
roi.sfv = [];
[~,geometryRebuilt] = roi.update_roi('initialize');
verifyTrue(testCase,geometryRebuilt);
expected = roi.patchH.Vertices;
verifyNotEmpty(testCase,expected);
verifyEqual(testCase,roi.patchH.Faces,roi.sfv.faces);
roi.breathelife(false);

roi.threshold = 0.75;

verifyNotEqual(testCase,roi.patchH.Vertices,expected);
verifyEqual(testCase,roi.patchH.Vertices,roi.sfv.vertices);


function testColorModeChangeKeepsGeometry(testCase)
[roi,fig] = fixture(testCase);
roi.binary = false;
roi.colormap = jet(64);
roi.nii.img(1:3,3:6,3:6) = 1; % Include caps with values above the isolevel.
roi.update_roi('all');
roi.breathelife(false);
expected = roi.sfv;
listenForMeshChanges(testCase,roi,fig);

roi.usesolidcolor = false;

verifyEqual(testCase,roi.sfv,expected);
verifyEqual(testCase,getappdata(fig,'meshUpdates'),0);
verifyEqual(testCase,roi.patchH.FaceColor,'interp');
verifySize(testCase,roi.patchH.FaceVertexCData,[size(expected.vertices,1),3]);

roi.usesolidcolor = true;
verifyEqual(testCase,roi.patchH.FaceColor,roi.color);
verifyEqual(testCase,getappdata(fig,'meshUpdates'),0);


function testConstantIntensityHasValidColormapIndices(testCase)
[roi,~] = fixture(testCase);
roi.binary = false;
roi.usesolidcolor = false;
roi.colormap = jet(16);
roi.update_roi('all');

verifyTrue(testCase,all(isfinite(roi.patchH.FaceVertexCData(:))));
verifySize(testCase,roi.patchH.FaceVertexCData,[size(roi.sfv.vertices,1),3]);


function testBinaryColorModeRemainsSolid(testCase)
[roi,fig] = fixture(testCase);
roi.update_roi('initialize');
roi.breathelife(false);
listenForMeshChanges(testCase,roi,fig);

roi.usesolidcolor = false;

verifyEqual(testCase,roi.patchH.FaceColor,roi.color);
verifyEqual(testCase,getappdata(fig,'meshUpdates'),0);


function testRuntimeListenerStateIsNotSerialized(testCase)
[roi,~] = fixture(testCase);
roi.update_roi('initialize');
roi.breathelife(false);
saved = roi.saveobj;
atlasSaved = ea_roi2struct(roi);

verifyFalse(testCase,isfield(saved,'listenersInitialized'));
verifyFalse(testCase,isfield(saved,'contextToggleH'));
verifyFalse(testCase,isfield(atlasSaved,'listenersInitialized'));
verifyFalse(testCase,isfield(atlasSaved,'contextToggleH'));


function [roi,fig] = fixture(testCase)
fig = figure('Visible','off');
roi = ea_roi;
testCase.addTeardown(@() delete(fig));
testCase.addTeardown(@() delete(roi));
roi.niftiFilename = 'cached_rendering_test.nii';
roi.name = 'cached_rendering_test';
roi.Tag = 'cached_rendering_test_right';
roi.color = [0.2,0.4,0.6];
roi.binary = true;
roi.smooth = 0;
roi.hullsimplify = 0;
roi.threshold = 0.5;
roi.nii.img = zeros(8,8,8);
roi.nii.img(3:6,3:6,3:6) = 1;
roi.nii.mat = eye(4);
roi.fv.vertices = [0,0,0;1,0,0;0,1,0;0,0,1];
roi.fv.faces = [1,2,3;1,2,4;1,3,4;2,3,4];
roi.sfv = roi.fv;
roi.plotFigureH = fig;
roi.htH = uitoolbar(fig);


function listenForMeshChanges(testCase,roi,fig)
setappdata(fig,'meshUpdates',0);
verticesListener = addlistener(roi.patchH,'Vertices','PostSet',...
    @(~,~) incrementCount(fig,'meshUpdates'));
facesListener = addlistener(roi.patchH,'Faces','PostSet',...
    @(~,~) incrementCount(fig,'meshUpdates'));
testCase.addTeardown(@() delete(verticesListener));
testCase.addTeardown(@() delete(facesListener));


function incrementCount(fig,key)
setappdata(fig,key,getappdata(fig,key)+1);
