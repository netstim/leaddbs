function tests = test_ea_spec_atlas_uniform
% Constant-color atlas patches must retain their appearance without copies
% of an identical RGB triplet and unused scalar-opacity data per vertex.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repo = fileparts(fileparts(mfilename('fullpath')));
previousPath = path;
addpath(repo);
testCase.addTeardown(@() path(previousPath));
end

function testUniformPathKeepsMeshAndScalarAppearance(testCase)
[fig, patchH, mesh] = fixture(testCase);
color = [0.2,0.4,0.6];
map = jet(17);
colormap(fig,map);
setappdata(fig,'meshChanges',0);
verticesListener = addlistener(patchH,'Vertices','PostSet',@(~,~) increment(fig));
facesListener = addlistener(patchH,'Faces','PostSet',@(~,~) increment(fig));
testCase.addTeardown(@() delete(verticesListener));
testCase.addTeardown(@() delete(facesListener));

ea_spec_atlas(patchH,'palhy.mat',map,1,color);

verifyEqual(testCase,patchH.Vertices,mesh.vertices);
verifyEqual(testCase,patchH.Faces,mesh.faces);
verifyEqual(testCase,getappdata(fig,'meshChanges'),0);
verifyEqual(testCase,patchH.FaceColor,color);
verifyEqual(testCase,patchH.FaceAlpha,0.5);
verifyEqual(testCase,patchH.FaceLighting,'gouraud');
verifyEqual(testCase,patchH.EdgeColor,'none');
verifyEqual(testCase,patchH.LineStyle,'--');
verifyEqual(testCase,patchH.LineWidth,0.1);
verifyEqual(testCase,patchH.SpecularColorReflectance,1);
verifyEqual(testCase,patchH.SpecularExponent,3);
verifyEqual(testCase,patchH.SpecularStrength,0.3);
verifyEqual(testCase,patchH.DiffuseStrength,0.4);
verifyLessThanOrEqual(testCase,size(patchH.FaceVertexCData,1),1);
verifyLessThanOrEqual(testCase,size(patchH.FaceVertexAlphaData,1),1);
verifyEqual(testCase,colormap(fig),map);
end

function testFourArgumentPathRetainsInterpolatedColorsAndLegacyAlpha(testCase)
[fig,patchH,mesh] = fixture(testCase);
colors = [0.2,0.4,0.6;0.3,0.5,0.7;0.4,0.6,0.8;0.5,0.7,0.9];
patchH.FaceVertexCData = colors;
colormap(fig,jet(17));

ea_spec_atlas(patchH,'palhy.mat',jet(17),1);

verifyEqual(testCase,patchH.FaceColor,'interp');
verifyEqual(testCase,patchH.FaceVertexCData,colors);
verifyEqual(testCase,patchH.FaceVertexAlphaData,0.3*ones(4,1));
verifyEqual(testCase,patchH.FaceAlpha,0.5);
verifyEqual(testCase,colormap(fig),gray(17));
verifyEqual(testCase,patchH.Vertices,mesh.vertices);
verifyEqual(testCase,patchH.Faces,mesh.faces);
end

function testUniformPathPreservesSpecialCommandStyles(testCase)
[fig,~,mesh] = fixture(testCase);
color = [0.2,0.4,0.6];
commands = {'labeling','mask.nii','cortex.nii','vat'};
properties = {'FaceAlpha','EdgeColor','FaceLighting','LineStyle',...
    'LineWidth','SpecularColorReflectance','SpecularExponent',...
    'SpecularStrength','DiffuseStrength'};
for k = 1:numel(commands)
    original = patch('Parent',fig.CurrentAxes,'Faces',mesh.faces,...
        'Vertices',mesh.vertices,'FaceVertexCData',repmat(color,4,1));
    uniform = patch('Parent',fig.CurrentAxes,'Faces',mesh.faces,...
        'Vertices',mesh.vertices,'FaceColor',color);
    ea_spec_atlas(original,commands{k},gray(17),1);
    ea_spec_atlas(uniform,commands{k},gray(17),1,color);
    verifyEqual(testCase,get(uniform,properties),get(original,properties));
    verifyLessThanOrEqual(testCase,size(uniform.FaceVertexAlphaData,1),1);
    verifyLessThanOrEqual(testCase,size(uniform.FaceVertexCData,1),1);
    if strcmp(commands{k},'vat')
        verifyEqual(testCase,uniform.FaceColor,[0.9,0.1,0.1]);
    else
        verifyEqual(testCase,uniform.FaceColor,color);
    end
end
end

function [fig,patchH,mesh] = fixture(testCase)
fig = figure('Visible','off');
testCase.addTeardown(@() delete(fig));
ax = axes('Parent',fig);
mesh.vertices = [0,0,0;1,0,0;0,1,0;0,0,1];
mesh.faces = [1,2,3;1,2,4;1,3,4;2,3,4];
patchH = patch('Parent',ax,'Faces',mesh.faces,'Vertices',mesh.vertices,...
    'FaceColor',[0.2,0.4,0.6]);
end

function increment(fig)
setappdata(fig,'meshChanges',getappdata(fig,'meshChanges')+1);
end
