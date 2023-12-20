function [Template,space,new_vx2mm,voxsize,neworigin_vx] = ea_createTemplateSpace(varargin)
TemplateCenter = varargin{1};
TemplateVXDimensions = varargin {2};
voxsize = varargin{3};
templatepath = varargin{4};
if length(varargin) == 5
    templatename = varargin{5};
else
    templatename = 'template.nii';
end
if ~isequal(size(TemplateCenter),[3,1])
    TemplateCenter=TemplateCenter';
    if ~isequal(size(TemplateCenter),[3,1])
        warning('TemplateCenter does not appear to be a 3d coordinate?')
        keyboard
    end
end
%% VTA04_GenerateTemplateSpace
% Generates a Templatespace with TemplateCenterin mm MNI coordinates and 
% voxel dimensions of TemplateVXDimensions.
%
% Voxels are identical to those of the MNI-template (org_space).
%
% Functions returns an empty space, voxsize and origin which can be used to
% create a Nifti using nii = make_nifti(space,voxsize,origin) which has
% exact overlap with the MNI template.

% TemplateCenter = [-13 -13 5]';
% TemplateVXDimensions = [80 80 80];

newcorner_mm = TemplateCenter - (((TemplateVXDimensions/2) .* voxsize)');
neworigin_vx = (-newcorner_mm)'./voxsize;
new_vx2mm = [   voxsize(1)     0        0          newcorner_mm(1);...
                0           voxsize(2)  0          newcorner_mm(2);...
                0              0      voxsize(3)   newcorner_mm(3);...
                0              0        0                   1       ];
space = single(zeros(TemplateVXDimensions));

%% build nii
cleannii=struct;
cleannii.fname = [templatepath,filesep,templatename];
cleannii.dim = TemplateVXDimensions;
cleannii.dt = 16;
cleannii.pinfo = [1;0;352];
cleannii.mat = new_vx2mm;
cleannii.n = [1,1];
cleannii.descrip = '';
cleannii.img=space;
cleannii.voxsize = voxsize;
cleannii.volnum = 1;
ea_write_nii(cleannii);
end