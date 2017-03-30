function [mask_matrix,mask_indices] = ea_gethull(grayfilename,LR,varargin)
% Modified from Hermes method
%
% Inputs -
%   grayfilename    - t1_class = gray matter outline
%   outputdir       - save location
%   varargin(1)     - smooth level
%   varargin(2)     - br_cutoff
%
% function [mask_matrix,mask_indices]=get_hull(FWHM for smoothing, brain_cutoff)
% default FWHM for smoothing = 13
% default brain_cutoff = 0.3 (of smoothed images)
%
% Sample:
%       grayfilename = '/Users/edk26/Dropbox/Research/Programs/electrode_imaging/Fluoro_Imaging/Subjects/RZ_FS/mri/t1_class.nii'
%       outputdir='hull';
%       [~,~] = ea_gethull(grayfilename,outputdir,3,21,.3);
% 
% Dependencies - hollow_brain
%
% this function extracts a surface map of the brain from the gray matter
% input (in the structure "data"), thresholding at value "threshold"
% make sure that data has the fields designated in the
%
%     Copyright (C) 2009  D. Hermes & K.J. Miller, Dept of Neurology and Neurosurgery, University Medical Center Utrecht
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
% use spm_select to get image file 
data.grayfilename=grayfilename;

brain_info=spm_vol([data.grayfilename]); 
[g]=spm_read_vols(brain_info);

if isequal(LR,'r')
    % select right brain:
    g(g==2)=0; 
    g(g==3)=0; 
elseif isequal(LR,'l') % select left brain:
    g(g==41)=0; 
    g(g==42)=0; 
else
end
g(g>0)=1;

if length(varargin)==0
    sm_lvl=6; %smoothing parameter for rendering gray
    br_cutoff=.1; %cutoff for thresholding
elseif length(varargin)==1
    sm_lvl=varargin{1};
    br_cutoff=.1; %cutoff for thresholding
elseif length(varargin)==2
    sm_lvl=varargin{1};
    br_cutoff=varargin{2};
end

%lightly smooth gray matter
a=g;
a=smooth3(a,'gaussian',[sm_lvl sm_lvl sm_lvl]);
a=smooth3(a,'box',[sm_lvl sm_lvl sm_lvl]);

%identifies "enclosed points" for later removal
a=a>br_cutoff; %clear g w %combination of grey and white matter
a=hollow_brain(a);

bwa=bwlabeln(a);
% [x,y,z]=ind2sub(size(bwa),...
%     find(bwa==2));
brainsize=length(a(:));
size4surface=[brainsize/10 brainsize/500]; %set required size for surface
for k=1:max(max(max(bwa)))
    if length(find(bwa==k))<size4surface(1) && length(find(bwa==k))>size4surface(2)
        a=bwa==k;
        disp('nice surface found');
        break;
    end
end
if max(max(max(a)))>1
    disp('no good surface representation found, change size of surface in get_mask');
end
%%%%
% outputdir= spm_select(1,'dir','select directory to save surface');
a=double(a);

[x, y, z] = ind2sub(size(a),find(a));
gs=([x y z]*brain_info.mat(1:3,1:3)')+repmat(brain_info.mat(1:3,4),1,length(x))'; 

mask_indices    =([x y z]*brain_info.mat(1:3,1:3)')+repmat(brain_info.mat(1:3,4),1,length(x))'; 
mask_matrix     = a;
