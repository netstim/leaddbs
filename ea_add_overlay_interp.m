function cuts=ea_add_overlay_interp(filename,boundboxmm,cuts,tracor,options)
% This function overlays atlas data over 2d-slice views to be exported by
% LEAD-DBS. The function is called from ea_writeplanes.m
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn
keyboard
    set(0,'CurrentFigure',cuts)

    % load/generate atlas_index.mat
    
    
    slice=ea_sample_slice(filename,tracor,
    
planemm=[boundboxmm{1}(:),boundboxmm{2}(:),boundboxmm{3}(:),ones(length(boundboxmm{1}(:)),1)];
V=spm_vol(filename);
planevx=V.mat\planemm';

planevx=planevx(1:3,:);
planesample=meshgrid(planevx(1,:),planevx(2,:));
planeval=spm_sample_vol(V,planevx(1,:),planevx(2,:),planevx(3,:),1);

    
    
    
    %axis off
    
    % save table information
    save([options.earoot,'atlases',filesep,options.atlasset,filesep,'atlas_index.mat'],'atlases');


function sides=detsides(opt)

switch opt
    case 1 % left hemispheric atlas
        sides=1;
    case 2 % right hemispheric atlas
        sides=2;
    case 3
        sides=1:2;
    case 4
        sides=1:2;
    case 5
        sides=1; % midline
        
end

function in=ea_intersecdim(tracor)

switch tracor
    case 1
        in=3;
    case 2
        in=2;
    case 3
        in=1;
end

function pl=ea_planesdim(tracor)

switch tracor
    case 1
        pl=[1,2];
    case 2
        pl=[1,3];
    case 3
        pl=[2,3];
end


function fslice=ea_zeroframe(slice)

fslice=zeros(size(slice)+[2,2]);
fslice(2:end-1,2:end-1)=slice;


