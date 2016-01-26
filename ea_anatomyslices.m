function ea_anatomyslices(resultfig,togglestates,options)
% input xyz in mm coordinates.
% this function plots an anatomical slice to the 3D viewer of lead-dbs.
% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn


set(0, 'currentfigure', resultfig);  % for figures
atlases=getappdata(resultfig,'atlases');
togglestates.xyztransparencies=double(togglestates.xyztransparencies/100);

xsliceplot=getappdata(resultfig,'xsliceplot');
try delete(xsliceplot); end
ysliceplot=getappdata(resultfig,'ysliceplot');
try delete(ysliceplot); end

zsliceplot=getappdata(resultfig,'zsliceplot');
try delete(zsliceplot); end

V=getappdata(resultfig,'V');
inverted=getappdata(resultfig,'inverted');
if isempty(inverted)
    inverted=0;
end
options.d2.writeatlases=1;
templateused=getappdata(resultfig,'templateused');
if ~strcmp(templateused,togglestates.template) || isempty(V) % reload image(s)
    clear V
    switch togglestates.template
        case 'MNI-Template'
            V{1}=spm_vol([options.earoot,'templates',filesep,'mni_hires.nii']);
            
        case 'Patient Post-OP'
            if options.native
                V{1}=spm_vol([options.root,options.patientname,filesep,options.prefs.tranii_unnormalized]);
                try
                    V{2}=spm_vol([options.root,options.patientname,filesep,options.prefs.cornii_unnormalized]);
                catch
                    try
                        V{2}=spm_vol([options.root,options.patientname,filesep,options.prefs.tranii_unnormalized]);
                    end
                end
                
                try
                    V{3}=spm_vol([options.root,options.patientname,filesep,options.prefs.sagnii_unnormalized]);
                catch
                    try
                        V{3}=spm_vol([options.root,options.patientname,filesep,options.prefs.tranii_unnormalized]);
                    end
                end
            
            else
            % load tra
            try
                V{1}=spm_vol([options.root,options.patientname,filesep,options.prefs.gtranii]);
            catch
                V{1}=spm_vol([options.root,options.patientname,filesep,options.prefs.tranii]);
            end
            
            % load cor
            try
                V{2}=spm_vol([options.root,options.patientname,filesep,options.prefs.gcornii]);
            catch
                try
                V{2}=spm_vol([options.root,options.patientname,filesep,options.prefs.cornii]);
                end
            end
            
            % load sag
            try
                V{3}=spm_vol([options.root,options.patientname,filesep,options.prefs.gsagnii]);
            catch
                try
                 V{3}=spm_vol([options.root,options.patientname,filesep,options.prefs.sagnii]);
                end
            end
            end
        case 'Patient Pre-OP'
            if options.native
                V{1}=spm_vol([options.root,options.patientname,filesep,options.prefs.prenii_unnormalized]);
            else
                
                % load tra
                try
                    V{1}=spm_vol([options.root,options.patientname,filesep,options.prefs.gprenii]);
                catch
                    V{1}=spm_vol([options.root,options.patientname,filesep,options.prefs.prenii]);
                end
            end
            
        case 'Choose...'
            
            V{1}=spm_vol(togglestates.customfile);
                        
    end
    setappdata(resultfig,'templateused',togglestates.template); % refresh used template.
end



if ~inverted==togglestates.tinvert
    inverted=togglestates.tinvert;
else
    
end
    
togglestates.xyzmm=[togglestates.xyzmm';1];
try
xyzv= V{1}.mat \ togglestates.xyzmm;
catch
    keyboard
end
xyzv=round(xyzv(1:3)); % now in voxel coordinates.

%colormap gray



if togglestates.xyztoggles(1)

    
    usesag=(length(V)>2)*2; % check if explicit saggital volume is available
    if inverted
        [~,slice]=ea_writeplanes(options, togglestates.xyzmm(1),3,V{1+usesag},'off', 0,atlases);
        slice=flipdim(permute(double(slice),[2,1,3]),2);
        %slice=flipdim(slice,1);
    else
        [xx,yy,zz]=meshgrid(xyzv(1),1:0.5:V{1+usesag}.dim(2),1:0.5:V{1+usesag}.dim(3));
        slice=spm_sample_vol(V{1+usesag},xx,yy,zz,1);
    end
    
    %slice=ea_invert(slice,inverted);
    imin=proxy_slice(slice,togglestates,1);
    clear bb
    bb(1,:)=[xyzv(1),V{1+usesag}.dim(2),V{1+usesag}.dim(3),1]; % upper left point of image in voxels
    bb(2,:)=[xyzv(1),0,V{1+usesag}.dim(3),1];
    bb(3,:)=[xyzv(1),V{1+usesag}.dim(2),0,1];
    bb(4,:)=[xyzv(1),0,0,1];
    
    bb(:,1:3)=bb(:,1:3);
    
    bb=V{1+usesag}.mat*bb'; % in mm
    bb=bb(1:3,:)';
    %zsliceplot=imsurf(imin,ulp,[1,0,0],[0,0,-1],scale);
    
    xsliceplot=surface('XData',[min(bb(:,1)) max(bb(:,1));min(bb(:,1)) max(bb(:,1))],...
        'YData',[min(bb(:,2)) min(bb(:,2));max(bb(:,2)), max(bb(:,2))],...
        'ZData',[min(bb(:,3)) max(bb(:,3)); min(bb(:,3)) max(bb(:,3))],...
        'CData',imin(:,:,1:3),...
        'FaceColor','texturemap','AlphaDataMapping','none','FaceAlpha',togglestates.xyztransparencies(1),'EdgeColor','none');
    %set(xsliceplot,'SpecularColorReflectance',0)
    set(xsliceplot,'SpecularStrength',0.1)
    bbmm{1}=linspace(bb(1,1),bb(4,1),20);
    bbmm{2}=linspace(bb(1,2),bb(2,2),20);
    bbmm{3}=linspace(bb(1,3),bb(3,3),20);
    %ea_add_overlay_3d(bbmm,resultfig,3,options);
    
    %surface('XData',[max(ulp(:,1)),min(ulp(:,1))],'YData',[min(ulp(:,2)),max(ulp(:,2))],'ZData',[min(ulp(:,3)),max(ulp(:,3))],'CData',imin(:,:,1:3));
    %catch
    %    disp('Z-Volume cut out of bounds.');
    %end
    
end

if togglestates.xyztoggles(2)
    
    % check whether second nii is being used:
    usecor=length(V)>1; % check if explicit coronar volume is available
    if inverted
        [~,slice]=ea_writeplanes(options, togglestates.xyzmm(2),2,V{1+usecor},'off', 0,atlases);
        slice=flipdim(permute(double(slice),[2,1,3]),2);
    else
        [xx,yy,zz]=meshgrid(1:0.5:V{1+usecor}.dim(1),xyzv(2),1:0.5:V{1+usecor}.dim(3));
        slice=flipud(squeeze((spm_sample_vol(V{1+usecor},squeeze(xx),squeeze(yy),squeeze(zz),2))));
    end
    %slice=ea_invert(slice,inverted);
    
    %slice=flipud(squeeze(double(nii{1+usecor}.img(:,xyzv(2),:))));
    imin=proxy_slice(slice,togglestates,2);
    clear bb
    bb(1,:)=[V{1+usecor}.dim(1),xyzv(2),V{1+usecor}.dim(3),1]; % upper left point of image in voxels
    bb(2,:)=[0,xyzv(2),V{1+usecor}.dim(3),1]; % upper left point of image in voxels
    bb(3,:)=[V{1+usecor}.dim(1),xyzv(2),0,1]; % upper left point of image in voxels
    bb(4,:)=[0,xyzv(2),0,1]; % upper left point of image in voxels
    
    bb(:,1:3)=bb(:,1:3);
    
    bb=V{1+usecor}.mat*bb'; % in mm
    bb=bb(1:3,:)';
    
    
    ysliceplot=surface('XData',[min(bb(:,1)) min(bb(:,1));max(bb(:,1)) max(bb(:,1))],...
        'YData',[min(bb(:,2)) max(bb(:,2));min(bb(:,2)), max(bb(:,2))],...
        'ZData',[min(bb(:,3)) max(bb(:,3)); min(bb(:,3)) max(bb(:,3))],...
        'CData',imin(:,:,1:3),...
        'FaceColor','texturemap','AlphaDataMapping','none','FaceAlpha',togglestates.xyztransparencies(2),'EdgeColor','none');
    set(ysliceplot,'SpecularStrength',0.1)
    set(ysliceplot,'DiffuseStrength',0.5)
    set(ysliceplot,'AmbientStrength',1)

bbmm{1}=linspace(bb(1,1),bb(4,1),20);
bbmm{2}=linspace(bb(1,2),bb(2,2),20);
bbmm{3}=linspace(bb(1,3),bb(3,3),20);
%ea_add_overlay_3d(bbmm,resultfig,2,options);
    %ysliceplot=imsurf(imin,ulp,[0,1,0],[-1,0,0],scale);
    %catch
    %    disp('Y-Volume cut out of bounds.');
    %end
end

if togglestates.xyztoggles(3)
    %try
    %imsurf(repmat(uint8((squeeze(fliplr(nii.img(:,:,graphopts.volxyz(3)))'*255)/max(nii.img(:)))),[1,1,2]),[size(nii.img,1)+0.5,size(nii.img,2)+0.5,graphopts.volxyz(3)],[0,0,1],[-1,0,0],1)
    %original: imin=repmat(uint8((squeeze(fliplr(nii.img(:,:,xyzv(1)))'*255)/max(nii.img(:)))),[1,1,4]);
    
    if inverted
        [~,slice]=ea_writeplanes(options, togglestates.xyzmm(3),1,V{1},'off', 0,atlases);
        if V{1}.mat(1)<0
        slice=flip(permute(double(slice),[2,1,3]),2);
        else
        slice=permute(double(slice),[2,1,3]);
        end
    else
        [xx,yy,zz]=meshgrid(1:0.2:V{1}.dim(1),1:0.2:V{1}.dim(2),xyzv(3));
        if V{1}.mat(1)<0
        slice=flip(spm_sample_vol(V{1},xx,yy,zz,4)',1);
        else
        slice=spm_sample_vol(V{1},xx,yy,zz,4)';            
        end
    end
    %slice=ea_invert(slice,inverted);
    %slice=flipud(squeeze(double((nii{1}.img(:,:,xyzv(3))))));
    imin=proxy_slice(slice,togglestates,3);
        clear bb
    bb(1,:)=[V{1}.dim(1),V{1}.dim(2),xyzv(3),1]'; % upper left point of image in voxels
    bb(2,:)=[0,V{1}.dim(2),xyzv(3),1]'; % upper left point of image in voxels
    bb(3,:)=[V{1}.dim(1),0,xyzv(3),1]'; % upper left point of image in voxels
    bb(4,:)=[0,0,xyzv(3),1]'; % upper left point of image in voxels

    bb(:,1:3)=bb(:,1:3);
    bb(:,1)=bb(:,1);
    bb=V{1}.mat*bb'; % in mm
    bb=bb(1:3,:)';
    
    zsliceplot=surface('XData',[min(bb(:,1)) min(bb(:,1));max(bb(:,1)) max(bb(:,1))],...
   'YData',[min(bb(:,2)) max(bb(:,2));min(bb(:,2)), max(bb(:,2))],...
   'ZData',[min(bb(:,3)) max(bb(:,3));min(bb(:,3)) max(bb(:,3))],...
   'CData',imin(:,:,1:3),...
    'FaceColor','texturemap','AlphaDataMapping','none','FaceAlpha',togglestates.xyztransparencies(3),'EdgeColor','none');    
set(zsliceplot,'SpecularStrength',0)
set(zsliceplot,'DiffuseStrength',0.5)
set(zsliceplot,'AmbientStrength',0.3)

    

%xsliceplot=imsurf(imin,bb,[0,0,1],[-1,0,0],scale);
    bbmm{1}=linspace(bb(1,1),bb(4,1),20);
    bbmm{2}=linspace(bb(1,2),bb(2,2),20);
    bbmm{3}=linspace(bb(1,3),bb(3,3),20);
    %ea_add_overlay_3d(bbmm,resultfig,1,options);
    %catch
    %    disp('X-Volume cut out of bounds.');
    %end
end

% store data in figure
setappdata(resultfig,'xsliceplot',xsliceplot);
setappdata(resultfig,'ysliceplot',ysliceplot);
setappdata(resultfig,'zsliceplot',zsliceplot);
setappdata(resultfig,'V',V);
setappdata(resultfig,'inverted',inverted);




function slice=ea_invert(slice,flag)
if flag
slice=slice*-1+max(slice(:));
end




function imin=proxy_slice(slice,togglestates,dim)
maxv=max(slice(:));
minv=min(slice(:));

slice=slice-minv; % 0 smallest number.
slice=(slice/(maxv-minv))*255; % 255 highest number.
imin=repmat(uint8((((slice)))),[1,1,4]);
imin(:,:,4)=uint8(togglestates.xyztransparencies(dim));

