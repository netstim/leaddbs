function ea_anatomyslices(resultfig,togglestates,options)
% input xyz in mm coordinates.


set(0, 'currentfigure', resultfig);  % for figures

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

templateused=getappdata(resultfig,'templateused');
if ~strcmp(templateused,togglestates.template) % reload image(s)
    clear V
    switch togglestates.template
        case 'MNI-Template'
            V{1}=spm_vol([options.earoot,'templates',filesep,'mni_hires.nii']);
            
        case 'Patient Post-Op'
            
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
xyzv= V{1}.mat \ togglestates.xyzmm;
xyzv=round(xyzv(1:3)); % now in voxel coordinates.

%colormap gray



if togglestates.xyztoggles(1)
    %try
    %slice=(double(squeeze(nii{1}.img(xyzv(1),:,:))));
    
    usesag=(length(V)>2)*2; % check if explicit saggital volume is available
    
    
    [xx,yy,zz]=meshgrid(xyzv(1),1:0.5:V{1+usesag}.dim(2),1:0.5:V{1+usesag}.dim(3));
    slice=spm_sample_vol(V{1+usesag},xx,yy,zz,4);
    slice=ea_invert(slice,inverted);
    maxv=max(slice(:));
    minv=abs(min(slice(:)));
    
    slice=slice+minv; % 0 smallest number.
    slice=slice/maxv*255; % 255 highest number.
    imin=repmat(uint8((((slice)))),[1,1,4]);
    imin(:,:,4)=uint8(togglestates.xyztransparencies(1));
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
    
    %surface('XData',[max(ulp(:,1)),min(ulp(:,1))],'YData',[min(ulp(:,2)),max(ulp(:,2))],'ZData',[min(ulp(:,3)),max(ulp(:,3))],'CData',imin(:,:,1:3));
    %catch
    %    disp('Z-Volume cut out of bounds.');
    %end
    
end

if togglestates.xyztoggles(2)
    
    % check whether second nii is being used:
    usecor=length(V)>1; % check if explicit coronar volume is available
    
    [xx,yy,zz]=meshgrid(1:0.5:V{1+usecor}.dim(1),xyzv(2),1:0.5:V{1+usecor}.dim(3));
    slice=flipud(squeeze((spm_sample_vol(V{1+usecor},squeeze(xx),squeeze(yy),squeeze(zz),2))));
    slice=ea_invert(slice,inverted);
    
    %slice=flipud(squeeze(double(nii{1+usecor}.img(:,xyzv(2),:))));
    maxv=max(slice(:));
    minv=abs(min(slice(:)));
    
    slice=slice+minv; % 0 smallest number.
    slice=slice/maxv*255; % 255 highest number.
    imin=repmat(uint8((((slice)))),[1,1,4]);
    
    
    
    %try
    %imin=repmat(uint8((squeeze((nii{1+usecor}.img(:,xyzv(2),:)+minv)*255)/(minv+maxv))),[1,1,4]);
    imin(:,:,4)=uint8(togglestates.xyztransparencies(2));
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
    
    %ysliceplot=imsurf(imin,ulp,[0,1,0],[-1,0,0],scale);
    %catch
    %    disp('Y-Volume cut out of bounds.');
    %end
end

if togglestates.xyztoggles(3)
    %try
    %imsurf(repmat(uint8((squeeze(fliplr(nii.img(:,:,graphopts.volxyz(3)))'*255)/max(nii.img(:)))),[1,1,2]),[size(nii.img,1)+0.5,size(nii.img,2)+0.5,graphopts.volxyz(3)],[0,0,1],[-1,0,0],1)
    %original: imin=repmat(uint8((squeeze(fliplr(nii.img(:,:,xyzv(1)))'*255)/max(nii.img(:)))),[1,1,4]);
      [xx,yy,zz]=meshgrid(1:0.2:V{1}.dim(1),1:0.2:V{1}.dim(2),xyzv(3));
    slice=fliplr(spm_sample_vol(V{1},xx,yy,zz,4))';
    slice=ea_invert(slice,inverted);
    %slice=flipud(squeeze(double((nii{1}.img(:,:,xyzv(3))))));
    maxv=max(slice(:));
    minv=abs(min(slice(:)));
    
    slice=slice+minv; % 0 smallest number.
    slice=slice/maxv*255; % 255 highest number.
    imin=repmat(uint8(((slice))),[1,1,4]);
    

    %imin=repmat(uint8((squeeze((flipud(nii{1}.img(:,:,xyzv(3)))+minv)*255)/(minv+maxv))),[1,1,4]);
    imin(:,:,4)=uint8(togglestates.xyztransparencies(3)); % transparency
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
    %xsliceplot=imsurf(imin,bb,[0,0,1],[-1,0,0],scale);
    
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





