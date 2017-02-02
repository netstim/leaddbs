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
ysliceplot=getappdata(resultfig,'ysliceplot');
zsliceplot=getappdata(resultfig,'zsliceplot');

%% Parse togglestates
if ~isempty(xsliceplot)
switch togglestates.cutview
    case 'xcut'
        set(xsliceplot,'Visible','on')
        set(ysliceplot,'Visible','off')
        set(zsliceplot,'Visible','off')
        ea_settransparency(resultfig,togglestates)
        setappdata(resultfig,'xsliceplot',xsliceplot);
        setappdata(resultfig,'ysliceplot',ysliceplot);
        setappdata(resultfig,'zsliceplot',zsliceplot);
        return
    case 'ycut'
        set(xsliceplot,'Visible','off')
        set(ysliceplot,'Visible','on')
        set(zsliceplot,'Visible','off')
        ea_settransparency(resultfig,togglestates)
        setappdata(resultfig,'xsliceplot',xsliceplot);
        setappdata(resultfig,'ysliceplot',ysliceplot);
        setappdata(resultfig,'zsliceplot',zsliceplot);
        return
    case 'zcut'
        set(xsliceplot,'Visible','off')
        set(ysliceplot,'Visible','off')
        set(zsliceplot,'Visible','on')
        ea_settransparency(resultfig,togglestates)
        setappdata(resultfig,'xsliceplot',xsliceplot);
        setappdata(resultfig,'ysliceplot',ysliceplot);
        setappdata(resultfig,'zsliceplot',zsliceplot);
        return
    case '3d'
        set(xsliceplot,'Visible','on')
        set(ysliceplot,'Visible','on')
        set(zsliceplot,'Visible','on')
        ea_settransparency(resultfig,togglestates)
        setappdata(resultfig,'xsliceplot',xsliceplot);
        setappdata(resultfig,'ysliceplot',ysliceplot);
        setappdata(resultfig,'zsliceplot',zsliceplot);
        return
end
end

%% Render slices
V=getappdata(resultfig,'V');
inverted=getappdata(resultfig,'inverted');
if isempty(inverted)
    inverted=0;
end
options.d2.writeatlases=1;
templateused=getappdata(resultfig,'templateused');


mcr=ea_checkmacaque(options);


if ~isfield(options,'native')
    options.native=0;
end

if ~strcmp(templateused,togglestates.template) || isempty(V) % reload image(s)
    clear V
    [V1,V2,V3]=ea_assignbackdrop(togglestates.template,options,'Patient',options.native);

    V{1}=V1; V{2}=V2; V{3}=V3;
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

% balance the contrast
[balanced, cmap] = ea_autocontrast(double(V{1}.private.dat),2.5);

if togglestates.xyztoggles(1)

    
    usesag=(length(V)>2)*2; % check if explicit saggital volume is available
    if inverted
        [~,slice]=ea_writeplanes(options, togglestates.xyzmm(1),3,V{1+usesag},'off', 0,atlases);
        
        if V{1}.mat(11)>0
            slice=flip(permute(double(slice),[2,1,3]),2);
        else
            slice=permute(double(slice),[2,1,3]);
        end
        
        if V{1}.mat(6)<0
            slice=flip(slice,1);
        end

    end
    
    xsliceplot=slice3i(balanced,V{1+usesag}.mat,1,xyzv(1));
    
end

if togglestates.xyztoggles(2)
    
    % check whether second nii is being used:
    usecor=length(V)>1; % check if explicit coronal volume is available
    if inverted
        [~,slice]=ea_writeplanes(options, togglestates.xyzmm(2),2,V{1+usecor},'off', 0,atlases);
            slice=permute(double(slice),[2,1,3]);
            slice=flip(flip(slice,1),2);
    else
    end
    
    ysliceplot=slice3i(balanced,V{1+usecor}.mat,2,xyzv(2));
    
end

if togglestates.xyztoggles(3)
    
    if inverted
        [~,slice]=ea_writeplanes(options, togglestates.xyzmm(3),1,V{1},'off', 0,atlases);
        if V{1}.mat(1)>0
        slice=flip(permute(double(slice),[2,1,3]),2);
        else
        slice=permute(double(slice),[2,1,3]);
        end
        
        if V{1}.mat(6)<0
            slice=flip(slice,1);
        end
        
    else
    end
     
    zsliceplot=slice3i(balanced,V{1}.mat,3,xyzv(3));    
    
end

colormap(cmap);

% store data in figure
setappdata(resultfig,'xsliceplot',xsliceplot);
setappdata(resultfig,'ysliceplot',ysliceplot);
setappdata(resultfig,'zsliceplot',zsliceplot);
ea_settransparency(resultfig,togglestates)
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

    
