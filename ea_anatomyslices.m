function ea_anatomyslices(resultfig,togglestates,options,controlhandles)
% input xyz in mm coordinates.
% this function plots an anatomical slice to the 3D viewer of lead-dbs.
% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ~exist('controlhandles', 'var')
    controlhandles = [];
end

set(0, 'currentfigure', resultfig);  % for figures
atlases=getappdata(resultfig,'atlases');
togglestates.xyztransparencies=double(togglestates.xyztransparencies/100);

xsliceplot=getappdata(resultfig,'xsliceplot');
ysliceplot=getappdata(resultfig,'ysliceplot');
zsliceplot=getappdata(resultfig,'zsliceplot');

%% Parse togglestates
if ~isempty(xsliceplot) && ~togglestates.refreshcuts
    switch togglestates.cutview
        case 'xcut'
            set(xsliceplot,'Visible','on')
            set(ysliceplot,'Visible','off')
            set(zsliceplot,'Visible','off')
            ea_settransparency(resultfig,togglestates)
            setappdata(resultfig,'xsliceplot',xsliceplot);
            setappdata(resultfig,'ysliceplot',ysliceplot);
            setappdata(resultfig,'zsliceplot',zsliceplot);
        case 'ycut'
            set(xsliceplot,'Visible','off')
            set(ysliceplot,'Visible','on')
            set(zsliceplot,'Visible','off')
            ea_settransparency(resultfig,togglestates)
            setappdata(resultfig,'xsliceplot',xsliceplot);
            setappdata(resultfig,'ysliceplot',ysliceplot);
            setappdata(resultfig,'zsliceplot',zsliceplot);
        case 'zcut'
            set(xsliceplot,'Visible','off')
            set(ysliceplot,'Visible','off')
            set(zsliceplot,'Visible','on')
            ea_settransparency(resultfig,togglestates)
            setappdata(resultfig,'xsliceplot',xsliceplot);
            setappdata(resultfig,'ysliceplot',ysliceplot);
            setappdata(resultfig,'zsliceplot',zsliceplot);
        case '3d'
            set(xsliceplot,'Visible','on')
            set(ysliceplot,'Visible','on')
            set(zsliceplot,'Visible','on')
            ea_settransparency(resultfig,togglestates)
            setappdata(resultfig,'xsliceplot',xsliceplot);
            setappdata(resultfig,'ysliceplot',ysliceplot);
            setappdata(resultfig,'zsliceplot',zsliceplot);
    end
end

if ~isempty(xsliceplot) && ~isequal(togglestates.xyztoggles,[1 1 1]) && strcmp(togglestates.cutview,'3d')
    if togglestates.xyztoggles(1)
        set(xsliceplot,'Visible','on');
    else
        set(xsliceplot,'Visible','off');
    end
    if togglestates.xyztoggles(2)
        set(ysliceplot,'Visible','on');
    else
        set(ysliceplot,'Visible','off');
    end
    if togglestates.xyztoggles(3)
        set(zsliceplot,'Visible','on');
    else
        set(zsliceplot,'Visible','off');
    end    
    ea_settransparency(resultfig,togglestates)
    setappdata(resultfig,'xsliceplot',xsliceplot);
    setappdata(resultfig,'ysliceplot',ysliceplot);
    setappdata(resultfig,'zsliceplot',zsliceplot);
end

if (~togglestates.refreshcuts) && (~togglestates.refreshview)
    return
end

%% Render slices
V=getappdata(resultfig,'V');
inverted=getappdata(resultfig,'inverted');
if isempty(inverted)
    inverted=0;
end
options.d2.writeatlases=1;

if ~isfield(options,'native')
    options.native=0;
end

if togglestates.refreshcuts % reload image(s)
    clear V
    if strcmp(togglestates.template, 'Choose...')
        togglestates.template = togglestates.customfile;
    end
    [V1, V2, V3] = ea_assignbackdrop(togglestates.template,options,'Patient',options.native);
    if ~isfield(V1,'img') % image supplied
        V{1} = nifti(V1.fname);
        V{2} = nifti(V2.fname);
        V{3} = nifti(V3.fname);
    else
        V{1}=V1; V{2}=V2; V{3}=V3;        
    end
    setappdata(resultfig,'templateused',togglestates.template); % refresh used template.
end

if ~inverted==togglestates.tinvert
    inverted=togglestates.tinvert;    
end

if ~isempty(controlhandles)
    if get(controlhandles.slicepopup,'Value')==1
        togglestates.xyzmm=[togglestates.xyzmm,1]';

        try
            xyzv= V{1}.mat \ togglestates.xyzmm;
        catch
            keyboard
        end

        xyzv=round(xyzv(1:3)); % now in voxel coordinates.
        %keyboard
    elseif get(controlhandles.slicepopup,'Value')==2
        xyzv = togglestates.xyzmm;
        % xyzv= V{1}.mat * togglestates.xyzmm;
    end
else % direct call from script.
    xyzv = round(V{1}.mat \ [togglestates.xyzmm,1]');
end

% balance the contrast
% if togglestates.refreshcuts
% [balanced,colormap] = ea_autocontrast(double(V{1}.dat),2.5);
% end

usesag=(length(V)>2)*2; % check if explicit saggital volume is available
if isa(V{1+usesag},'nifti') % memory mapped, nifti function
    xsliceplot=slice3i(resultfig,V{1+usesag}.dat,V{1+usesag}.mat,1,xyzv(1),controlhandles);
else
    xsliceplot=slice3i(resultfig,V{1+usesag}.img,V{1+usesag}.mat,1,xyzv(1),controlhandles);
end

usecor=length(V)>1; % check if explicit coronal volume is available
if isa(V{1+usecor},'nifti') % memory mapped, nifti function
    ysliceplot=slice3i(resultfig,V{1+usecor}.dat,V{1+usecor}.mat,2,xyzv(2),controlhandles);
else
    ysliceplot=slice3i(resultfig,V{1+usecor}.img,V{1+usecor}.mat,2,xyzv(2),controlhandles);
end

if isa(V{1},'nifti') % memory mapped, nifti function
    zsliceplot=slice3i(resultfig,V{1}.dat,V{1}.mat,3,xyzv(3),controlhandles);
else
    zsliceplot=slice3i(resultfig,V{1}.img,V{1}.mat,3,xyzv(3),controlhandles);
end

%colormap(cmap);

% Set template popup and slice toggle properly, useful when scripting
if isempty(controlhandles)
    awin = getappdata(resultfig, 'awin');
    if isvalid(awin) % check if not deleted (closed)
        handle = guidata(awin);
        [~, idx] = ismember(togglestates.template, handle.templatepopup.String);
        if idx
            handle.templatepopup.Value = idx;
        else
            ea_cprintf('CmdWinWarnings', 'Template not found in the list!\n');
        end
        handle.xtoggle.Value = togglestates.xyztoggles(1);
        handle.ytoggle.Value = togglestates.xyztoggles(2);
        handle.ztoggle.Value = togglestates.xyztoggles(3);
    end
end

% store data in figure
setappdata(resultfig,'xsliceplot',xsliceplot);
setappdata(resultfig,'ysliceplot',ysliceplot);
setappdata(resultfig,'zsliceplot',zsliceplot);
%ea_settransparency(resultfig,togglestates)
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

