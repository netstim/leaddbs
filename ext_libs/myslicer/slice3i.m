function h = slice3i(resultfig, vol, I2X, slicedim, sliceidx, controlhandles)
% Display a slice from a volume in 3-D
% h = slice3(vol, I2X, slicedim, sliceidx, handle)
%
% Vol is either a scalar or RGB volume, e.g. N x M x K or N x M x K x 3.
% I2X is a transformation matrix from volume coordinates to xyz 3-D world
% coordinates, similar to the transform used by image3.
%
%     [x y z 1]' = I2X * [i j k 1]'
%
% slicedim  is 1, 2 or 3 and corresponds to i, j and k.
% sliceidx  the index of the slice along slicedim
% handle    an optional handle to a previous slice to reuse
% h         the handle to the created slice.
%
% Example:
%
% load mri;
% T = [1 0 0 0;0 1 0 0;0 0 2.5 0];
% h1 = slice3i(squeeze(D),T,1,64);
% h2 = slice3i(squeeze(D),T,2,64);
% h3 = slice3i(squeeze(D),T,3,14);
% colormap gray(88);
% view(30,45); axis equal; axis vis3d;
%
% SEE ALSO: image3, imagesc3
%
% Author: Anders Brun, anders@cb.uu.se (2009)
% ==========================================

switch slicedim
    case 1
        h=getappdata(resultfig,'xsliceplot');
    case 2

        h=getappdata(resultfig,'ysliceplot');
    case 3

        h=getappdata(resultfig,'zsliceplot');

end
if ~exist('controlhandles','var')
    controlhandles='';
end

h = update_slice(vol, I2X, slicedim, sliceidx, h, controlhandles, resultfig);

% set up gui
gui.handle = h;
gui.vol = vol;
gui.I2X = I2X;
gui.slicedim = slicedim;
gui.sliceidx = sliceidx;
set(gui.handle,'ButtonDownFcn',{@startmovit, controlhandles, resultfig});
set(h,'UserData',gui);


function startmovit(src,evnt,controlhandles,resultfig)
% Unpack gui object
gui = get(src,'UserData');
%turn off mouse pointer
%set(gcf,'PointerShapeCData',nan(16,16));
%set(gcf,'Pointer','custom');

thisfig = gcbf();
gui.startray = get(gca,'CurrentPoint');
gui.startidx = gui.sliceidx;
% Store gui object
set(src,'UserData',gui);
set(thisfig,'WindowButtonMotionFcn',{@movit,controlhandles,resultfig});
set(thisfig,'WindowButtonUpFcn',@stopmovit);
set(thisfig,'UserData',src);


function movit(src,evnt,controlhandles,resultfig)
% Unpack gui object
gui = get(get(gcf,'UserData'),'UserData');
% Some safetymeasures
try
    if isequal(gui.startray,[])
        return
    end
end

% Do "smart" positioning of the markers...
nowray = get(gca,'CurrentPoint');

% Project rays on slice-axis
s = gui.I2X(1:3,gui.slicedim);
a = gui.startray(1,:)';
b = gui.startray(2,:)';
alphabeta = pinv([s'*s, -s'*(b-a);(b-a)'*s, -(b-a)'*(b-a)])*[s'*a, (b-a)'*a]';
pstart = alphabeta(1)*s;
alphastart = alphabeta(1);
a = nowray(1,:)';
b = nowray(2,:)';
alphabeta = pinv([s'*s, -s'*(b-a);(b-a)'*s, -(b-a)'*(b-a)])*[s'*a, (b-a)'*a]';
pnow = alphabeta(1)*s;
alphanow = alphabeta(1);
slicediff = alphanow-alphastart;

gui.sliceidx = gui.startidx+slicediff;
gui.sliceidx = min(max(1,gui.sliceidx),size(gui.vol,gui.slicedim));
update_slice(gui.vol, gui.I2X, gui.slicedim, gui.sliceidx, gui.handle, controlhandles, resultfig);
drawnow;

% Store gui object
set(get(gcf,'UserData'),'UserData',gui);


function stopmovit(src,evnt)
thisfig = gcbf();
%set(gcf,'Pointer','arrow');
set(thisfig,'WindowButtonUpFcn','');
set(thisfig,'WindowButtonMotionFcn','');
gui = get(get(gcf,'UserData'),'UserData');
gui.startray = [];
%gui.sliceidx = gui.nextsliceidx
set(get(gcf,'UserData'),'UserData',gui);
set(gcf,'UserData',[]);
drawnow;


function h = update_slice(vol, I2X, slicedim, sliceidx, handle, controlhandles, resultfig)

if ndims(vol) == 3         %Scalar mode
elseif ndims(vol) == 4     %RGB mode
else
	error('Only scalar and RGB images supported')
end

% Create the slice
if slicedim == 3 % k
    ij2xyz = I2X(:,[1 2]);
    ij2xyz(:,3) = I2X*[0 0 sliceidx 1]';
    if round(sliceidx)<1
    	sliceidx=1;
    elseif round(sliceidx)>size(vol,3)
    	sliceidx=size(vol,3);
    end
    sliceim = squeeze(vol(:,:,round(sliceidx),:));
elseif slicedim == 2 % j
  ij2xyz = I2X(:,[1 3]);
  ij2xyz(:,3) = I2X*[0 sliceidx 0 1]';
  if round(sliceidx)<1
      sliceidx=1;
  elseif round(sliceidx)>size(vol,2)
      sliceidx=size(vol,2);
  end
  sliceim = squeeze(vol(:,round(sliceidx),:,:));
elseif slicedim == 1 % i
  ij2xyz = I2X(:,[2 3]);
  ij2xyz(:,3) = I2X*[sliceidx 0 0 1]';
  if round(sliceidx)<1
      sliceidx=1;
  elseif round(sliceidx)>size(vol,1)
      sliceidx=size(vol,1);
  end
  sliceim = squeeze(vol(round(sliceidx),:,:,:));
else
    error('Slicedim should be 1, 2 or 3')
end
c=1; o=0;
sc=getappdata(resultfig,'slidecontrast');

if ~isempty(sc) % add contrast from user GUI
    c=c+sc.c;
    o=o+sc.o;
end

% if ~isempty(strfind(controlhandles.templatepopup.String{controlhandles.templatepopup.Value},'BigBrain'))
%     sliceim=ea_contrast(sliceim,c,o)*64;
% end

if length(size(sliceim))==2
    sliceim=ea_contrast(sliceim,c,o)*length(gray);
else
    sliceim=uint8(ea_contrast(single(sliceim),c,o)*255);
end

resdivs=1; % could increase to 2 but would render a bit slow.
if size(sliceim,3)==1
    sliceim=interp2(sliceim,resdivs);
else
    resdivs=0;
end

ij2xyz(:,1:2)=ij2xyz(:,1:2)/(2^resdivs);
try
    ea_update_anatomycontrol(sliceidx,slicedim,I2X,controlhandles);
end

if isempty(handle)
	h = image3(sliceim,ij2xyz);
else
	h = image3(sliceim,ij2xyz,handle);
end


function ea_update_anatomycontrol(sliceidx,slicedim,mat,controlhandles)

slicecoord=[ones(4,1)];
slicecoord(slicedim)=sliceidx;
slicecoord=mat*slicecoord;
switch slicedim
    case 1
        slicehdl='xval';
    case 2
        slicehdl='yval';
    case 3
        slicehdl='zval';
end

if get(controlhandles.slicepopup,'Value')==1
    set(controlhandles.(slicehdl),'String',sprintf('%1.1f',slicecoord(slicedim)));
elseif get(controlhandles.slicepopup,'Value')==2
    set(controlhandles.(slicehdl),'String',sprintf('%1.0f',sliceidx));
end
