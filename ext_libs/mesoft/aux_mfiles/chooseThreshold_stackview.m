function threshold = chooseThreshold_stackview(img1,img2)

threshold = [];

showfun = @imshow;


if any(size(img1) ~= size(img2)) && any((2*size(img1)) ~= size(img2)),        
    errordlg('Data and Mask dimensions not consistent');
    return;
end;

if all((2*size(img1)) == size(img2))
    img.back = zeros(size(img1)*2);
    img.back(1:2:end,1:2:end,1:2:end) = img1;
    img.back(2:2:end,1:2:end,1:2:end) = img1;
    img.back(1:2:end,2:2:end,1:2:end) = img1;
    img.back(1:2:end,1:2:end,2:2:end) = img1;
    img.back(2:2:end,2:2:end,1:2:end) = img1;
    img.back(2:2:end,1:2:end,2:2:end) = img1;
    img.back(1:2:end,2:2:end,2:2:end) = img1;
    img.back(2:2:end,2:2:end,2:2:end) = img1;
else,
    img.back = img1;
end;
img.mask = img2;



maxm = max(img2(:))+eps;
minm = min(img2(:))-eps;

default = (maxm+minm)/2;

%f = figure('Position',[360,500,750,585]);    
f = figure;
set(f,'userdata',nan);
set(f,'toolbar','none','menubar','none','name','Choose a Threshold','numbertitle','off','dockcontrols','off','windowstyle','modal','closerequestfcn',{@okbut_Callback,gcbo,[],[],f,0})
clf;

h = uicontrol('Style','slider','String','Surf','Position',[0,0,300,25],'Tag','theslider','Value',0.5,'Callback',{@resolutionslider_Callback,gcbo,[],[],img,showfun});
botPanel = uibuttongroup('Units','pixels','Position',[400 0  200 30],'SelectionChangeFcn',{@radiob_Callback,gcbo,img},'Tag','butgroup');
hX = uicontrol('Style','radiobutton','String','X','Position',[10,1,40,25],'Tag','1','Parent',botPanel);
hY = uicontrol('Style','radiobutton','String','Y','Position',[60,1,40,25],'Tag','2','Parent',botPanel);
hZ = uicontrol('Style','radiobutton','String','Z','Position',[110,1,40,25],'Tag','3','Parent',botPanel);

coord = uicontrol('Style','text','String','x:1','Tag','coordinfo','Position',[300,0,100,25]);

uicontrol('Style','text','String','Choose Threshold','Position',[0,240,90,35]);
uicontrol('Style','edit','String',default,'Position',[0,210,90,25],'tag','editthreshold','Callback',{@editthres_Callback,gcbo,[],[],img,showfun});
uicontrol('Style','pushbutton','String','OK','Position',[0,140,30,20],'Tag','ok_but','Callback',{@okbut_Callback,gcbo,[],[],f,1});
uicontrol('Style','pushbutton','String','Cancel','Position',[40,140,50,20],'Tag','cancel_but','Callback',{@okbut_Callback,gcbo,[],[],f,0});

ht = uicontrol('Style','slider','String','Surf','Position',[0,170,90,25],'Tag','theslider_threshold','Value',default,'Callback',{@thresholdslider_Callback,gcbo,[],[],img,showfun},'max',maxm,'min',minm, ...
        'SliderStep',[1/1000 1/10]);

set(h,'SliderStep',[1/size(img1,1) 0.5]);
showslice(img,showfun);
set(gca,'position',[0.23 0.11 0.675 0.815])
waitfor(f,'userdata');
threshold = get(f,'userdata');
delete(f);


 
function varargout = okbut_Callback(h, eventdata, handles, varargin)
  
    h = varargin{3};
    ok = varargin{4};
    if ok 
       hs = findobj('tag','theslider_threshold');
       val = get(hs,'Value');
       set(h,'userdata',val);
    else
       set(h,'userdata',[]);
    end;
   
 
% --------------------------------------------------------------------
function varargout = resolutionslider_Callback(h, eventdata, handles, varargin)
  
    img = varargin{3};
    showslice(img,varargin{4});
   
function varargout = editthres_Callback(h, eventdata, handles, varargin)
 
    hs = findobj('tag','theslider_threshold');
    set(hs,'Value',str2num(get(h,'String')));

    img = varargin{3};
    showslice(img,varargin{4});
     
    
function varargout = thresholdslider_Callback(h, eventdata, handles, varargin)
  
    ht = findobj('tag','editthreshold');
    set(ht,'String',get(h,'Value'));
    
    img = varargin{3};
    showslice(img,varargin{4});
  
     
% --------------------------------------------------------------------
function varargout = radiob_Callback(h, eventdata, handles, varargin)
    img = varargin{1};
     
    radiob = get(h,'SelectedObject');
    selection = str2num(get(radiob,'Tag'));    
    
    slider = findobj(gcf,'Tag','theslider');
    set(slider,'SliderStep',[1/size(img.back,selection) 0.5]);
    cb = get(slider,'Callback');
    
    showslice(img,cb{6});
    
function showslice(img,showfun)

    slider = findobj(gcf,'Tag','theslider');
    slider2 = findobj(gcf,'Tag','theslider_threshold');
    pos = get(slider,'value');
    threshold = get(slider2,'value');
    botpanel = findobj(gcf,'Tag','butgroup');
    radiob = get(botpanel,'SelectedObject');
    selection = str2num(get(radiob,'Tag'));    
    idxpos = round(pos*(size(img.back,selection)-1))+1;
    cinfo = findobj(gcf,'Tag','coordinfo');
    
    switch selection,
        case 3,
            toshow(:,:,1) = squeeze(img.back(:,:,idxpos,:));
            toshow(:,:,2) = squeeze(img.mask(:,:,idxpos,:)>threshold);
            toshow(:,:,3) = 0;
            set(cinfo,'String',sprintf('%i/%i',idxpos,size(img,3)));    
        case 2,
           toshow(:,:,1) = flipdim(permute(squeeze(img.back(:,idxpos,:,:)),[2 1 3]),1);
           toshow(:,:,2) = flipdim(permute(squeeze(img.mask(:,idxpos,:,:)),[2 1 3]),1)>threshold;
           toshow(:,:,3) = 0;            
            set(cinfo,'String',sprintf('%i/%i',idxpos,size(img,2)));    
        case 1,
            toshow(:,:,1) = flipdim(permute(squeeze(img.back(idxpos,:,:,:)),[2 1 3]),1);
            toshow(:,:,2) = flipdim(permute(squeeze(img.mask(idxpos,:,:,:)),[2 1 3]),1)>threshold;
            toshow(:,:,3) = 0;
            set(cinfo,'String',sprintf('%i/%i',idxpos,size(img,1)));    
    end;
    toshow(:,:,1) = 2*toshow(:,:,1) / (eps+max(max(toshow(:,:,1))));
    toshow(:,:,2) = 0.8*toshow(:,:,2) / (eps+max(max(toshow(:,:,2))));
    showfun(toshow);
          
            
            
            
            
            