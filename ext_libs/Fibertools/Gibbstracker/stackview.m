function stackview(imgs,showfun)

if nargin == 1,
    showfun = @imshow;
end;


if iscell(imgs),
    if length(imgs) > 3,
        display('maximum three volumes');
        return;
    end;

    for k = 1:length(imgs),
        if ~isempty(imgs{k}),
            sz = size(imgs{k});
            break;
        end;
    end;

    img = zeros([sz 3]);
    for k = 1:length(imgs),
        if isempty(imgs{k}),
            img(:,:,:,k) = zeros(sz);
        else
            img(:,:,:,k) = imgs{k};
        end;
    end;
else
    img = imgs;
end;






%f = figure('Position',[360,500,750,585]);    
f = gcf;
clf;

h = uicontrol('Style','slider','String','Surf','Position',[0,0,300,25],'Tag','theslider','Value',0.5,'Callback',{@resolutionslider_Callback,gcbo,[],[],img,showfun});
botPanel = uibuttongroup('Units','pixels','Position',[400 0  200 30],'SelectionChangeFcn',{@radiob_Callback,gcbo,img},'Tag','butgroup');
hX = uicontrol('Style','radiobutton','String','X','Position',[10,1,40,25],'Tag','1','Parent',botPanel);
hY = uicontrol('Style','radiobutton','String','Y','Position',[60,1,40,25],'Tag','2','Parent',botPanel);
hZ = uicontrol('Style','radiobutton','String','Z','Position',[110,1,40,25],'Tag','3','Parent',botPanel);
coord = uicontrol('Style','text','String','x:1','Tag','coordinfo','Position',[300,0,100,25]);
set(h,'SliderStep',[1/size(img,1) 0.5]);
showslice(img,showfun);


% --------------------------------------------------------------------
function varargout = resolutionslider_Callback(h, eventdata, handles, varargin)
  
    img = varargin{3};
    showslice(img,varargin{4});
    
% --------------------------------------------------------------------
function varargout = radiob_Callback(h, eventdata, handles, varargin)
    img = varargin{1};
     
    radiob = get(h,'SelectedObject');
    selection = str2num(get(radiob,'Tag'));    
    
    slider = findobj(gcf,'Tag','theslider');
    set(slider,'SliderStep',[1/size(img,selection) 0.5]);
    cb = get(slider,'Callback');
    
    showslice(img,cb{6});
    
function showslice(img,showfun)

    slider = findobj(gcf,'Tag','theslider');
    pos = get(slider,'value');
    botpanel = findobj(gcf,'Tag','butgroup');
    radiob = get(botpanel,'SelectedObject');
    selection = str2num(get(radiob,'Tag'));    
    idxpos = round(pos*(size(img,selection)-1))+1;
    cinfo = findobj(gcf,'Tag','coordinfo');
    
    switch selection,
        case 3,
            showfun(squeeze(img(:,:,idxpos,:)));
            set(cinfo,'String',sprintf('%i/%i',idxpos,size(img,3)));    
        case 2,
            showfun(flipdim(permute(squeeze(img(:,idxpos,:,:)),[2 1 3]),1));
            set(cinfo,'String',sprintf('%i/%i',idxpos,size(img,2)));    
        case 1,
            showfun(flipdim(permute(squeeze(img(idxpos,:,:,:)),[2 1 3]),1)); 
            set(cinfo,'String',sprintf('%i/%i',idxpos,size(img,1)));    
    end;
          
            
            
            
            
            