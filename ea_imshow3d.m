function ea_imshow3d(nii)
%
% This function displays a NifTi image. 
%
% SYNTAX: 
%       imshow3d('anat_t1.nii')
%       imshow3d('Users/user/Images/anat_t1.nii')
%       imshow3d(anat_t1)
%       imshow3d(anat_t1.img)
%
% Example1:
%       nii = load_nii('t1.nii')
%       imshow3d(nii)
%
% Example2:
%       anat_t1 = load_nii('Users/user/Images/anat_t1.nii')
%       imshow3d(anat_t1.img)
%
% INPUT:
%   nii: 
%        case 1) filename: full filename to .nii image or name of file in 
%                current diretory
%        case 2) nii: structure loaded into matlab workspace from load_nii.m 
%        case 3) img: NxMx3 image in workspace, e.g. nii.img from load_nii.m
%
%
% Part of this file is copied and modified from: 
% imshow3Dfull by Maysam Shahedi (mshahedi@gmail.com), September 22, 2016
% Available for Download at:
% https://www.mathworks.com/matlabcentral/fileexchange/47463-imshow3dfull--3d-imshow-in-3-views-
% __________________________________________________________________________________
% Copyright (C) 2023 Brigham and Women's Hospital, Boston, MA
% Ari Kappel, MD
%% 
%%%%%%%%% USER SETS PREFERENCES %%%%%%%%%%%%%
% Set Prefs
WLPref = 1024*32; % Default WLPref is 1024
ZoomCoef = 0.05;  % Default ZoomCoef is 0.05
toggleZoom = 'on';
toggleAll=  'on'; % All except slider and Finetune

toggleSlider = 'off';
toggleFineTune = 'off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
if ischar(nii)
    filename = nii;
    % check if file exists
    if ~exist(filename,'file')
        sprintf('%sr', ['File does not exist: Cannot find ' filename])
        [filename,path] = uigetfile();
        filename = fullfile(path,filename);
    end
    try
        nii = load_nii(filename);
    catch
        nii = load_untouch_nii(filename);
    end
    nii = nii.img;
    
elseif isstruct(nii)
    if isfield(nii,'img')
        nii = nii.img;
    else
        error('Error in input file')
    end
end

% Check for 3d image
% dim = size(squeeze(nii));
if length(size(squeeze(nii)))<3
    error('Error: Please choose a 3-dimensional image')
end
    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isequal(toggleFineTune,'on')
    lengthUI = 530;
else
    toggleFineTune = 'off';
    lengthUI = 530-80;
end
if strcmp(toggleAll,'off')
    toggleSlider = 'off';
    toggleFineTune = 'off';
    lengthUI = 530-80;
end

try
    figtit = filename;
catch
    figtit = '';
end

%% Open figure;
isp=figure('color','k','Name',figtit,'NumberTitle','off');
% isp=figure('color','k','Name',figtit,'NumberTitle','off','MenuBar','none','DockControls','off','ToolBar','none');
% set(gcf,'MenuBar','figure','ToolBar','auto')

sno = size(nii);
sno = size(nii);  % image size
sno_a = sno(3);  % number of axial slices
S_a = round(sno_a/2);
sno_s = sno(2);  % number of sagittal slices
S_s = round(sno_s/2);
sno_c = sno(1);  % number of coronal slices
S_c = round(sno_c/2);
S = S_a;
sno = sno_a;

global InitialCoord;

MinV = 0;
MaxV = nanmax(nii(:)); %dbs_nanmax
LevV = (double( MaxV) + double(MinV)) / 2;
Win = double(MaxV) - double(MinV);
WLAdjCoe = (Win + 1)/WLPref;
FineTuneC = [1 1/16];    % Regular/Fine-tune mode coefficients


if isa(nii,'uint8')
    MaxV = uint8(Inf);
    MinV = uint8(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/WLPref;
elseif isa(nii,'uint16')
    MaxV = uint16(Inf);
    MinV = uint16(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/WLPref;
elseif isa(nii,'uint32')
    MaxV = uint32(Inf);
    MinV = uint32(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/WLPref;
elseif isa(nii,'uint64')
    MaxV = uint64(Inf);
    MinV = uint64(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/WLPref;
elseif isa(nii,'int8')
    MaxV = int8(Inf);
    MinV = int8(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/WLPref;
elseif isa(nii,'int16')
    MaxV = int16(Inf);
    MinV = int16(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/WLPref;
elseif isa(nii,'int32')
    MaxV = int32(Inf);
    MinV = int32(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/WLPref;
elseif isa(nii,'int64')
    MaxV = int64(Inf);
    MinV = int64(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/WLPref;
elseif isa(nii,'logical')
    MaxV = 0;
    MinV = 1;
    LevV =0.5;
    Win = 1;
    WLAdjCoe = 0.1;
else
    MaxV = 0;
    MinV = 1;
    LevV =0.5;
    Win = 1;
    WLAdjCoe = 0.1;
end

ImgAx = nii;
if verLessThan('matlab', '8')
    ImgSg = flipdim(permute(nii, [3 1 2 4]),1);   % Sagittal view image
    ImgCr = flipdim(permute(nii, [3 2 1 4]),1);   % Coronal view image
else
    ImgSg = flip(permute(nii, [3 1 2 4]),1);   % Sagittal view image
    ImgCr = flip(permute(nii, [3 2 1 4]),1);   % Coronal view image
end

View = 'A';

SFntSz = 9;
LFntSz = 10;
WFntSz = 10;
VwFntSz = 10;
LVFntSz = 9;
WVFntSz = 9;
BtnSz = 10;
ChBxSz = 10;

% if (nargin < 2)
    [Rmin Rmax] = WL2R(Win, LevV);
% elseif numel(disprange) == 0
%     [Rmin Rmax] = WL2R(Win, LevV);
% else
%     LevV = (double(disprange(2)) + double(disprange(1))) / 2;
%     Win = double(disprange(2)) - double(disprange(1));
%     WLAdjCoe = (Win + 1)/WLPref;
%     [Rmin Rmax] = WL2R(Win, LevV);
% end


% Get/Set Window Handles
% set(gcf,'MenuBar','figure','ToolBar','auto')
% [xStart yStart width height]
set(gcf,'Position',[125   250   868   750]) %[125   250   1038   750]
FigPos = get(gcf,'Position');
startUI = (FigPos(3)-lengthUI)/2;
S_Pos = [startUI 38 lengthUI 20]; % [273 39 493 20]
Stxt_Pos = [startUI+lengthUI-71 FigPos(4)-29 71 14]; %[711 721 55 14]

UInow = startUI;    Wtxt_Pos    = [UInow 20 60 20]; % 255
UInow = UInow +55;  Wval_Pos    = [UInow 20 60 20]; % 55; 310
UInow = UInow +65;  Ltxt_Pos    = [UInow 20 45 20]; % 65; 375
UInow = UInow +40;  Lval_Pos    = [UInow 20 60 20]; % 40; 415
UInow = UInow +65;  Btn_Pos     = [UInow 20 80 20]; %70; 485
UInow = UInow +85; Vwtxt_Pos    = [UInow 20 35 20]; % 110; 595
UInow = UInow +40;  VAxBtn_Pos  = [UInow 20 15 20]; % 40
UInow = UInow +20;  VSgBtn_Pos  = [UInow 20 15 20]; % 20
UInow = UInow +20;  VCrBtn_Pos  = [UInow 20 15 20]; % 20
UInow = UInow +20;  VReBtn_Pos  = [UInow 20 35 20]; % 20
UInow = UInow +40;  ChBx_Pos    = [UInow 20 80 20];   % 30
% Total Length: 450+80

% Set UI handles
if sno > 1
    shand = uicontrol('Visible',toggleSlider,'Style', 'slider','Min',1,'Max',sno,'Value',S,'SliderStep',[1/(sno-1) 10/(sno-1)],'Position', S_Pos,'Callback', {@SliceSlider, nii});
    stxthand = uicontrol('Visible',toggleAll,'Style', 'text','Position', Stxt_Pos,'String',sprintf('Slice# %d / %d',S, sno), 'BackgroundColor', 'k', 'ForegroundColor','w','FontSize', SFntSz);
else
    stxthand = uicontrol('Visible',toggleAll,'Style', 'text','Position', Stxt_Pos,'String','2D image', 'BackgroundColor', 'k','ForegroundColor', 'w', 'FontSize', SFntSz);
end  
wtxthand = uicontrol('Visible',toggleAll,'Style', 'text','Position', Wtxt_Pos,'String','Window: ', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', WFntSz);
wvalhand = uicontrol('Visible',toggleAll,'Style', 'edit','Position', Wval_Pos,'String',sprintf('%6.0f',Win), 'BackgroundColor', [1 1 1], 'FontSize', WVFntSz,'Callback', @WinLevChanged);
ltxthand = uicontrol('Visible',toggleAll,'Style', 'text','Position', Ltxt_Pos,'String','Level: ', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', LFntSz);
lvalhand = uicontrol('Visible',toggleAll,'Style', 'edit','Position', Lval_Pos,'String',sprintf('%6.0f',LevV), 'BackgroundColor', [1 1 1], 'FontSize', LVFntSz,'Callback', @WinLevChanged);
Btnhand = uicontrol('Visible',toggleAll,'Style', 'pushbutton','Position', Btn_Pos,'String','Auto W/L', 'FontSize', BtnSz, 'Callback' , @AutoAdjust);
Vwtxthand = uicontrol('Visible',toggleAll,'Style', 'text','Position', Vwtxt_Pos,'String','View: ', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', LFntSz);
VAxBtnhand = uicontrol('Visible',toggleAll,'Style', 'pushbutton','Position', VAxBtn_Pos,'String','A', 'FontSize', BtnSz, 'Callback' , @AxialView);
VSgBtnhand = uicontrol('Visible',toggleAll,'Style', 'pushbutton','Position', VSgBtn_Pos,'String','S', 'FontSize', BtnSz, 'Callback' , @SagittalView);
VCrBtnhand = uicontrol('Visible',toggleAll,'Style', 'pushbutton','Position', VCrBtn_Pos,'String','C', 'FontSize', BtnSz, 'Callback' , @CoronalView);
VReBtnhand = uicontrol('Visible',toggleAll,'Style', 'pushbutton','Position', VReBtn_Pos,'String','Reset', 'FontSize', BtnSz, 'Callback' , @ResetView);
ChBxhand = uicontrol('Visible',toggleAll,'Visible',toggleFineTune,'Style', 'checkbox','Position', ChBx_Pos,'String','Fine Tune', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', ChBxSz);

% Set MainImage
MainImage = 1;
XImage=1:size(nii,1);
YImage=1:size(nii,2);

% SHOW IMAGE
% set(gcf,'Position',[134   232   959   742])% to use medium window size; % dbs_maximize(isp); %to maximize window
try % image toolbox
    hdl_im=imshow(squeeze(nii(XImage,YImage,S,MainImage)), [Rmin Rmax]);
catch
    hdl_im=imagesc(squeeze(nii(XImage,YImage,S,MainImage)), [Rmin Rmax]);
end;
%
axis image
AutoAdjust

origPos = get(gca,'Position');
set(get(gca,'children'),'cdata',squeeze(nii(:,:,S,:)))
set (gcf, 'WindowScrollWheelFcn', @mouseScroll);
set (gcf, 'ButtonDownFcn', @mouseClick);
set(get(gca,'Children'),'ButtonDownFcn', @mouseClick);
set(gcf,'WindowButtonUpFcn', @mouseRelease)
set(gcf,'ResizeFcn', @figureResized)

% for i = 1:size(nii,3)
% level{i} = graythresh(nii(:,:,i));
% BW{i} = imbinarize(nii(:,:,i),level{i});
% end

% RI = imref2d(size(nii));
% imshow(BW{i},RI)


% -=< Figure resize callback function >=-
    function figureResized(object, eventdata)
        FigPos = get(gcf,'Position');
        startUI = (FigPos(3)-lengthUI)/2; 
        S_Pos = [startUI 38 lengthUI 20]; %[startUI 38 (startUI+530-startUI) 20];
        Stxt_Pos = [startUI+lengthUI-71 FigPos(4)-29 71 14];
        if sno > 1
            set(shand,'Position', S_Pos);
        end
        UInow = startUI;    Wtxt_Pos    = [UInow 20 60 20]; % 255
        UInow = UInow +55;  Wval_Pos    = [UInow 20 60 20]; % 55; 310
        UInow = UInow +65;  Ltxt_Pos    = [UInow 20 45 20]; % 65; 375
        UInow = UInow +40;  Lval_Pos    = [UInow 20 60 20]; % 40; 415
        UInow = UInow +65;  Btn_Pos     = [UInow 20 80 20]; %70; 485
        UInow = UInow +85;  Vwtxt_Pos    = [UInow 20 35 20]; % 110; 595
        UInow = UInow +40;  VAxBtn_Pos  = [UInow 20 15 20]; % 40
        UInow = UInow +20;  VSgBtn_Pos  = [UInow 20 15 20]; % 20
        UInow = UInow +20;  VCrBtn_Pos  = [UInow 20 15 20]; % 20
        UInow = UInow +20;  VReBtn_Pos  = [UInow 20 35 20]; % 20
        UInow = UInow +40;  ChBx_Pos    = [UInow 20 80 20];   % 30
        % Total Length: 450
        
        set(wtxthand,'Position', Wtxt_Pos);
        set(wvalhand,'Position', Wval_Pos);
        set(ltxthand,'Position', Ltxt_Pos);
        set(lvalhand,'Position', Lval_Pos);
        set(Btnhand,'Position', Btn_Pos);
        set(ChBxhand,'Position', ChBx_Pos);
        set(Vwtxthand,'Position', Vwtxt_Pos);
        set(VAxBtnhand,'Position', VAxBtn_Pos);
        set(VSgBtnhand,'Position', VSgBtn_Pos);
        set(VCrBtnhand,'Position', VCrBtn_Pos);
        set(VReBtnhand,'Position', VReBtn_Pos);
        set(stxthand,'Position', Stxt_Pos);
        
    end

% -=< Slice slider callback function >=-
    function SliceSlider (hObj,event, Img)
        S = round(get(hObj,'Value'));
        set(get(gca,'children'),'cdata',squeeze(Img(:,:,S,:)))
        caxis([Rmin Rmax])
        if sno > 1
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            set(stxthand, 'String', '2D image');
        end
    end

% -=< Mouse scroll wheel callback function >=-
    function mouseScroll (object, eventdata)
        UPDN = eventdata.VerticalScrollCount;
        S = S - UPDN;
        if (S < 1)
            S = 1;
        elseif (S > sno)
            S = sno;
        end
        if sno > 1
            set(shand,'Value',S);
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            set(stxthand, 'String', '2D image');
        end
        set(get(gca,'children'),'cdata',squeeze(nii(:,:,S,:)))
    end

% -=< Mouse button released callback function >=-
    function mouseRelease (object,eventdata)
        set(gcf, 'WindowButtonMotionFcn', '')
    end

% -=< Mouse click callback function >=-
    function mouseClick (object, eventdata)
        MouseStat = get(gcbf, 'SelectionType');
        if strcmp(MouseStat,'normal')        %   Left CLICK
            InitialCoord = get(0,'PointerLocation');
            set(gcf, 'WindowButtonMotionFcn', @WinLevAdj);
        elseif (MouseStat(1) == 'a')   && strcmp(toggleZoom,'on')     %   RIGHT CLICK OPTION TO ZOOM
            InitialCoord = get(0,'PointerLocation');
            set(gcf, 'WindowButtonMotionFcn', @ZoomAdj);
        end
    end

% -=< Window and level mouse adjustment >=-
    function WinLevAdj(varargin)
        PosDiff = get(0,'PointerLocation') - InitialCoord;

        Win = Win + PosDiff(1) * WLAdjCoe * FineTuneC(get(ChBxhand,'Value')+1);
        LevV = LevV - PosDiff(2) * WLAdjCoe * FineTuneC(get(ChBxhand,'Value')+1);
        if (Win < 1)
            Win = 1;
        end

        [Rmin, Rmax] = WL2R(Win,LevV);
        caxis([Rmin, Rmax])
        set(lvalhand, 'String', sprintf('%6.0f',LevV));
        set(wvalhand, 'String', sprintf('%6.0f',Win));
        InitialCoord = get(0,'PointerLocation');
    end

% -=< Zoom mouse adjustment >=-
    function ZoomAdj(varargin)
            % (2) = up[+]/down[-] , right[+]/left[-]
            % hdl_im = axes('position',[0.13,0.11,0.755,0.815]);
        PosDiff = get(0,'PointerLocation') - InitialCoord;
        ImgPos = get(hdl_im.Parent,'Position');
        if ~isequal(PosDiff(2),0)
            ImgPos(1) = ImgPos(1) - ZoomCoef/PosDiff(2)*0.5;
            ImgPos(2) = ImgPos(2) - ZoomCoef/PosDiff(2)*0.5;
            ImgPos(3) = ImgPos(3) + ZoomCoef/PosDiff(2)*1;
            ImgPos(4) = ImgPos(4) + ZoomCoef/PosDiff(2)*1;
        else
            return
        end
        set(hdl_im.Parent,'Position',ImgPos)
        
        InitialCoord = get(0,'PointerLocation');
    end

% -=< Window and level text adjustment >=-
    function WinLevChanged(varargin)

        LevV = str2double(get(lvalhand, 'string'));
        Win = str2double(get(wvalhand, 'string'));
        if (Win < 1)
            Win = 1;
        end

        [Rmin, Rmax] = WL2R(Win,LevV);
        caxis([Rmin, Rmax])
    end

% -=< Window and level to range conversion >=-
    function [Rmn Rmx] = WL2R(W,L)
        Rmn = L - (W/2);
        Rmx = L + (W/2);
        if (Rmn >= Rmx)
            Rmx = Rmn + 1;
        end
    end

% -=< Window and level auto adjustment callback function >=-
    function AutoAdjust(object,eventdata)
        Win = double(max(nii(:))-min(nii(:)));
        Win (Win < 1) = 1;
        LevV = double(min(nii(:)) + (Win/2));
        [Rmin, Rmax] = WL2R(Win,LevV);
        caxis([Rmin, Rmax])
        set(lvalhand, 'String', sprintf('%6.0f',LevV));
        set(wvalhand, 'String', sprintf('%6.0f',Win));
    end

% -=< Axial view callback function >=-
    function AxialView(object,eventdata)
        if View == 'S'
            S_s = S;
        elseif View == 'C'
            S_c = S;
        end            
        View = 'A';
        
        nii = ImgAx;
        S = S_a;
        sno = sno_a;
        cla(gcf);
        
        hdl_im = imshow(squeeze(nii(:,:,S,:)), [Rmin Rmax]);
        % hdl_im=imshow(squeeze(nii(XImage,YImage,S,MainImage)), [Rmin Rmax]);
        % hdl_im = axes('position',hdl_im.Parent);
        % hdl_im = axes('position',[0.13,0.11,0.755,0.815]);
        
        if sno > 1
            shand = uicontrol('Visible',toggleSlider,'Style', 'slider','Min',1,'Max',sno,'Value',S,'SliderStep',[1/(sno-1) 10/(sno-1)],'Position', S_Pos,'Callback', {@SliceSlider, nii});
            stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String',sprintf('Slice# %d / %d',S, sno), 'BackgroundColor', 'k','ForegroundColor', 'w', 'FontSize', SFntSz);
        else
            stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String','2D image', 'BackgroundColor', 'k','ForegroundColor', 'w', 'FontSize', SFntSz);
        end
        
        caxis([Rmin Rmax])
        if sno > 1
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            set(stxthand, 'String', '2D image');
        end
        
        set(get(gca,'children'),'cdata',squeeze(nii(:,:,S,:)))
        set (gcf, 'ButtonDownFcn', @mouseClick);
        set(get(gca,'Children'),'ButtonDownFcn', @mouseClick);
    end

% -=< Sagittal view callback function >=-
    function SagittalView(object,eventdata)
        if View == 'A'
            S_a = S;
        elseif View == 'C'
            S_c = S;
        end            
        View = 'S';

        nii = ImgSg;
        S = S_s;
        sno = sno_s;
        %         hdl_im=imshow(squeeze(nii(XImage,YImage,S,MainImage)), [Rmin Rmax]);
        %         cla(hdl_im);
        %         hdl_im = axes('position',[0.13,0.11,0.755,0.815]);
        hdl_im = imshow(squeeze(nii(:,:,S,:)), [Rmin Rmax]);
        axis normal

        if sno > 1
            shand = uicontrol('Visible',toggleSlider,'Style', 'slider','Min',1,'Max',sno,'Value',S,'SliderStep',[1/(sno-1) 10/(sno-1)],'Position', S_Pos,'Callback', {@SliceSlider, nii});
            stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String',sprintf('Slice# %d / %d',S, sno), 'BackgroundColor', 'k','ForegroundColor', 'w', 'FontSize', SFntSz);
        else
            stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String','2D image', 'BackgroundColor', 'k','ForegroundColor', 'w', 'FontSize', SFntSz);
        end
        
        caxis([Rmin Rmax])
        if sno > 1
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            set(stxthand, 'String', '2D image');
        end

        set(get(gca,'children'),'cdata',squeeze(nii(:,:,S,:)))
        set (gcf, 'ButtonDownFcn', @mouseClick);
        set(get(gca,'Children'),'ButtonDownFcn', @mouseClick);

    end

% -=< Coronal view callback function >=-
    function CoronalView(object,eventdata)
        if View == 'A'
            S_a = S;
        elseif View == 'S'
            S_s = S;
        end            
        View = 'C';
        
        nii = ImgCr;
        S = S_c;
        sno = sno_c;
        %         cla(hdl_im);
        %         hdl_im = axes('position',[0.13,0.11,0.755,0.815]);
        hdl_im = imshow(squeeze(nii(:,:,S,:)), [Rmin Rmax]);
        axis normal
        
        if sno > 1
            shand = uicontrol('Visible',toggleSlider,'Style', 'slider','Min',1,'Max',sno,'Value',S,'SliderStep',[1/(sno-1) 10/(sno-1)],'Position', S_Pos,'Callback', {@SliceSlider, nii});
            stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String',sprintf('Slice# %d / %d',S, sno), 'BackgroundColor', 'k','ForegroundColor', 'w', 'FontSize', SFntSz);
        else
            stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String','2D image', 'BackgroundColor', 'k','ForegroundColor', 'w', 'FontSize', SFntSz);
        end
        
        caxis([Rmin Rmax])
        if sno > 1
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            set(stxthand, 'String', '2D image');
        end

        set(get(gca,'children'),'cdata',squeeze(nii(:,:,S,:)))
        set (gcf, 'ButtonDownFcn', @mouseClick);
        set(get(gca,'Children'),'ButtonDownFcn', @mouseClick);
    end
% -=< Reset view callback function >=-
    function ResetView(object,eventdata)
              
        nii = ImgAx;
        S = S_a;
        sno = sno_a;
%         cla(gcf);
        set(gca,'Position',origPos)
        hdl_im=imshow(squeeze(nii(XImage,YImage,S,MainImage)), [Rmin Rmax]);
        
        if sno > 1
            shand = uicontrol('Visible',toggleSlider,'Style', 'slider','Min',1,'Max',sno,'Value',S,'SliderStep',[1/(sno-1) 10/(sno-1)],'Position', S_Pos,'Callback', {@SliceSlider, nii});
            stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String',sprintf('Slice# %d / %d',S, sno), 'BackgroundColor', 'k','ForegroundColor', 'w', 'FontSize', SFntSz);
        else
            stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String','2D image', 'BackgroundColor', 'k','ForegroundColor', 'w', 'FontSize', SFntSz);
        end
        
        caxis([Rmin Rmax])
        if sno > 1
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            set(stxthand, 'String', '2D image');
        end

        set(get(gca,'children'),'cdata',squeeze(nii(:,:,S,:)))
        set (gcf, 'ButtonDownFcn', @mouseClick);
        set(get(gca,'Children'),'ButtonDownFcn', @mouseClick);
    end

end
