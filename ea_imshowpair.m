function  ea_imshowpair(Img, options, addstring, callingfunction)
% this function is based on IMSHOW3DFULL by Maysam Shahedi and supports
% truecolor images. Windowed view is adapted from MAGNIFY by Rick Hindman.
%
% Todd Herrington, 2016-03-16

if ~exist('callingfunction','var')
   callingfunction='normalization dbs';
end

switch callingfunction
    case 'coregistration'
        wiresIX=3:5;
        gridIX=nan;
    case {'normalization dbs', 'normalization connectome'}
        wiresIX=3;
        gridIX=4;
end

if nargin == 1
    figtit = '';
elseif nargin == 2
    figtit = options.patientname;
elseif nargin >= 3
    figtit = [options.patientname,', ',addstring];
end

isp=figure('color','k','Name',figtit,'NumberTitle','off','MenuBar','none','DockControls','off','ToolBar','none');

% bind drag-drop event
dndobj = ea_bind_dragndrop(isp, @DropFcn, @DropFcn);

ea_maximize(isp);
Img=single(Img);
sno = size(Img);  % image size
sno_a = sno(3);  % number of axial slices
S_a = round(sno_a/2);
sno_s = sno(2);  % number of sagittal slices
S_s = round(sno_s/2);
sno_c = sno(1);  % number of coronal slices
S_c = round(sno_c/2);
S = S_a;
sno = sno_a;
PostOpView=0;
PostOpLoaded={''};

global InitialCoord;

MinV = 0;
MaxV = ea_nanmax(Img(:));
LevV = (double( MaxV) + double(MinV)) / 2;
Win = double(MaxV) - double(MinV);
WLAdjCoe = (Win + 1)/1024;
FineTuneC = [1 1/16];    % Regular/Fine-tune mode coefficients

if isa(Img,'uint8')
    MaxV = uint8(Inf);
    MinV = uint8(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'uint16')
    MaxV = uint16(Inf);
    MinV = uint16(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'uint32')
    MaxV = uint32(Inf);
    MinV = uint32(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'uint64')
    MaxV = uint64(Inf);
    MinV = uint64(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int8')
    MaxV = int8(Inf);
    MinV = int8(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int16')
    MaxV = int16(Inf);
    MinV = int16(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int32')
    MaxV = int32(Inf);
    MinV = int32(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int64')
    MaxV = int64(Inf);
    MinV = int64(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'logical')
    MaxV = 0;
    MinV = 1;
    LevV =0.5;
    Win = 1;
    WLAdjCoe = 0.1;
end

ImgO = Img; % ImgO will never be permuted. ImgCr and ImgSg won't be used anymore but generated on the fly via ImgO.
%ImgCr = flip(permute(Img, [3 1 2 4]),1);   % Coronal view image
%ImgSg = flip(permute(Img, [3 2 1 4]),1);   % Sagittal view image

ImgZ=0; % zoomed or unzoomed state
ImgZax{1}=1:size(ImgO,1); % axial zoomed out boundingboxes
ImgZax{2}=1:size(ImgO,2); % axial zoomed out boundingboxes
ImgZax{3}=round(size(ImgO,1)/4):round((size(ImgO,1)/4)*3); % axial zoomed boundingboxes
ImgZax{4}=round(size(ImgO,2)/4):round((size(ImgO,2)/4)*3); % axial zoomed boundingboxes
ImgZcr{1}=1:size(ImgO,3); % coronal zoomed out boundingboxes
ImgZcr{2}=1:size(ImgO,1); % coronal zoomed out boundingboxes
ImgZcr{3}=round(size(ImgO,3)/4):round((size(ImgO,3)/4)*3); % coronal zoomed boundingboxes
ImgZcr{4}=round(size(ImgO,1)/4):round((size(ImgO,1)/4)*3); % coronal zoomed boundingboxes
ImgZsg{1}=1:size(ImgO,3); % saggital zoomed out boundingboxes
ImgZsg{2}=1:size(ImgO,2); % saggital zoomed out boundingboxes
ImgZsg{3}=round(size(ImgO,3)/4):round((size(ImgO,3)/4)*3); % saggital zoomed boundingboxes
ImgZsg{4}=round(size(ImgO,2)/4):round((size(ImgO,2)/4)*3); % saggital zoomed boundingboxes

View = 'A';

SFntSz = 9;
LFntSz = 10;
WFntSz = 10;
VwFntSz = 10;
LVFntSz = 9;
WVFntSz = 9;
BtnSz = 10;
ChBxSz = 10;

magnFactor = .1;    % standard magnification value that is used then left mouse button is clicked for the first time

[Rmin Rmax] = WL2R(Win, LevV);

hdl_im = axes('position',[0,0,1,1]);
set(0,'CurrentFigure',isp);
MainImage = 1;
XImage=1:size(Img,1);
YImage=1:size(Img,2);

try % image toolbox
    ImHndl=imshow(squeeze(Img(XImage,YImage,S,MainImage)), [Rmin Rmax]);
catch
    ImHndl=imagesc(squeeze(Img(XImage,YImage,S,MainImage)), [Rmin Rmax]);
end

showhelptext(callingfunction);

FigPos = get(gcf,'Position');
S_Pos = [50 20 uint16(FigPos(3)-150)+1 20];
Stxt_Pos = [50 90 uint16(FigPos(3)-100)+1 15];
Wtxt_Pos = [20 20 60 20];
Wval_Pos = [75 20 60 20];
Ltxt_Pos = [140 20 45 20];
Lval_Pos = [180 20 60 20];
BtnStPnt = uint16(FigPos(3)-210)+1;
if BtnStPnt < 360
    BtnStPnt = 360;
end
Btn_Pos = [BtnStPnt 20 80 20];
ChBx_Pos = [BtnStPnt+90 20 100 20];
Vwtxt_Pos = [255 20 35 20];
VAxBtn_Pos = [490 20 15 20];
VSgBtn_Pos = [510 20 15 20];
VCrBtn_Pos = [530 20 15 20];

if sno > 1
%    shand = uicontrol('Style', 'slider','Min',1,'Max',sno,'Value',S,'SliderStep',[1/(sno-1) 10/(sno-1)],'Position', S_Pos,'Callback', {@SliceSlider, Img});
%    stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String',sprintf('Slice# %d / %d',S, sno), 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
else
%    stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String','2D image', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
end
%ltxthand = uicontrol('Style', 'text','Position', Ltxt_Pos,'String','Level: ', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', LFntSz);
%wtxthand = uicontrol('Style', 'text','Position', Wtxt_Pos,'String','Window: ', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', WFntSz);
%lvalhand = uicontrol('Style', 'edit','Position', Lval_Pos,'String',sprintf('%6.0f',LevV), 'BackgroundColor', [1 1 1], 'FontSize', LVFntSz,'Callback', @WinLevChanged);
%wvalhand = uicontrol('Style', 'edit','Position', Wval_Pos,'String',sprintf('%6.0f',Win), 'BackgroundColor', [1 1 1], 'FontSize', WVFntSz,'Callback', @WinLevChanged);
%Btnhand = uicontrol('Style', 'pushbutton','Position', Btn_Pos,'String','Auto W/L', 'FontSize', BtnSz, 'Callback' , @AutoAdjust);
%ChBxhand = uicontrol('Style', 'checkbox','Position', ChBx_Pos,'String','Fine Tune', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', ChBxSz);
%Vwtxthand = uicontrol('Style', 'text','Position', Vwtxt_Pos,'String','View: ', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', LFntSz);
VAxBtnhand = uicontrol('Style', 'pushbutton','Position', VAxBtn_Pos,'String','A', 'FontSize', BtnSz, 'Callback' , @AxialView);
VSgBtnhand = uicontrol('Style', 'pushbutton','Position', VSgBtn_Pos,'String','S', 'FontSize', BtnSz, 'Callback' , @SagittalView);
VCrBtnhand = uicontrol('Style', 'pushbutton','Position', VCrBtn_Pos,'String','C', 'FontSize', BtnSz, 'Callback' , @CoronalView);

set (gcf, 'WindowScrollWheelFcn', @mouseScroll);
set (gcf, 'ButtonDownFcn', @mouseClick);
set(ImHndl,'ButtonDownFcn', @mouseClick);
set(gcf,'WindowButtonUpFcn', @mouseRelease)
set(gcf,'ResizeFcn', @figureResized)
set(gcf,'WindowButtonMotionFcn', @ButtonMotionCallback)
set(gcf,'KeyPressFcn', @KeyPressCallback);


% -=< Figure resize callback function >=-
    function figureResized(object, eventdata)
        FigPos = get(gcf,'Position');
        S_Pos = [50 45 uint16(FigPos(3)-100)+1 20];
        Stxt_Pos = [50 65 uint16(FigPos(3)-100)+1 15];
        BtnStPnt = uint16(FigPos(3)-210)+1;
        if BtnStPnt < 360
            BtnStPnt = 360;
        end
        Btn_Pos = [BtnStPnt 20 80 20];
        ChBx_Pos = [BtnStPnt+90 20 100 20];
        if sno > 1
           % set(shand,'Position', S_Pos);
        end
        %set(stxthand,'Position', Stxt_Pos);
        %set(ltxthand,'Position', Ltxt_Pos);
        %set(wtxthand,'Position', Wtxt_Pos);
        %set(lvalhand,'Position', Lval_Pos);
        %set(wvalhand,'Position', Wval_Pos);
        %set(Btnhand,'Position', Btn_Pos);
        %set(ChBxhand,'Position', ChBx_Pos);
        %set(Vwtxthand,'Position', Vwtxt_Pos);
        %set(VAxBtnhand,'Position', VAxBtn_Pos);
        %set(VSgBtnhand,'Position', VSgBtn_Pos);
        %set(VCrBtnhand,'Position', VCrBtn_Pos);
    end

% -=< Slice slider callback function >=-
    function SliceSlider (hObj,event, Img)
        S = round(get(hObj,'Value'));
        set(ImHndl,'cdata',squeeze(Img(:,:,S,:)))
        caxis([Rmin Rmax])
        if sno > 1
            %set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            %set(stxthand, 'String', '2D image');
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
         %   set(shand,'Value',S);
%            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
 %           set(stxthand, 'String', '2D image');
        end
        set(ImHndl,'cdata',squeeze(Img(XImage,YImage,S,MainImage)))
    end

% -=< Mouse button released callback function >=-
    function mouseRelease (object,eventdata)
        set(gcbf, 'WindowButtonMotionFcn', '')
        H = get(object,'UserData');
        if (isempty(H)) % presumed RIGHT CLICK
            % do nothing
        else % presumed LEFT CLICK
            f1 = H(1); a1 = H(2); a2 = H(3);
            set(a1, 'Color',get(a2,'Color'));
            set(f1, ...
                'UserData',[], ...
                'Pointer','arrow', ...
                'CurrentAxes',a1);
            if ~strcmp(get(f1,'SelectionType'),'alt'),
                delete(a2);
            end;
        end
    end

% -=< Mouse click callback function >=-
    function mouseClick (object, eventdata)
        MouseStat = get(gcbf, 'SelectionType');
        if (MouseStat(1) == 'a')        %   RIGHT CLICK
            InitialCoord = get(0,'PointerLocation');
            %             set(gcf, 'WindowButtonMotionFcn', @WinLevAdj);
        else    % assumed LEFT CLICK
            i1 = object;
            a1 = get(i1,'Parent');
            f1 = get(a1,'Parent');

            a2 = copyobj(a1,f1, 'legacy');
            % helptext is also copied, we delete it here so it is not
            % duplicated when when magnifying
            delete(findall(a2,'Type','text')); 
            
            i2 = get(a2,'Children');
            try
            i2=i2(2);
            end

            set(f1, ...
                'UserData',[f1,a1,a2], ...
                'Pointer','crosshair', ...
                'CurrentAxes',a2);
            set(a2, ...
                'UserData',[1,magnFactor], ...  %magnification, frame size
                'Color',get(a1,'Color'), ...
                'Box','on');
            set(i2,'CData',Img(XImage,YImage,S,2));
            set(f1,'WindowButtonMotionFcn', @ButtonMotionCallback)
            ButtonMotionCallback(f1);

        end
    end

    function ButtonMotionCallback(object,eventdata)
        H = get(object,'UserData');
        if ~isempty(H)
            f1 = H(1); a1 = H(2); a2 = H(3);
            a2_param = get(a2,'UserData');
            f_pos = get(f1,'Position');
            a1_pos = get(a1,'Position');

            [f_cp, a1_cp] = pointer2d(f1,a1);

            set(a2,'Position',[(f_cp./f_pos(3:4)) 0 0]+a2_param(2)*a1_pos(3)*[-1 -1 2 2]);
            a2_pos = get(a2,'Position');

            set(a2,'XLim',a1_cp(1)+(1/a2_param(1))*(a2_pos(3)/a1_pos(3))*diff(get(a1,'XLim'))*[-0.5 0.5]);
            set(a2,'YLim',a1_cp(2)+(1/a2_param(1))*(a2_pos(4)/a1_pos(4))*diff(get(a1,'YLim'))*[-0.5 0.5]);
        end;
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
        Win = double(max(Img(:))-min(Img(:)));
        Win (Win < 1) = 1;
        LevV = double(min(Img(:)) + (Win/2));
        [Rmin, Rmax] = WL2R(Win,LevV);
        caxis([Rmin, Rmax])
        set(lvalhand, 'String', sprintf('%6.0f',LevV));
        set(wvalhand, 'String', sprintf('%6.0f',Win));
    end

    function KeyPressCallback(object,eventdata)
        H = get(gcf,'UserData');

        if ~isempty(H)
            f1 = H(1); a1 = H(2); a2 = H(3);
            if (strcmp(get(f1,'CurrentCharacter'),'<') | strcmp(get(f1,'CurrentCharacter'),','))
                magnFactor = magnFactor/1.2;
            elseif (strcmp(get(f1,'CurrentCharacter'),'>') | strcmp(get(f1,'CurrentCharacter'),'.'))
                magnFactor = magnFactor*1.2;
            end
            set(a2,'UserData',[1, magnFactor]);
        end

        if (strcmp(eventdata.Key,'leftarrow') || strcmp(eventdata.Key,'downarrow'))
            ev = []; ev.VerticalScrollCount = -1;
            mouseScroll (gcf, ev);
        elseif (strcmp(eventdata.Key,'rightarrow') || strcmp(eventdata.Key,'uparrow'))
            ev = []; ev.VerticalScrollCount = 1;
            mouseScroll (gcf, ev);
        elseif (strcmpi(eventdata.Key,'c'))
            CoronalView([]);
        elseif (strcmpi(eventdata.Key,'a'))
            AxialView([]);
        elseif (strcmpi(eventdata.Key,'s'))
            SagittalView([]);
        elseif (strcmpi(eventdata.Key,'P') && strcmp(callingfunction,'normalization dbs'))
            SwitchPostop;
        elseif ~isnan(str2double(eventdata.Character))
            SwitchModality(eventdata.Key,eventdata.Modifier)
        elseif (strcmpi(eventdata.Key,'x'))
            if MainImage(1)==1
                MainImage=wiresIX;
            elseif MainImage(1)==wiresIX(1)
                MainImage=1;
            elseif MainImage(1)==gridIX
                MainImage=wiresIX;
            end
            set(ImHndl,'cdata',squeeze(Img(XImage,YImage,S,MainImage)));
                        
        elseif (strcmpi(eventdata.Key,'g'))
            if size(Img,4)==4 && strfind(callingfunction,'normalization') % only do if grid is available.
                if MainImage(1)==1
                    MainImage=gridIX;
                elseif MainImage(1)==gridIX
                    MainImage=1;
                elseif MainImage(1)==wiresIX
                    MainImage=gridIX;
                end
                set(ImHndl,'cdata',squeeze(Img(XImage,YImage,S,MainImage)));
            end
        elseif (strcmpi(eventdata.Key,'z')) % toggles zoom in/out
            ImgZ=~ImgZ;
            switch View
                case 'A'
                    AxialView([]);
                case 'C'
                    CoronalView([]);
                case 'S'
                    SagittalView([]);
            end
        end;
        ButtonMotionCallback(gcf);
    end

    function SwitchModality(numkey,Mod)
        if ~strfind(callingfunction,'normalization') % this on
            return
        end
        %[options] = ea_assignpretra(options);
        ea_busyaction('on',gcf,'normcheck');
        PostOpView=0;
        PostOpLoaded={''};
        if strcmp(Mod,'alt') % load templates
            SwitchTemplateMod(numkey)
        elseif isempty(Mod)
            anchor=SwitchPatientMod(numkey);
            SwitchTemplateMod(anchor)
        end
        ea_busyaction('off',gcf,'normcheck');
    end

    function SwitchPostop(update)
        if ~strfind(callingfunction,'normalization') % this on
            return
        end

        if exist('update','var')
            if ismember(View,PostOpLoaded) % update only, check if Post-OP acquisition is already loaded.
                return
            end
        end

        ea_busyaction('on',gcf,'normcheck');

        switch options.modality
            case 1 % MR
                try
                    switch View
                        case 'A'
                            pt=ea_load_nii([options.root,options.patientname,filesep,options.prefs.gtranii]);
                            PostOpLoaded={'A'};
                        case 'C'
                            pt=ea_load_nii([options.root,options.patientname,filesep,options.prefs.gcornii]);
                            PostOpLoaded={'C'};
                        case 'S'
                            pt=ea_load_nii([options.root,options.patientname,filesep,options.prefs.gsagnii]);
                            PostOpLoaded={'S'};
                    end
                catch % fallback to transversal acquisition
                    try
                        pt=ea_load_nii([options.root,options.patientname,filesep,options.prefs.gtranii]);
                        PostOpLoaded={'A','C','S'};
                    catch
                        ea_error('No normalized postoperative acquisition seems available.');
                    end
                end
            case 2 % CT
                try
                    pt=ea_load_nii([options.root,options.patientname,filesep,'tp_',options.prefs.gctnii]);
                    PostOpLoaded={'A','C','S'};
                catch
                    ea_tonemapct(options, 'norm');
                    try
                        pt=ea_load_nii([options.root,options.patientname,filesep,'tp_',options.prefs.gctnii]);
                    catch
                        ea_error('No normalized postoperative acquisition seems available.');
                    end
                end
        end
        if options.modality==2 % only do windowing for MR
            pt.img=(pt.img-min(pt.img(:)))/(max(pt.img(:)));
            pt.img(pt.img>0.5) = 0.5;
            pt.img=(pt.img-min(pt.img(:)))/(max(pt.img(:)));
        end
        ImgO(:,:,:,1)=pt.img;

        w=load([ea_space(options),'wires.mat']);
        w.wires=single(w.wires);
        w.wires=w.wires/255;
        w.wires=w.wires.*0.2;
        w.wires=w.wires+0.8;

        pt.img=pt.img.*w.wires; %shows white wires, if commented out, normalizations are shown without the wires, useful if other templates than the MNI are used to normalize images to
        ImgO(:,:,:,3)=pt.img;

        switch View
            case 'S'
                Img = flip(permute(ImgO, [3 2 1 4]),1);   % Sagittal view image;
            case 'C'
                Img = flip(permute(ImgO, [3 1 2 4]),1);   % Coronal view image;
            case 'A'
                Img = ImgO;
        end
        PostOpView=1;
        if ~MainImage==wiresIX
            set(1,'cdata',squeeze(Img(XImage,YImage,S,MainImage)))
        else
            try
                set(ImHndl,'cdata',squeeze(Img(XImage,YImage,S,MainImage)))
            catch
                keyboard
            end
        end
        SwitchTemplateMod('2') % set to T2 if available.

        [~, fname] = ea_niifileparts(pt.fname);
        switch options.modality
            case 1
                set(gcf, 'Name', regexprep(get(gcf, 'Name'),'(?<=& ).*', ['Postoperative MRI (',fname,')']));
            case 2
                set(gcf, 'Name', regexprep(get(gcf, 'Name'),'(?<=& ).*', ['Postoperative (tonemapped) CT (',fname,')']));
        end

        ea_busyaction('off',gcf,'normcheck');
    end

    function SwitchTemplateMod(numkey)
        if isempty(numkey)
            return
        end

        options.sd=load([ea_space,'ea_space_def.mat'],'spacedef');
        if str2double(numkey)>length(options.sd.spacedef.templates)
            ea_busyaction('off',gcf,'normcheck');
            return
        end
        if str2double(numkey)==0
            if options.sd.spacedef.hasfa
                whichtemplate='fa.nii';
            else
                ea_busyaction('off',gcf,'normcheck');
                return
            end
        else
            whichtemplate=options.sd.spacedef.templates{str2double(numkey)};
        end
        tnii=ea_load_nii([ea_space,whichtemplate]);
        tnii.img(:)=zscore(tnii.img(:));
        tnii.img=(tnii.img-min(tnii.img(:)))/(max(tnii.img(:))-min(tnii.img(:)));
        ImgO(:,:,:,2)=tnii.img;
        switch View
            case 'S'
                Img = flip(permute(ImgO, [3 2 1 4]),1);   % Sagittal view image;
            case 'C'
                Img = flip(permute(ImgO, [3 1 2 4]),1);   % Coronal view image;
            case 'A'
                Img = ImgO;
        end
        [~, fname] = fileparts(whichtemplate);
        set(gcf, 'Name', regexprep(get(gcf, 'Name'),'(?<=MNI )\w+', upper(fname)));
    end

    function anchor=SwitchPatientMod(numkey)
        anchor=numkey; % just pass on numkey as default.

        % load wires first
        w=load([ea_space(options),'wires.mat']);
        w.wires=single(w.wires);
        w.wires=w.wires/255;
        w.wires=w.wires.*0.2;
        w.wires=w.wires+0.8;

        [options,presentfiles] = ea_assignpretra(options);
        if str2double(numkey)==1
            pt=ea_load_nii([options.root,options.patientname,filesep,options.prefs.gprenii]);
            sd=load([ea_space,'ea_space_def.mat']);
            [~,anchor]=ismember(options.primarytemplate,sd.spacedef.templates); % in this case, "1" could e.g. be the T2 as well, so pass on the T2 in this case to also load the correct template.
            anchor=num2str(anchor);
        elseif str2double(numkey)>length(presentfiles)
            ea_busyaction('off',gcf,'normcheck');
            anchor = [];
            return
        elseif str2double(numkey)==0
            if exist([options.root,options.patientname,filesep,'gl',options.prefs.fa2anat],'file')
                pt=ea_load_nii([options.root,options.patientname,filesep,'gl',options.prefs.fa2anat]);
            else
                ea_busyaction('off',gcf,'normcheck');
                return
            end
        elseif ~isnan(str2double(numkey))  % numkey is the index of the image in presentfiles cell
            directory=[options.root,options.patientname,filesep];
            toload=[directory,'gl',presentfiles{str2double(numkey)}];
            if ~exist(toload,'file')
                from{1}=[directory,presentfiles{str2double(numkey)}];
                to{1}=toload;
                ea_apply_normalization_tofile(options,from,to,directory,0);
            end
            pt=ea_load_nii(toload);
        else    % numkey is actually the full path of the image (used by dragndrop)
            pt=ea_load_nii(numkey);
            if any(pt.dim~=size(w.wires))
                msgbox(sprintf('The file you selected seems unnormalized!\nWill try to apply the normalization now.'),'Warning','warn');
                [fpath, fname, fext] = ea_niifileparts(numkey);
                from{1}=numkey;
                to{1}=[fileparts(fpath), filesep, 'gl', fname, fext];
                ea_apply_normalization_tofile(options,from,to,fileparts(fpath),0);
                pt=ea_load_nii(to{1});
            end
        end

        pt.img=(pt.img-min(pt.img(:)))/(max(pt.img(:)));
        pt.img(pt.img>0.5) = 0.5;
        pt.img=(pt.img-min(pt.img(:)))/(max(pt.img(:)));
        ImgO(:,:,:,1)=pt.img;

        pt.img=pt.img.*w.wires; % shows white wires, if commented out, normalizations are shown without the wires, useful if other templates than the MNI are used to normalize images to
        ImgO(:,:,:,3)=pt.img;

        switch View
            case 'S'
                Img = flip(permute(ImgO, [3 2 1 4]),1);   % Sagittal view image;
            case 'C'
                Img = flip(permute(ImgO, [3 1 2 4]),1);   % Coronal view image;
            case 'A'
                Img = ImgO;
        end

        if ~MainImage==wiresIX
            set(ImHndl,'cdata',squeeze(Img(XImage,YImage,S,MainImage)))
        else
            set(ImHndl,'cdata',squeeze(Img(XImage,YImage,S,MainImage)))
        end

        [~, fname] = ea_niifileparts(pt.fname);
        set(gcf, 'Name', regexprep(get(gcf, 'Name'),'(?<=& ).*', ['Preoperative MRI (',fname,')']));
    end

    % DragnDrop callback function
    function DropFcn(~, eventdata)
        switch eventdata.DropType
            case 'file' % drag file into imshowpair window
                SwitchPatientMod(eventdata.Data{1});
            case 'string'   % % drag string (full path of the file) into imshowpair window
                [fpath, fname] = ea_niifileparts(eventdata.Data);
                fullpath = ea_regexpdir(fileparts(fpath), ['^',fname,'(\.nii|\.nii\.gz)?$'],0);
                if ~isempty(fullpath)
                    SwitchPatientMod(eventdata.Data);
                end
        end
    end

% -=< Axial view callback function >=-
    function AxialView(object,eventdata)
        if View == 'S'
            S_s = S;
        elseif View == 'C'
            S_c = S;
        end
        View = 'A';
        Img = ImgO;
        S = S_a;
        sno = sno_a;
        cla(hdl_im);
        hdl_im = axes('position',[0,0,1,1]);

        if ~ImgZ
            XImage=ImgZax{1};
            YImage=ImgZax{2};
        else
            XImage=ImgZax{3};
            YImage=ImgZax{4};
        end

        try % image toolbox
            ImHndl=imshow(squeeze(Img(XImage,YImage,S,1)), [Rmin Rmax]);
        catch
            ImHndl=imagesc(squeeze(Img(XImage,YImage,S,1)), [Rmin Rmax]);
        end
        showhelptext(callingfunction);

        if sno > 1
        %    shand = uicontrol('Style', 'slider','Min',1,'Max',sno,'Value',S,'SliderStep',[1/(sno-1) 10/(sno-1)],'Position', S_Pos,'Callback', {@SliceSlider, Img});
            %stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String',sprintf('Slice# %d / %d',S, sno), 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
        else
            %stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String','2D image', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
        end

        caxis([Rmin Rmax])
        if sno > 1
            %set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            %set(stxthand, 'String', '2D image');
        end
        if PostOpView && options.modality==1
            SwitchPostop('update');
        else
            set(ImHndl,'cdata',squeeze(Img(XImage,YImage,S,MainImage)))
        end
        set(ImHndl,'ButtonDownFcn', @mouseClick);
    end

% -=< Sagittal view callback function >=-
    function SagittalView(object,eventdata)
        if View == 'A'
            S_a = S;
        elseif View == 'C'
            S_c = S;
        end
        View = 'S';

        if ~ImgZ
            XImage=ImgZsg{1};
            YImage=ImgZsg{2};
        else
            XImage=ImgZsg{3};
            YImage=ImgZsg{4};
        end

        Img = flip(permute(ImgO, [3 2 1 4]),1);   % Sagittal view image;
        S = S_s;
        sno = sno_s;
        cla(hdl_im);
        hdl_im = axes('position',[0,0,1,1]);

        try % image toolbox
            ImHndl=imshow(squeeze(Img(XImage,YImage,S,1)), [Rmin Rmax]);
        catch
            ImHndl=imagesc(squeeze(Img(XImage,YImage,S,1)), [Rmin Rmax]);
        end
        showhelptext(callingfunction);

        if sno > 1
          %  shand = uicontrol('Style', 'slider','Min',1,'Max',sno,'Value',S,'SliderStep',[1/(sno-1) 10/(sno-1)],'Position', S_Pos,'Callback', {@SliceSlider, Img});
            %stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String',sprintf('Slice# %d / %d',S, sno), 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
        else
            %stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String','2D image', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
        end

        caxis([Rmin Rmax])
        if sno > 1
%            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
%            set(stxthand, 'String', '2D image');
        end
        if PostOpView && options.modality==1
            SwitchPostop('update');
        else
            set(ImHndl,'cdata',squeeze(Img(XImage,YImage,S,MainImage)));
        end
        set(ImHndl,'ButtonDownFcn', @mouseClick);
    end

% -=< Coronal view callback function >=-
    function CoronalView(object,eventdata)
        if View == 'A'
            S_a = S;
        elseif View == 'S'
            S_s = S;
        end
        View = 'C';


        if ~ImgZ
            XImage=ImgZcr{1};
            YImage=ImgZcr{2};
        else
            XImage=ImgZcr{3};
            YImage=ImgZcr{4};
        end

        Img = flip(permute(ImgO, [3 1 2 4]),1);   % Coronal view image;
        S = S_c;
        sno = sno_c;
        cla(hdl_im);
        hdl_im = axes('position',[0,0,1,1]);

        try % image toolbox
            ImHndl=imshow(squeeze(Img(XImage,YImage,S,1)), [Rmin Rmax]);
        catch
            ImHndl=imagesc(squeeze(Img(XImage,YImage,S,1)), [Rmin Rmax]);
        end
        showhelptext(callingfunction);

        if sno > 1
          %  shand = uicontrol('Style', 'slider','Min',1,'Max',sno,'Value',S,'SliderStep',[1/(sno-1) 10/(sno-1)],'Position', S_Pos,'Callback', {@SliceSlider, Img});
            %stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String',sprintf('Slice# %d / %d',S, sno), 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
        else
            %stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String','2D image', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
        end

        caxis([Rmin Rmax])

        if PostOpView && options.modality==1
            SwitchPostop('update');
        else
            set(ImHndl,'cdata',squeeze(Img(XImage,YImage,S,MainImage)))
        end
        set(ImHndl,'ButtonDownFcn', @mouseClick);
    end
end

function showhelptext(callingfunction)
    hold on
    if strcmp(callingfunction,'normalization dbs')
        helptext=text(5,5,{'Click to show reference image','Use </> to decrease/increase box size while clicking','Arrow keys / Mouse wheel: Scroll through image','',...
            '1,2,...: Show preoperative acquisitions [FA=0]','Alt+1,2,...: Switch between available templates [FA=0]','P: Show postoperative acquisitions','',...
            'Z: Zoom in/out','X: Hybrid view on/off','A: Axial view','C: Coronal view','S: Saggital view',},...
            'Color','w','HorizontalAlignment','left','VerticalAlignment','top');
    elseif strcmp(callingfunction,'normalization connectome')
        helptext=text(5,5,{'Click to show reference image','Use </> to decrease/increase box size while clicking','Arrow keys / Mouse wheel: Scroll through image','',...
            '1,2,...: Show preoperative acquisitions [FA=0]','Alt+1,2,...: Switch between available templates [FA=0]','',...
            'Z: Zoom in/out','X: Hybrid view on/off','A: Axial view','C: Coronal view','S: Saggital view',},...
            'Color','w','HorizontalAlignment','left','VerticalAlignment','top');
    else
        helptext=text(5,5,{'Click to show reference image','Use </> to decrease/increase box size while clicking','Arrow keys / Mouse wheel: Scroll through image','',...
            'Z: Zoom in/out','X: Hybrid view on/off','A: Axial view','C: Coronal view','S: Saggital view'},...
            'Color','w','HorizontalAlignment','left','VerticalAlignment','top');
    end
end

% Included for completeness (usually in own file)
function [fig_pointer_pos, axes_pointer_val] = pointer2d(fig_hndl,axes_hndl)
    %
    %pointer2d(fig_hndl,axes_hndl)
    %
    %	Returns the coordinates of the pointer (in pixels)
    %	in the desired figure (fig_hndl) and the coordinates
    %       in the desired axis (axes coordinates)
    %
    % Example:
    %  figure(1),
    %  hold on,
    %  for i = 1:1000,
    %     [figp,axp]=pointer2d;
    %     plot(axp(1),axp(2),'.','EraseMode','none');
    %     drawnow;
    %  end;
    %  hold off

    % Rick Hindman - 4/18/01

    if (nargin == 0)
        fig_hndl = gcf;
        axes_hndl = gca;
    elseif (nargin == 1)
        axes_hndl = get(fig_hndl,'CurrentAxes');
    end

    set(fig_hndl,'Units','pixels');

    pointer_pos = get(0,'PointerLocation');	%pixels {0,0} lower left
    fig_pos = get(fig_hndl,'Position');	%pixels {l,b,w,h}

    fig_pointer_pos = pointer_pos - fig_pos([1,2]);
    set(fig_hndl,'CurrentPoint',fig_pointer_pos);

    if (isempty(axes_hndl))
        axes_pointer_val = [];
    elseif (nargout == 2)
        axes_pointer_line = get(axes_hndl,'CurrentPoint');
        axes_pointer_val = sum(axes_pointer_line)/2;
    end;
end
% -=< Maysam Shahedi (mshahedi@gmail.com), April 19, 2013>=-
