function  ea_imshowpair(Img, options, title, callingfunction, helptext)
% this function is based on IMSHOW3DFULL by Maysam Shahedi and supports
% truecolor images. Windowed view is adapted from MAGNIFY by Rick Hindman.
%
% Todd Herrington, 2016-03-16

if ~exist('callingfunction','var')
   callingfunction='normalization dbs';
end

if ~exist('helptext', 'var')
   helptext = '';
end

switch callingfunction
    case 'coregistration'
        wiresIX=3:5;
    case {'normalization dbs', 'normalization connectome'}
        wiresIX=3;
end

if nargin == 1
    figtit = '';
elseif nargin == 2
    figtit = options.subj.subjId;
elseif nargin >= 3
    figtit = strjoin(cellstr(strvcat({options.subj.subjId, title}))', ', ');
end

isp=figure('color','k','Name',figtit,'NumberTitle','off','MenuBar','none','DockControls','off','ToolBar','none');

% bind drag-drop event
ea_bind_dragndrop(isp, @DropFcn, @DropFcn);

isp.WindowState = 'maximized';

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

if isa(Img,'uint8')
    MaxV = uint8(Inf);
    MinV = uint8(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
elseif isa(Img,'uint16')
    MaxV = uint16(Inf);
    MinV = uint16(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
elseif isa(Img,'uint32')
    MaxV = uint32(Inf);
    MinV = uint32(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
elseif isa(Img,'uint64')
    MaxV = uint64(Inf);
    MinV = uint64(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
elseif isa(Img,'int8')
    MaxV = int8(Inf);
    MinV = int8(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
elseif isa(Img,'int16')
    MaxV = int16(Inf);
    MinV = int16(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
elseif isa(Img,'int32')
    MaxV = int32(Inf);
    MinV = int32(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
elseif isa(Img,'int64')
    MaxV = int64(Inf);
    MinV = int64(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
elseif isa(Img,'logical')
    LevV =0.5;
    Win = 1;
end

ImgO = Img; % ImgO will never be permuted. ImgCr and ImgSg won't be used anymore but generated on the fly via ImgO.

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

BtnSz = 10;

magnFactor = .1;    % standard magnification value that is used then left mouse button is clicked for the first time

[Rmin, Rmax] = WL2R(Win, LevV);

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

showhelptext(callingfunction, helptext);

FigPos = get(gcf,'Position');
BtnStPnt = uint16(FigPos(3)-210)+1;
if BtnStPnt < 360
    BtnStPnt = 360;
end
VAxBtn_Pos = [490 20 15 20];
VSgBtn_Pos = [510 20 15 20];
VCrBtn_Pos = [530 20 15 20];

uicontrol('Style', 'pushbutton','Position', VAxBtn_Pos,'String','A', 'FontSize', BtnSz, 'Callback' , @AxialView);
uicontrol('Style', 'pushbutton','Position', VSgBtn_Pos,'String','S', 'FontSize', BtnSz, 'Callback' , @SagittalView);
uicontrol('Style', 'pushbutton','Position', VCrBtn_Pos,'String','C', 'FontSize', BtnSz, 'Callback' , @CoronalView);

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
        BtnStPnt = uint16(FigPos(3)-210)+1;
        if BtnStPnt < 360
            BtnStPnt = 360;
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
            if ~strcmp(get(f1,'SelectionType'),'alt')
                delete(a2);
            end
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
        end
    end

% -=< Window and level to range conversion >=-
    function [Rmn Rmx] = WL2R(W,L)
        Rmn = L - (W/2);
        Rmx = L + (W/2);
        if (Rmn >= Rmx)
            Rmx = Rmn + 1;
        end
    end

    function KeyPressCallback(object,eventdata)
        H = get(gcf,'UserData');

        if ~isempty(H)
            f1 = H(1); a1 = H(2); a2 = H(3);
            if (strcmp(get(f1,'CurrentCharacter'),'<') || strcmp(get(f1,'CurrentCharacter'),','))
                magnFactor = magnFactor/1.2;
            elseif (strcmp(get(f1,'CurrentCharacter'),'>') || strcmp(get(f1,'CurrentCharacter'),'.'))
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
            end
            set(ImHndl,'cdata',squeeze(Img(XImage,YImage,S,MainImage)));
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
        end
        ButtonMotionCallback(gcf);
    end

    function SwitchModality(numkey,Mod)
        if ~contains(callingfunction,'normalization') % this on
            return
        end
        ea_busyaction('on',gcf,'normcheck');
        PostOpView = 0;
        PostOpLoaded = {''};
        if strcmp(Mod,'alt') % load templates
            SwitchTemplateMod(numkey)
        elseif isempty(Mod)
            anchor = SwitchPatientMod(numkey);
            SwitchTemplateMod(anchor)
        end
        ea_busyaction('off',gcf,'normcheck');
    end

    function SwitchPostop(update)
        if ~contains(callingfunction,'normalization') % this on
            return
        end

        if exist('update','var')
            if ismember(View,PostOpLoaded) % update only, check if Post-OP acquisition is already loaded.
                return
            end
        end

        ea_busyaction('on',gcf,'normcheck');

        switch options.subj.postopModality
            case 'MRI'
                try
                    switch View
                        case 'A'
                            pt=ea_load_nii(options.subj.norm.anat.postop.ax_MRI);
                            PostOpLoaded={'A'};
                        case 'C'
                            pt=ea_load_nii(options.subj.norm.anat.postop.cor_MRI);
                            PostOpLoaded={'C'};
                        case 'S'
                            pt=ea_load_nii(options.subj.norm.anat.postop.sag_MRI);
                            PostOpLoaded={'S'};
                    end
                catch % fallback to transversal acquisition
                    try
                        pt=ea_load_nii(options.subj.norm.anat.postop.ax_MRI);
                        PostOpLoaded={'A','C','S'};
                    catch
                        ea_error('Normalized post-op image seems not available.');
                    end
                end
            case 'CT'
                try
                    pt=ea_load_nii(options.subj.norm.anat.postop.tonemapCT);
                    PostOpLoaded={'A','C','S'};
                catch
                    ea_tonemapct(options, 'norm');
                    try
                        pt=ea_load_nii(options.subj.norm.anat.postop.tonemapCT);
                    catch
                        ea_error('No normalized postoperative acquisition seems available.');
                    end
                end
        end

        if strcmp(options.subj.postopModality, 'CT') % only do windowing for CT
            pt.img=(pt.img-min(pt.img(:)))/(max(pt.img(:)));
            pt.img(pt.img>0.5) = 0.5;
            pt.img=(pt.img-min(pt.img(:)))/(max(pt.img(:)));
        end

        ImgO(:,:,:,1)=pt.img;

        load([ea_space,'wires.mat'], 'wires');
        wires=single(wires);
        wires=wires/255;
        wires=wires.*0.2;
        wires=wires+0.8;

        pt.img=pt.img.*wires; %shows white wires, if commented out, normalizations are shown without the wires, useful if other templates than the MNI are used to normalize images to
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

        switch options.subj.postopModality
            case 'MRI'
                set(gcf, 'Name', regexprep(get(gcf, 'Name'),'(?<=& ).*', 'Postoperative MRI'));
            case 'CT'
                set(gcf, 'Name', regexprep(get(gcf, 'Name'),'(?<=& ).*', 'Postoperative (tonemapped) CT'));
        end

        ea_busyaction('off',gcf,'normcheck');
    end

    function SwitchTemplateMod(numkey)
        if isempty(numkey)
            return;
        else
            numkey = str2double(numkey);
        end

        spacedef = ea_getspacedef;
        if numkey > length(spacedef.templates)
            ea_busyaction('off',gcf,'normcheck');
            return
        elseif (numkey==0 || numkey==length(spacedef.templates)) && spacedef.hasfa
            % When "0" is pressed, choose FA as the backdrop template
            % Also suppose that FA is the last space template
            template = 'fa';
        else
            template = spacedef.templates{numkey};
        end

        tnii = ea_load_nii([ea_space, template, '.nii']);
        tnii.img(:) = zscore(tnii.img(:));
        tnii.img = (tnii.img-min(tnii.img(:)))/(max(tnii.img(:))-min(tnii.img(:)));
        ImgO(:,:,:,2) = tnii.img;

        switch View
            case 'S'
                Img = flip(permute(ImgO, [3 2 1 4]),1);   % Sagittal view image;
            case 'C'
                Img = flip(permute(ImgO, [3 1 2 4]),1);   % Coronal view image;
            case 'A'
                Img = ImgO;
        end

        set(gcf, 'Name', regexprep(get(gcf, 'Name'),'(?<=MNI )\w+', upper(template)));
    end

    function anchor = SwitchPatientMod(numkey)
        anchor = numkey; % just pass on numkey as default.
        numkey = str2double(numkey);

        % load wires first
        load([ea_space,'wires.mat'], 'wires');
        wires = single(wires);
        wires = wires/255;
        wires = wires.*0.2;
        wires = wires+0.8;

        preopCoregImages = struct2cell(options.subj.coreg.anat.preop);

        if numkey>length(preopCoregImages)
            ea_busyaction('off',gcf,'normcheck');
            anchor = [];
            return;
        elseif numkey==0
            if exist([options.root,options.patientname,filesep,'gl',options.prefs.fa2anat],'file')
                pt = ea_load_nii([options.root,options.patientname,filesep,'gl',options.prefs.fa2anat]);
            else
                ea_busyaction('off',gcf,'normcheck');
                return
            end
        elseif numkey==1
            pt = ea_load_nii(options.subj.norm.anat.preop.(options.subj.AnchorModality));
            spacedef = ea_getspacedef;
            anchor = find(ismember(spacedef.templates, options.primarytemplate), 1);
            anchor = num2str(anchor);
        elseif ~isnan(numkey)  % numkey is the index of the preopCoregImages cell
            normAnchor = options.subj.norm.anat.preop.(options.subj.AnchorModality);
            modality = ea_getmodality(preopCoregImages{numkey});
            normImage = strrep(normAnchor, options.subj.AnchorModality, modality);
            if ~isfile(normImage)
                ea_apply_normalization_tofile(options, preopCoregImages{numkey}, normImage, 0);
            end
            pt = ea_load_nii(normImage);
        else % isnan, numkey is actually the full path of the image (used by dragndrop)
            imagePath = anchor; % Take the original numkey (char)
            pt = ea_load_nii(anchor);
            if any(pt.dim~=size(wires))
                msgbox(sprintf('The file you selected seems unnormalized!\nWill try to apply the normalization now.'), 'Warning', 'warn');
                normAnchor = options.subj.norm.anat.preop.(options.subj.AnchorModality);
                modality = ea_getmodality(imagePath);
                normImage = strrep(normAnchor, options.subj.AnchorModality, modality);
                if ~isfile(normImage)
                    ea_apply_normalization_tofile(options, imagePath, normImage, 0);
                end
                pt = ea_load_nii(normImage);
            end
        end

        pt.img = (pt.img-min(pt.img(:)))/(max(pt.img(:)));
        pt.img(pt.img>0.5) = 0.5;
        pt.img = (pt.img-min(pt.img(:)))/(max(pt.img(:)));
        ImgO(:,:,:,1) = pt.img;

        pt.img = pt.img .* wires; % shows white wires, if commented out, normalizations are shown without the wires, useful if other templates than the MNI are used to normalize images to
        ImgO(:,:,:,3) = pt.img;

        switch View
            case 'S'
                Img = flip(permute(ImgO, [3 2 1 4]),1);   % Sagittal view image;
            case 'C'
                Img = flip(permute(ImgO, [3 1 2 4]),1);   % Coronal view image;
            case 'A'
                Img = ImgO;
        end

        set(ImHndl,'cdata',squeeze(Img(XImage,YImage,S,MainImage)))

        modality = ea_getmodality(pt.fname);
        set(gcf, 'Name', regexprep(get(gcf, 'Name'),'(?<=& ).*', ['Preoperative MRI (',modality,')']));
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
        showhelptext(callingfunction, helptext);

        caxis([Rmin Rmax])
        if PostOpView && strcmp(options.subj.postopModality, 'MRI')
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
        showhelptext(callingfunction, helptext);

        caxis([Rmin Rmax])

        if PostOpView && strcmp(options.subj.postopModality, 'MRI')
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
        showhelptext(callingfunction, helptext);

        caxis([Rmin Rmax])

        if PostOpView && strcmp(options.subj.postopModality, 'MRI')
            SwitchPostop('update');
        else
            set(ImHndl,'cdata',squeeze(Img(XImage,YImage,S,MainImage)))
        end
        set(ImHndl,'ButtonDownFcn', @mouseClick);
    end
end

function showhelptext(callingfunction, helptext)
    hold on
    if ~exist('helptext', 'var') || isempty(helptext)
        switch callingfunction
            case 'normalization dbs'
                text(5,5,{'Click to show reference image','Use </> to decrease/increase box size while clicking','Arrow keys / Mouse wheel: Scroll through image','',...
                    '1,2,...: Show preoperative acquisitions [FA=0]','Alt+1,2,...: Switch between available templates [FA=0]','P: Show postoperative acquisitions','',...
                    'Z: Zoom in/out','X: Hybrid view on/off','A: Axial view','C: Coronal view','S: Saggital view',},...
                    'Color','w','HorizontalAlignment','left','VerticalAlignment','top');
            case 'normalization connectome'
                text(5,5,{'Click to show reference image','Use </> to decrease/increase box size while clicking','Arrow keys / Mouse wheel: Scroll through image','',...
                    '1,2,...: Show preoperative acquisitions [FA=0]','Alt+1,2,...: Switch between available templates [FA=0]','',...
                    'Z: Zoom in/out','X: Hybrid view on/off','A: Axial view','C: Coronal view','S: Saggital view',},...
                    'Color','w','HorizontalAlignment','left','VerticalAlignment','top');
            otherwise
                text(5,5,{'Click to show reference image','Use </> to decrease/increase box size while clicking','Arrow keys / Mouse wheel: Scroll through image','',...
                    'Z: Zoom in/out','X: Hybrid view on/off','A: Axial view','C: Coronal view','S: Saggital view'},...
                    'Color','w','HorizontalAlignment','left','VerticalAlignment','top');
        end
    else
        text(5, 5, helptext, 'Color', 'w', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
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
    end
end
% -=< Maysam Shahedi (mshahedi@gmail.com), April 19, 2013>=-
