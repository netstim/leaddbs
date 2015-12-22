function mesoGT_tool(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  mesoGT_tool(commandstr1,varargin1,commandstr2,varargin2,.....)
% %
% %
% %  Marco Reisert, Medical Physics, University Hospital Freiburg
% %  03/10
% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  All functionalities of the tool can be accessed via the GUI or the
% %  commands given below.
% %
% %  Note: This program is memory intensive. We recommend to have at
% %        least 8 Gigabyte of free memory.
% %
% % -----------------------------------------------------------------
% % commandstr == 'loadData'
% %     mesoGT_tool('loadData',type,dwifile,gradfile,maskfile,threshold)
% %     loads dwi-data to the internal datastructure of the tracker.
% %      - type == 'mat'
% %           * dwifile , a mrstruct file containing q-space data
% %           * gradfile , [], gradinfo is in dwifile.user.bTensor
% %           * maskfile , a mrstruct or nifti containing a WM-mask
% %           * threshold , if maskfile is a prbablitic segmentation, the threshold for creating the mask
% %      - type == 'nii'
% %           * dwifile, a 4D-nifti, or a cellarray of 3D-niftis containing q-space data
% %           * gradfile, a cellarray of two files containing gradient directions and b.values (FSL style)
% %           * maskfile, a nifti containing a WM-mask
% %           * threshold, dito
% % commandstr == 'setData'
% %     mesoGT_tool('setData',data,tensor,mask,vox,edges,name)
% %     - data , a [w h d N] array containg the dwi data
% %     - tensor , a [3 3 N] array containing q-space information such that trace(tensor(:,:,k)) = bval(k)
% %     - mask , a [w h d] array containing the WM-mask
% %     - vox , a [3 1] with voxelsizes in mm
% %     - edges , a [4 4] matrix with voxel->worldcoord info (homogenous coords)
% %     - name , a name
% % -----------------------------------------------------------------
% % commandstr == 'loadHCP'
% %     mesoGT_tool('loadHCP',subjfolder)
% %     loads a hcp subject. Assumes <subjfolder>/T1w/wmparc.nii.gz
% %     and <subjfolder>/T1w/Diffusion/data.nii.gz(bvals,bvecs) is present
% %     You need at least 16 GB of free RAM
% % ----------------------------------------------------------------
% % commandstr == 'loadFTR'
% %     mesoGT_tool('setFTR',filename)
% %      - filename  , filename of ftr-file (the fibertracks) formerly created with mesoGT_tool
% %
% %
% % commandstr == 'setFTRname'
% %     mesoGT_tool('setFTR',filename)
% %      - filename  , name of file the ftr is going to be saved
% % ----------------------------------------------------------------
% % commandstr == 'start'
% %     start tracking
% %
% % commandstr == 'reset'
% %     reset tracking state
% %
% % commandstr == 'saveFD'
% %     save fiber densities
% %
% % ----------------------------------------------------------------
% %
% % commandstr == 'setparam'
% %     mesoGT_tool('setparam',parametername,value);
% %       possible parameternames are shown with 'showParameters'
% %
% % commandstr == 'showParameters'
% %     show current parameters
% %
% %
% % commandstr == 'reactivate'
% %    reactive uicontrols after crash
% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



hmain = findobj('Tag','fiberGT_main');
if isempty(hmain),
    hmain = figure('Name','MesoTracker','Tag','fiberGT_main','MenuBar','none','NumberTitle','off', ...
                   'resize','off','CloseRequestFcn',@my_closereq);
    poshm = get(hmain,'Position');
    sz = [800 550];
    set(hmain,'Position',[(poshm(1:2)-(sz-poshm(3:4))) sz]);


    [params firstParams] = GTdefaults;


     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     pos = [5 300 130 20];
     twid = [140 0 -50 0];
     thei = -[0 25 0 0];
     uicontrol('Style','text','String','Parameters','Position',pos + [0 -235 130 260],'HorizontalAlignment','left','FontWeight','bold');
     for k = 1:length(firstParams),
         uicontrol('Style','text','String',firstParams(k).name,'Position',pos+[5 0 0 0],'Tag','fiberGT_text','HorizontalAlignment','left');
         uicontrol('Style','edit','String',num2str(getfield(params,firstParams(k).tag)),'Position',pos + twid,'Tag',['fiberGT_edit_' firstParams(k).tag] ,'BackGroundColor',[1 1 1]);
         pos = pos + thei;
     end;

    pos = [10 35 150 25];
    uicontrol('Style','pushbutton','String','More Parameters','Position',pos + [0 0 0 0],'Tag','fiberGT_moreparams','Callback',{@moreparams_Callback,gcbo,[],[]});
    uicontrol('Style','pushbutton','String','More Statistics','Position',pos + [140 0 0 0],'Tag','fiberGT_morestats','Callback',{@morestats_Callback,gcbo,[],[]});

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pos = [15 410 100 18];
    uicontrol('Style','text','String','Fiber Data','Position',pos + [-10 -60 160 80],'HorizontalAlignment','left','FontWeight','bold');

    twid = [70 0 70 0];
    uicontrol('Style','text','String','FTR name','Position',pos,'Tag','fiberGT_text','HorizontalAlignment','left');
    uicontrol('Style','pushbutton','String','default_FTR','Position',pos + twid,'Tag','fiberGT_editftrname','BackGroundColor',[1 1 1],'Callback',{@setftrname_Callback,gcbo,[],[],1});
    pos = pos + [0 -25 100 0]; %%[15 395 200 18];
    twid = [150 0 -145 0];

    fibthres = params.fibrange;
    uicontrol('Style','text','String','Fiber length [min;max]','Position',pos,'Tag','fiberGT_text','HorizontalAlignment','left');
    uicontrol('Style','edit','String',sprintf('[%i;%i]',fibthres(1),fibthres(2)),'Position',pos + twid,'Tag','fiberGT_editfiblength','BackGroundColor',[1 1 1]);
    uicontrol('Style','text','String','segs','Position',pos + twid + [60 -4 -30 0],'Tag','fiberGT_text','HorizontalAlignment','left');

    pos = pos + [0 -30 -150 8];%%[15 365 50 25];
    uicontrol('Style','pushbutton','String','Load','Position',pos + [0 0 0 0],'Tag','fiberGT_savefibres','Callback',{@loadftr_Callback,gcbo,[],[],1});
    uicontrol('Style','pushbutton','String','Save','Position',pos + [50 0 0 0],'Tag','fiberGT_savefibres','Callback',{@saveftr_Callback,gcbo,[],[],1});

    uicontrol('Style','pushbutton','String','Save PM','Position',pos + [100 0 20 0],'Tag','fiberGT_saveFD','Callback',{@savefd_Callback,gcbo,[],[]});

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pos = [15 490 100 25];

    uicontrol('Style','text','String','Diffusion MR-Signal','Position',pos + [-10 -30 680 50],'HorizontalAlignment','left','FontWeight','bold');

    uicontrol('Style','text','String','<empty>','Position',pos + [310,-25,350,25],'Tag','fiberGT_info','HorizontalAlignment','left','BackGroundColor',0.8*[1 1 1]);

    uicontrol('Style','pushbutton','String','Load Mat Data','Position',pos,'Tag','fiberGT_loadmat','Callback',{@loadmatdata_Callback,gcbo,[],[]},'ToolTip','Load mrStruct with HARDI data and gradient direction information');
    uicontrol('Style','pushbutton','String','Load Nii Data','Position',pos + [100 0 0 0 ],'Tag','fiberGT_loadnii','Callback',{@loadniidata_Callback,gcbo,[],[]},'ToolTip','Estimate Whitematter mask');

    uicontrol('Style','text','String','<empty>','Position',pos + [205 0 0 0],'Tag','fiberGT_hardistat','BackGroundColor',0.8*[1 1 1],'ButtonDownFcn',@buttonsignalinfo_Callback);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    uicontrol('Style','text','String','status: ready ...','Position',[5,5,300,20],'Tag','fiberGT_status','HorizontalAlignment','left');

    uicontrol('Style','pushbutton','String','Start Tracking','Position',[280,240,100,25],'Tag','fiberGT_startstop','Callback',{@startstop_Callback,gcbo,[],[]});
    uicontrol('Style','pushbutton','String','Reset State','Position',[680,240,100,25],'Tag','fiberGT_reset','Callback',{@reset_Callback,gcbo,[],[]});

    uicontrol('Style','text','String','Tracking Log','Position',[280 425 90 20],'HorizontalAlignment','left');
    h = uicontrol('Style','listbox','String',{'Welcome to mesoFT!'},'Units','pixels','Position',[280 270 500 160],'Tag','fiberGT_log','BackGroundColor',[1 1 1],'ForeGroundColor',[1 1 1]*0);

    hax = axes;
    set(hax,'Units','pixels');
    set(hax,'position',[330 50 200 150]);
    dataStruct.axeshandles = hax;
    title('particles vs connections');
    ylabel('#Conn./#Part.');
    xlabel('iteration');
    grid on;

    hax = axes;
    set(hax,'Units','pixels','YAxisLocation','right');
    set(hax,'position',[560 50 200 150]);
    dataStruct.axeshandles(2) = hax;
    title('fiber length distribution');
    xlabel('length (mm)');
    ylabel('log(1+#fibers)');
    grid on;


    sinterpfile = fullfile(tempdir,'mesoGT',['sinterp' num2str(params.sphericalDiscNumber)]);
    try
        dataStruct.sphereInterpolation = load(sinterpfile);
    catch
        display('computing spherical interpolation LUTs');
        if isempty(dir(fullfile(tempdir,'mesoGT'))),
            mkdir(fullfile(tempdir,'mesoGT'));
        end;
        [pathtofgt dummy] = fileparts(mfilename('fullpath'));
        dwidirs = load(fullfile(pathtofgt,'dwidirections.mat'));
        sinterp = sphereInterpolLUT(getfield(dwidirs,['dirs',num2str(params.sphericalDiscNumber)])',true);
        save(sinterpfile,'-struct','sinterp');
        dataStruct.sphereInterpolation = sinterp;
    end;



    dataStruct.state = single([]);
    dataStruct.vfmap = single([]);
    dataStruct.currentiteration = 1;
    dataStruct.conratio = [];
    dataStruct.lengthhistogram = [];
    dataStruct.updownratio = [];
    dataStruct.mask_fname = '';
    dataStruct.edges = [];
    dataStruct.name = [];
    dataStruct.signal_fname = [];
    dataStruct.signal_type = [];
    dataStruct.vox = [];
    dataStruct.numit = 0;
    dataStruct.totti = 0;
    dataStruct.fixedTissueMode = params.fixedTissueMode;
    dataStruct.params = params;
    dataStruct.ftr = [];

    set(hmain,'UserData',dataStruct);

else
    figure(hmain);
end;

if ~isempty(varargin),
    busy;
end;

%%%%%%%% process parameters
%try
    k = 1;
    while k <= length(varargin),
        if strcmp(varargin{k},'loadData'),
            type = varargin{k+1};
            fn = varargin{k+2};
            gradfile = varargin{k+3};
            fnmask = varargin{k+4};
            thres = varargin{k+5};
            loadData(type,fn,gradfile,fnmask,thres);
            reportstatus('ready ...');
            k = k + 5;

        elseif strcmp(varargin{k},'setData'),
            data = varargin{k+1};
            tensor = varargin{k+2};
            mask = varargin{k+3};
            vox = varargin{k+4};
            edges = varargin{k+5};
            name = varargin{k+6};
            data.dwifile = 'manual';
            data.dwi = data;
            data.gradfile = [];
            data.tensor = tensor;
            data.edges = edges;
            data.vox = vox;
            data.name = name
            data.WM.mask = mask;
            data.WM.threshold = 0.5;
            data.WM.file = 'nofile';
            setData(data,'fromMem');
            reportstatus('ready ...');
            k = k + 6;
       elseif strcmp(varargin{k},'loadHCP'),
            hcpsubject_folder = varargin{k+1};
            k = k + 1;
            if hcpsubject_folder(end) == '/';
                hcpsubject_folder = hcpsubject_folder(1:end-1);
            end;
            [path subname] = fileparts(hcpsubject_folder);
            wmfile = fullfile(hcpsubject_folder,'T1w','wmparc.nii.gz');
            dwifile = fullfile(hcpsubject_folder,'T1w','Diffusion','data.nii.gz');
            bvecfile = fullfile(hcpsubject_folder,'T1w','Diffusion','bvecs');
            bvalfile = fullfile(hcpsubject_folder,'T1w','Diffusion','bvals');
            wm = load_nii(wmfile);
            wm.img = double((wm.img>2500 & wm.img ~= 0 )| (wm.img >= 250 & wm.img <= 255));
            wmfiletmp = fullfile(tempdir,['sub_' subname '_wmtmp.nii']);
            save_nii(wm,wmfiletmp);
            ftrname = fullfile(hcpsubject_folder,'T1w','Diffusion','the_FTR.mat');
            mesoGT_tool('setFTRname',ftrname);
            mesoGT_tool('loadData','nii',dwifile,{bvecfile bvalfile},wmfiletmp,0.5);
        elseif strcmp(varargin{k},'loadFTR'),
            fname = varargin{k+1};
            ftr  = load(fname);
            if isempty(ftr),
                errordlg(sprintf('error: %s',err),'Loading FTR');
            else
                loadftr(ftr,fname);
            end;
            k = k + 1;

        elseif strcmp(varargin{k},'saveFD'),
            saveFD([],1);

       elseif strcmp(varargin{k},'computeCFD'),
           datastruct = get(MainHandle,'UserData');
           computeCFD(datastruct);

        elseif strcmp(varargin{k},'setFTRname');
            setFTRname(varargin{k+1});
            k = k + 1;

        elseif strcmp(varargin{k},'mexit');
            datastruct = get(MainHandle,'UserData');

            Pstruc = getParams;
            if not(isfield(datastruct,'original_bTensor')),
                ready;
                errordlg('You need to load some data before mexing');
                return;
            end;
            bTensor = datastruct.original_bTensor;
            thishashnum = createDWIHashnum(datastruct,Pstruc);

            display('mexing c-code');
            mexpcRJMCMC(thishashnum);


       elseif strcmp(varargin{k},'genPhantom');
            scheme =  varargin{k+1};
            nz =  varargin{k+2};
            k = k + 2;
            genPhantom(scheme,nz);

        elseif strcmp(varargin{k},'start');
            startstop_Callback(gcf,[],[]);

       elseif strcmp(varargin{k},'start_numit');
            numits = varargin{k+1};
            k = k + 1;
            datastruct = get(MainHandle,'UserData');
            datastruct.maxnumits = numits;
            set(MainHandle,'UserData',datastruct);
            startstop_Callback(gcf,[],[]);

        elseif strcmp(varargin{k},'reset');
            reset_Callback(gcf,[],[]);


        elseif strcmp(varargin{k},'showMask');
            showMask;

        elseif strcmp(varargin{k},'showParameters');
            p = getParams;
            fields = fieldnames(p);
            for k = 1:length(fields);
                val = p.(fields{k});
                fprintf('%s : ',fields{k});
                fprintf('%f ',val);
                fprintf('\n');
            end;

      elseif strcmp(varargin{k},'reactivate');
            ready;

      elseif strcmp(varargin{k},'setparam');
            para = varargin{k+1};
            val = varargin{k+2};
            datastruct = get(MainHandle,'UserData');
            datastruct.params = setfield(datastruct.params,para,val);
            set(MainHandle,'UserData',datastruct);
            updateParamsGUI;
            k = k + 2;
        else
            fprintf('Command %s unknown!\n',['<' varargin{k} '>']);
            fprintf('Type "help fiberGT_tool" for usage information.\n');
            break;
        end;
        k = k + 1;
    end;
ready;

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    utilities
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%% returns the handle to the fiberGT figure
function h = MainHandle
    h = findobj('Tag','fiberGT_main');

%%%%% set working states
function busy
    h = findobj(MainHandle,'Style','pushbutton','-or','Style','edit');
    set(h,'Enable','off');
    set(findobj(MainHandle,'Tag','fiberGT_morestats'),'Enable','on');
    set(MainHandle,'Pointer','watch');
    h = findobj('tag','mesoGT_statfigure');
    if not(isempty(h));
        set(h,'Pointer','watch');
    end;

    drawnow;

function ready
    h = findobj(MainHandle,'Style','pushbutton','-or','Style','edit');
    set(h,'Enable','on');
    set(MainHandle,'Pointer','arrow');
    h = findobj('tag','mesoGT_statfigure');
    if not(isempty(h));
        set(h,'Pointer','arrow');
    end;

    reportstatus('ready ...');
    drawnow;

%%%%% overrides closefun of mainfig
function my_closereq(src,evnt)
    if not(ishandle(539375677)) % we are not tracking
        delete(MainHandle);
    end;

%%%%% report message in statutsbar
function reportstatus(str)
    h = findobj('Tag','fiberGT_status');
    set(h,'String',sprintf('status: %s',str));
    drawnow;

%%%%% adds entry to logwidget
function addlog(str)
    h = findobj('Tag','fiberGT_log');
    entries = get(h,'String');
    set(h,'String',{entries{:} str});
    set(h,'Value',length(entries)+1);
    printTOstderr('---------------------------------------------------------------------------------------');
    printTOstderr(str);

%%%%% returns the log
function log = getlog
    h = findobj('Tag','fiberGT_log');
    log = get(h,'String');

%%% set log content
function log = setlog(log)
    h = findobj('Tag','fiberGT_log');
    set(h,'String',log);
    set(h,'Value',length(log));

function setFTRname(name)
    hftrname = findobj('Tag','fiberGT_editftrname');
    [path fn ext] = fileparts(name);
    set(hftrname,'UserData',path);
    set(hftrname,'String',fn);

function setftrname_Callback(h,eventdata,handles,varargin)
    hftrname = findobj('Tag','fiberGT_editftrname');

    path = get(hftrname,'UserData');
    name = get(hftrname,'String');


    [newname newpath] = uiputfile(fullfile(path,[name '.mat']),'Save FTR as');

    if not(newname==0),
        [pp newnamec nnext] = fileparts(newname);
        set(hftrname,'UserData',newpath);
        set(hftrname,'String',newnamec);
    end;

function  ftrname = getftrnameFromGUI
    hftrname = findobj('Tag','fiberGT_editftrname');
    path = get(hftrname,'UserData');
    name = get(hftrname,'String');
    ftrname = fullfile(path,name);





%%%%%%%%%% infobox update
function infobox_update
    h = findobj('Tag','fiberGT_info');
    datastruct = get(MainHandle,'UserData');
    str = '';
    if isfield(datastruct,'signal'),
        sz = size(datastruct.b0avg);
        str = sprintf('Patient: %s\nDimensions: [%i %i %i], %i DE\nElement size [%.2f %.2f %.2f] mm',datastruct.name,sz,datastruct.original_size(4),datastruct.vox(1),datastruct.vox(2),datastruct.vox(3));
    end;
    set(h,'String',str);


function buttonmaskinfo_Callback(src,event)
    showMask

function showMask
        ds =  get(MainHandle,'UserData');
        if isfield(ds,'spatialProbabilities'),
            ovsamp = max(floor(size(ds.spatialProbabilities)./size(ds.b0avg)));
            dta = 3*ds.b0avg / max(ds.b0avg(:));
            if ovsamp > 1,
                dta = imresize3D_diff(dta,size(dta)*ovsamp,'nearest');
            end;
            figure;
            ARGB(:,:,:,3) = 0*ds.spatialProbabilities;
            ARGB(:,:,:,2) = ds.spatialProbabilities;
            ARGB(:,:,:,1) = dta;
            stackview(ARGB);
            title('B0-image (green) overlayed with mask (red)');
            colormap('gray');
        else
            display('No mask has been loaded/estimated!');
        end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   precomputation of LUTs and mexing the C-code
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function  model2alphaLUT = precompLUTs(bTensor,sinterp,thishashnum)

       mainhandle = MainHandle;
       datastruct = get(mainhandle,'UserData');
       Pstruc = getParams;
       alpha = Pstruc.alpha;
       ordermax = Pstruc.ordermax;

       if isempty(dir(fullfile(tempdir,'mesoGT')))
             mkdir(fullfile(tempdir,'mesoGT'));
       end;

       thedir = fullfile(tempdir,'mesoGT',hash2Str(thishashnum));

       if isempty(dir(thedir))
             mkdir(thedir);
       end;

       display(sprintf('precomputing approx. model (#dirs=%i, alpha=%f, order=%i)',size(bTensor,3),alpha,ordermax(1)));
       halpha = findobj('Tag','fiberGT_editalpha');

       lutentry = createInteractionLUTs(bTensor,sinterp,Pstruc,datastruct.fixedTissueMode,thedir);

       save(fullfile(thedir,'maLUT.mat'),'-struct','lutentry');
       model2alphaLUT = lutentry.model2alphaLUT;

       display('mexing c-code');
       mexpcRJMCMC(thishashnum);


function mexpcRJMCMC(hashnum)

       params = getParams;
       compilerparams = params.compilerparams;

       if params.numcores > 1,
           compilerparams = ['-DPARALLEL_OPENMP ' compilerparams];
       end;

       if not(ispc),
          compilerparams = ['-DLINUX_MACHINE ' compilerparams];
       end;

       libs = params.libs;

       thedir = fullfile(tempdir,'mesoGT',hash2Str(hashnum));
       [pathtofgt dummy] = fileparts(mfilename('fullpath'));
       mexstr = sprintf('mex CXXFLAGS="\\$CXXFLAGS %s" %s -outdir %s -output %s -I%s %s',compilerparams,fullfile(pathtofgt,'ccode','pcRJMCMC.cpp') , thedir,fullfile(thedir,'pcRJMCMC'),thedir,libs);
       fprintf('%s\n',mexstr);
       eval(mexstr);



function str = hash2Str(hashnum)
    X = 'qwertzuiopasdfghjklyxcvbnmQWERTZUIOPASDFGHJKLYXCVBNM123456789';
    hashnum = abs(hashnum);
    str = [];
    for k = 1:length(hashnum),
        str = [str  X(1+fix(rem(hashnum(k)*length(X).^(-1:6),length(X))))];
    end;
    str = strrep(str,'q','');
    str = str(1:min(length(str),52));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    loading the DWI data and mask
%    and Preprocessing the data
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% callbacks for loading data
function varargout = loadmatdata_Callback(h,eventdata,handles,varargin)
   loadData('mat','');
function varargout = loadniidata_Callback(h,eventdata,handles,varargin)
   loadData('nii','');

%%%%%% load data (called from the user or by loadFTR)
function loadData(type,fn,gradinfo,fnmask,threshold)
        reportstatus('reading DWI data');
        busy;
        if isempty(fn),
            data = eval(['loadData_' type ';']);
        else
            data = eval(['loadData_' type '(fn,gradinfo,fnmask,threshold);']);
        end;
        if not(isempty(data))
            setData(data,type);
        end;
        ready;

%%%%% sets the loaded data and does the preprocessing
function setData(data,type)

        datastruct = get(MainHandle,'UserData');
        datastruct.edges = data.edges;
        datastruct.spatialProbabilities = single(data.WM.mask);
        datastruct.mask_fname =  data.WM.file;
        datastruct.mask_threshold = data.WM.threshold;
        datastruct.signal_fname = data.dwifile;
        datastruct.signal_type = type;
        datastruct.grad_fname = data.gradfile;
        datastruct.edges = data.edges;
        datastruct.vox = data.vox(1:3);
        datastruct.name = data.name;
        set(MainHandle,'UserData',datastruct);

        Pstruc = getParams;
        if not(isempty(Pstruc.maxbval))
            bval = data.tensor(1,1,:)+data.tensor(2,2,:)+data.tensor(3,3,:);
            idx = find(bval<Pstruc.maxbval+100);
            data.tensor = data.tensor(:,:,idx);
            data.dwi = data.dwi(:,:,:,idx);
        end;

        preprocessData(data.tensor,data.dwi,datastruct.spatialProbabilities);

        infobox_update ;
        if iscell(data.dwifile)
            filenamedwi = data.dwifile{1};
        else
            filenamedwi = data.dwifile;
        end;
        addlog(sprintf('DWI file %s loaded',filenamedwi));
        hstat = findobj('Tag','fiberGT_hardistat');
        set(hstat,'String','DWI loaded');
        [path fnamec] = fileparts(filenamedwi);
        setFTRname(fullfile(path,[strrep(fnamec,'_HARDI','') '_FTR']));
        h = findobj('tag','mesoGT_statfigure');
        if not(isempty(h)),
            close(h);
        end;
        ready;



%%%%% does the actual preprocessing
function preprocessData(bTensor,signal,mask)

       bTensor = bTensor/1000;

       datastruct = get(MainHandle,'UserData');

       reportstatus('computing b0 average');
       bval = squeeze(bTensor(1,1,:)+bTensor(2,2,:)+bTensor(3,3,:));
       b0idx = find(bval <= 0.1);
       datastruct.b0avg = single(mean(signal(:,:,:,b0idx),4));
       datastruct.original_bTensor = bTensor;
       datastruct.original_size = size(signal);


       % if mask has double the resolution of the data
       if any(size(mask(:,:,:,1)) ~= size(datastruct.b0avg)),
           osamp = round(size(mask,1)/size(datastruct.b0avg,1));
           ker = ones(osamp,osamp,osamp); ker = ker /sum(ker(:));
           mask = imfilter(double(mask),ker);
           st = ceil(osamp/2);
           mask = mask(st:osamp:end,st:osamp:end,st:osamp:end,1);
       end;

       mask = mask(:,:,:,1)>0;

       if datastruct.params.directional_propsal_distrib,
           for k = 1:size(bTensor,3),
               [U D] = eigs(bTensor(:,:,k));
               bDir(:,k) = U(:,1)*sqrt(D(1,1)*1000);
           end;

           ds = GLQball_DSI(signal,bDir,'lambda',5,'Nmax',10,'shell',4.5,'sinterp',datastruct.sphereInterpolation,'D',2);
           ds.signal = ds.signal.^2.*(ds.signal>0);
           ds.signal(:,:) = ds.signal(:,:) ./ (eps+repmat(sum(ds.signal(:,:)),[size(ds.signal,1) 1]));

           figure(98234);
           stackviewODF(ds,[0 0 1]);
       end

       reportstatus('normalizing by b0');
       signal = permute(signal,[4 1 2 3]);
       signal = signal(:,mask(:));
       for k = 1:size(signal,1),
           signal(k,:) = squeeze(signal(k,:))./ (eps+datastruct.b0avg(mask(:)))';
       end;
       signal(signal>1.5) = 1.5;
       signal = single(signal);

       reportstatus('preparing approximation model');
       Pstruc = getParams;
       alpha = Pstruc.alpha;
       ordermax = Pstruc.ordermax;
       thishashnum = createDWIHashnum(datastruct,Pstruc);
       try
            thedir = fullfile(tempdir,'mesoGT',hash2Str(thishashnum));
            mmluts = load(fullfile(thedir,'maLUT.mat'));
            model2alphaLUT = mmluts.model2alphaLUT;
       catch
           mmluts = [];
           hashnum = zeros(1,18);
       end
       if isempty(mmluts),
           model2alphaLUT = precompLUTs(bTensor,datastruct.sphereInterpolation,thishashnum);
       end;
       [W] = createWeightingScheme(bTensor,Pstruc.b_weighting);

       reportstatus('computing approximation coefficients');
       model2alphaLUT = single( model2alphaLUT);
       datastruct.meansignal = sum(W*signal(:,:));


       % the data
       datastruct.signal = reshape(reshape(model2alphaLUT,[size(model2alphaLUT,1)*size(model2alphaLUT,2) size(model2alphaLUT,3)]) * signal(:,:), ...
                                       [size(model2alphaLUT,1) size(model2alphaLUT,2) size(signal,2)]);


       if datastruct.params.directional_propsal_distrib,
           datastruct.signal(end+1,:,:) = ds.signal(:,mask(:));
       end;

       % S2(:,:,:,1) contains the q-space mean
       % S2(:,:,:,2) contains the idxmap used for sparse rep of signal
       datastruct.S2 = zeros(size(datastruct.b0avg),'single');
       datastruct.S2(mask) = datastruct.meansignal;
       idxmap = zeros(size(datastruct.b0avg),'single');
       idxmap(mask) = 1:sum(mask(:));
       datastruct.S2(:,:,:,2) = idxmap;


       set(MainHandle,'UserData',datastruct);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    load/save FTR
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% collect data for inclusion in FTR
function info = gatherFTRinfo(datastruct)
    info.fiberGT_FTRversion = 'V1.0';
    info.params = getParams;
    info.Log = getlog;
    info.edges = datastruct.edges;
    info.name = datastruct.name;
    info.signal_fname = datastruct.signal_fname;
    info.signal_type = datastruct.signal_type;
    info.mask_fname = datastruct.mask_fname;
    info.mask_threshold = datastruct.mask_threshold;
    info.grad_fname = datastruct.grad_fname;
    info.conratio = datastruct.conratio;
    info.lengthhistogram = datastruct.lengthhistogram;
    info.updownratio = datastruct.updownratio;
    info.currentiteration = datastruct.currentiteration;
    info.totti = datastruct.totti;


%%%%% callback for saving a ftr
function varargout = saveftr_Callback(h,eventdata,handles,varargin)
    busy;

    datastruct = get(MainHandle,'UserData');
    ftrname = getftrnameFromGUI;
    info = gatherFTRinfo(datastruct);
    fiblen = info.params.fibrange;

    reportstatus('building fibers');
    fibidx = BuildFibres(datastruct.state,fiblen);

    reportstatus('saving fibers');
    [ftr lens] = mesosaveFTR(ftrname,datastruct,fibidx,datastruct.vox,info,true);
    datastruct.lengthhistogram = lenhist(lens);
    updatePlots(datastruct);
    set(MainHandle,'UserData',datastruct);
    addlog(sprintf('%i fibers saved',length(fibidx)));
    ready;


%%%%% callback for loading a ftr
function varargout = loadftr_Callback(h,eventdata,handles,varargin)
    [files path] = uigetfile({'*_FTR.mat', 'All accepted filetypes';'*.*','All files'},'Select FTR to open');
    ftr = load(fullfile(path,files));
    if isempty(ftr),
        ready;
        return;
    end;
    loadftr(ftr,fullfile(path,files));

%%% % the actual loading procedure
function varargout = loadftr(ftr,fname)

    busy;

    wrongversion = 0;
    if isfield(ftr.trackParam,'fiberGT_FTRversion'),
        if ~strcmp(ftr.trackParam.fiberGT_FTRversion,'V1.0'),
            wrongversion = 1;
        end;
    else
        wrongversion = 1;
    end;

    if wrongversion == 1,
        errordlg('Not a valid mesoFT-FTR file!','Loading FTR');
        ready;
        return;
    end;


    type = ftr.trackParam.signal_type;
    signal_fname = ftr.trackParam.signal_fname;
    mask_fname = ftr.trackParam.mask_fname;
    mask_threshold = ftr.trackParam.mask_threshold;
    grad_fname = ftr.trackParam.grad_fname;


    if ~fileexists(mask_fname) | ~fileexists(signal_fname),
        uiwait(warndlg('Associated data do not exist!','Loading FTR','modal'));
        signal_fname = [];
        mask_fname = [];
    end;

    loadData(type,signal_fname,grad_fname,mask_fname,mask_threshold);

    setFTRname(fname);
    datastruct = get(MainHandle,'UserData');

    datastruct.conratio = ftr.trackParam.conratio;
    datastruct.lengthhistogram = ftr.trackParam.lengthhistogram;
    datastruct.state = ftr.user.P;
    datastruct.vfmap = ftr.user.vf;
    datastruct.currentiteration = ftr.trackParam.currentiteration;
    datastruct.totti = ftr.trackParam.totti;
    datastruct.ftr = ftr;
    if isfield(ftr.trackParam,'updownratio'),
        datastruct.updownratio = ftr.trackParam.updownratio;
    else,
        datastruct.updownratio = zeros(size(datastruct.conratio));
    end;
    setlog(ftr.trackParam.Log);
    datastruct.params = ftr.trackParam.params;
    set(MainHandle,'UserData',datastruct);

    updateParamsGUI;
    updatePlots(datastruct);

    ready;

%%%% callpack for saving paramteric maps
function varargout = savefd_Callback(h,eventdata,handles,varargin)
    saveFD([],1);

%%%% computing the maps and saving them as niftis
function saveFD(sel,M)
    busy;
    reportstatus('starting to build parametric maps');

    datastruct = get(MainHandle,'UserData');
    ftrname = getftrnameFromGUI;

    if ~isfield(datastruct,'state'),
        errordlg('No proper Tracking State!','Save Fiber Densities');
        ready;
        return;
    end

    if isempty(datastruct.state),
        errordlg('No fiber densities saved!','Save Fiber Densities');
        ready;
        return;
    end;

    info = gatherFTRinfo(datastruct);

    reportstatus('building fibers');
    fibidx = BuildFibres(datastruct.state,info.params.fibrange);
    ftr = mesosaveFTR(ftrname,datastruct,fibidx,datastruct.vox,info,false);

    reportstatus('building fiber densities');
    sz = size(datastruct.spatialProbabilities);
    [fpath fname ext] = fileparts(ftrname);
    maps = mesoftr2FDmaps(ftr);

    reportstatus('saving fiber densities');

    contrast = {'Daxon','Dextrapara','Dextraorth','turt','vfi','vfe','vfsw','fd','segcount','axondiameter'};

    for j = 1:length(contrast),
        mr = maps.mrProp;
        mr.dataAy = getfield(maps,contrast{j});
        save_mrstruct_as_nifti(mr,fullfile(fpath,[fname '_'  contrast{j} '.nii']));
        %mrstruct_to_nifti(mr,fullfile(fpath,[fname '_'  contrast{j} '.nii']));
    end;

    ready;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Parameter stuff
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%% get parametre values from GUI
function updateParams
    mainhandle = MainHandle;
    datastruct = get(mainhandle,'UserData');
    p = datastruct.params;
    fnames = fieldnames(p);
    for k = 1:length(fnames),
        h = findobj(gcf,'tag',['fiberGT_edit_' fnames{k}]);
        if not(isempty(h)),
            val = str2num(get(h,'String'));
            p = setfield(p, fnames{k},val);
        end;
    end;
    hfl = findobj('Tag','fiberGT_editfiblength');
    p.fibrange = eval(['[' get(hfl,'String') ']'])';
    datastruct.params = p;
    set(mainhandle,'UserData',datastruct);

function updateParamsGUI
    mainhandle = MainHandle;
    datastruct = get(mainhandle,'UserData');
    p = datastruct.params;
    fnames = fieldnames(p);
    for k = 1:length(fnames),
        h = findobj(gcf,'tag',['fiberGT_edit_' fnames{k}]);
        if not(isempty(h)),
            set(h,'String',num2str(getfield(p, fnames{k})));
        end;
    end;
    hfl = findobj('Tag','fiberGT_editfiblength');
    set(hfl,'String',sprintf('[%i;%i]',p.fibrange(1),p.fibrange(2)));



%%%%% gets parameter struct
function p = getParams
    updateParams;
    mainhandle = MainHandle;
    datastruct = get(mainhandle,'UserData');
    p = datastruct.params;


function varargout = moreparams_Callback(h,eventdata,handles,varargin)
    updateParams;
    mainhandle = MainHandle;
    datastruct = get(mainhandle,'UserData');
    datastruct.params  = editParamStruct(datastruct.params);
    set(mainhandle,'UserData',datastruct);
    updateParamsGUI;

function varargout = morestats_Callback(h,eventdata,handles,varargin)
    h = findobj('tag','mesoGT_statfigure');
    mainhandle = MainHandle;
    datastruct = get(mainhandle,'UserData');
    if isempty(h),
        ssds.contrast = [1 3 5 6];
        ssds.fibcoloring = 1;
        h = figure('Name','Tracking Stats','tag','mesoGT_statfigure','NumberTitle','off', 'userdata',ssds);
        set(h,'Toolbar','figure');set(h,'Menubar','none');
        pos = get(h,'position');
        sz = [1000 800];
        pos = [(pos(1:2)-(sz-pos(3:4))) sz];
        set(h,'position',pos);

        startstopstate = get(findobj('Tag','fiberGT_startstop'),'Enable');
        busy;
        showTractStats(datastruct);
        set(findobj('Tag','fiberGT_startstop'),'Enable',startstopstate);
        if not(ishandle(539375677)) % we are not tracking
            ready;
        end;
    end;


%%%%% checks data before tracking is started
function ret = checkdata(datastruct)

    ret = false;
    if ~isfield(datastruct,'spatialProbabilities'),
        errordlg('No White Matter Mask loaded','Data consistency');
        return
    end;

    if ~isfield(datastruct,'signal'),
        errordlg('No signal loaded','Data consistency');
        return
    end;

    if ~isfield(datastruct,'sphereInterpolation');
        errordlg('Sphere Interpolator not initialized','Data consistency');
        return
    end;


    ret = true;
    return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Tracking
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% callback for start tracking
function varargout = startstop_Callback(h,eventdata,handles,varargin)

    hbut = findobj('tag','fiberGT_startstop');
    if not(isempty(get(hbut,'userdata'))),
        hsb = get(hbut,'userdata');
        if ishandle(hsb),
            delete(hsb);
        end;
        set(hbut,'userdata',[]);
        set(hbut,'String','Start Tracking');
        return;
    end;

    busy;
    updateParams;

    datastruct = get(MainHandle,'UserData');
    if ~checkdata(datastruct),
        ready; return;
    end;

    Pstruc = getParams;

    szvf = size(datastruct.vfmap);
    szsig = size(datastruct.S2);
    if szvf(1) > 0,
        if any(szvf(1:3)/Pstruc.p_wid ~= szsig(1:3)),
            display('oversampling changed, please reset the state!');
            ready; return;
        end;
    end;

    % check whether LUTs and exectubale is there
    thishashnum = createDWIHashnum(datastruct,Pstruc);
    thedir = fullfile(tempdir,'mesoGT',hash2Str(thishashnum));
    mexfilename = 'pcRJMCMC';
    if isempty(dir(fullfile(thedir,[mexfilename '*']))),
        precompLUTs(datastruct.original_bTensor,datastruct.sphereInterpolation,thishashnum);
    end;


    % do it
    reportstatus('tracking');
    addlog('Tracking started');

    info.params = Pstruc;
    info.name = datastruct.name;
    datastruct = trackit(datastruct,Pstruc.Tstart,Pstruc.Tend,Pstruc.numstep,Pstruc.numiterations,info);
    set(MainHandle,'UserData',datastruct);
    ready;
    addlog('Tracking stopped');
    reportstatus('ready ...');



%%%%% tracking
function datastruct = trackit(datastruct,Tstart,Tend,ni,its,info)

    if its > 0,
        ipr = round(its/ni);
        lam_DE = 0;
    else
        ipr = 10^9;
        lam_DE = 1-10^its;
    end;

    lens = [];
    h = mesocreateStopButton;
    set(h,'visible','off');
    hb = findobj('tag','fiberGT_startstop');
    set(hb,'String','Stop Tracking');
    set(hb,'userdata',h);
    set(hb,'enable','on');


    sz = size(datastruct.spatialProbabilities);
    alpha = log(Tend/Tstart);
    vox = datastruct.vox;


    thishashnum = createDWIHashnum(datastruct,info.params);
    mexfilename = ['pcRJMCMC'];


    ftrname = getftrnameFromGUI;
    fiblen =  info.params.fibrange;
    info = gatherFTRinfo(datastruct);



    thedir = fullfile(tempdir,'mesoGT',hash2Str(thishashnum));
    mask = datastruct.spatialProbabilities;


    doexit = false;
    cnt = 1;


    ck = datastruct.currentiteration;


    while(1),

          datastruct.currentiteration = ck;

          tic

          T = Tstart * exp(alpha*(ck-1)/ni);
          paramStruct = datastruct.params;
          paramStruct.currentTemp = T;
          paramStruct.maxit = ipr;
          paramStruct.lam_DE = lam_DE;
          paramStruct.num_qspace_samples = size(datastruct.original_bTensor,3);
          paramStruct.GM_mask = double(size(datastruct.spatialProbabilities,4)>1);


          curdir = cd(thedir);
         % xx = dir([thedir '/pcRJMCMC.*']); fprintf('executable last change: %s\n',xx.date);
          eval(sprintf('[datastruct.state datastruct.vfmap stat clist] = %s(datastruct.state,datastruct.vfmap,datastruct.signal,mask,datastruct.S2,double(vox),paramStruct,datastruct.sphereInterpolation,h);',mexfilename));
          cd(curdir);
          neededtime = toc;

          concnt = (sum(datastruct.state(9,:) ~= -1) + sum(datastruct.state(10,:) ~= -1))/2;
          segcnt = size(datastruct.state,2);
          datastruct.numit = datastruct.numit + stat.iterations;
          datastruct.totti = datastruct.totti + neededtime;
          datastruct.clist = clist;

          fibidx = BuildFibres(datastruct.state,fiblen);

          addlog(sprintf('%i) p:%i #c:%i #f:%i T:%.0e ac:%.1f%% it:%.0e/%.0e ti:%.2fm/%.2fh \n',ck,segcnt,concnt,length(fibidx), T, 100*stat.accepted/(stat.accepted+stat.rejected), stat.iterations,  round(datastruct.numit), neededtime/60,datastruct.totti/3600));


          info = gatherFTRinfo(datastruct);
          [ftr lens] = mesosaveFTR(ftrname,datastruct,fibidx,datastruct.vox,info);

          datastruct.conratio(ck) = concnt/segcnt;
          datastruct.updownratio(ck) = stat.up/(stat.down+1);
          datastruct.lengthhistogram = lenhist(lens);
          datastruct.ftr = ftr;
          updatePlots(datastruct);

          hsf = findobj('tag','mesoGT_statfigure');
          if not(isempty(hsf));
              showTractStats(datastruct);
          end;



          if isfield(datastruct,'maxnumits'),
              if cnt >= datastruct.maxnumits,
                  datastruct.currentiteration = ck+1;
                  doexit = true;
              end;
          end;
          cnt = cnt + 1;


         if mesostopButtonPressed,
              addlog('interrupted by user');
              doexit = true;
          end;

          if doexit | ck >= ni
              break;
          end;

          ck = ck + 1;

    end;



    if ~mesostopButtonPressed,
        hbut = findobj('tag','fiberGT_startstop');
        if not(isempty(get(hbut,'userdata'))),
            hsb = get(hbut,'userdata');
            if ishandle(hsb),
                delete(hsb);
            end;
            set(hbut,'userdata',[]);
            set(hbut,'String','Start Tracking');
            return;
        end;
    end;




%%%%% clears tracking state
function varargout = reset_Callback(h,eventdata,handles,varargin)
    datastruct = get(MainHandle,'UserData');
    datastruct.state = single([]);
    datastruct.vfmap = single([]);
    datastruct.ftr = [];
    datastruct.currentiteration = 1;
    datastruct.numit = 0;
    datastruct.totti = 0;
    datastruct.conratio = [];
    datastruct.updownratio = [];
    cla(datastruct.axeshandles(1));
    cla(datastruct.axeshandles(2));
    set(MainHandle,'UserData',datastruct);
    addlog('Tracks deleted, state reseted');



function h = lenhist(lens)
h = [];
if ~isempty(lens),
        hi = histc(lens,1:1:max(lens));
        h = log(1+hi);
end;




%%%%% update tracking plots
function updatePlots(datastruct)

     if length(datastruct.conratio) > 1,
         n = length(datastruct.conratio);
        plot(datastruct.axeshandles(1),1:n,datastruct.conratio,'b'); %,1:n,datastruct.updownratio,'r');

        axis(datastruct.axeshandles(1),[1 length(datastruct.conratio) 0 1]);
        grid(datastruct.axeshandles(1),'on');
        title(datastruct.axeshandles(1),'particles vs connections');
        ylabel(datastruct.axeshandles(1),'#Conn./#Part.');
        xlabel(datastruct.axeshandles(1),'iteration');
      end;


      if ~isempty(datastruct.lengthhistogram),
        bar(datastruct.axeshandles(2),datastruct.lengthhistogram);
        grid(datastruct.axeshandles(2),'on');
        title(datastruct.axeshandles(2),'fiber length distribution');
        xlabel(datastruct.axeshandles(2),'length (mm)');
        ylabel(datastruct.axeshandles(2),'log(1+#fibers)');
        set(datastruct.axeshandles(2),'YAxisLocation','right');
      end;




function genPhantom(bTensor,nz)

           [dataAy themask] = genPhan(bTensor/1000,nz);
           data.dwifile = 'phantom';
           data.dwi = dataAy;
           data.gradfile = [];
           data.tensor = bTensor;
           data.edges = eye(4);
           data.vox = [2 2 2];
           data.name = 'phantom';
           data.WM.mask = themask;
           data.WM.threshold = 0.5;
           data.WM.file = 'nofile';
           setData(data,'fromMem');


















