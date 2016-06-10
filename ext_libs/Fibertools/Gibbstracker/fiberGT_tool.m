function ret = fiberGT_tool(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % fiberGT_tool(commandstr1,varargin1,commandstr2,varargin2,.....)
% %
% %  implements global fibertracking 
% %  based on M.Reisert et al. 'Global Reconstruction of Neuronal Fibers',
% %           MICCAI, Workshop DMFC 2009, pp 70-82
% %              and
% %           M.Reisert et al. 'Global Fiber Reconstruction becomes Practical',
% %                NeuroImage 2010
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
% %        least 2 Gigabyte of free memory (advise: crop your data).
% %
% % -----------------------------------------------------------------
% % commandstr == 'setHARDI'
% %     fiberGT_tool('setHARDI',bTensor,Signal,Voxelsize,name)
% %     loads hardi-data to the internal datastructure of the tracker.
% %      - bTensor   , direction info of size [3 3 N] (alternativly [3 N])      
% %      - Signal    , HARDI signal of size [m n k N]
% %      - Voxelsize , edge length of voxels in (millimeter) as [xlen ylen zlen]
% %      - name      , string with e.g. patient name
% %
% % commandstr == 'setDSI'
% %     fiberGT_tool('setDSI',bTensor,Signal,Voxelsize,name)
% %     loads dsi-data to the internal datastructure of the tracker.
% %      - bTensor   , direction info of size [3 3 N], contains the
% %                    gradient directions as a tensor, while the
% %                    eigenvalue corresponds to the b-value in s/mm²  (i.e. n*n' * b)
% %      - Signal    , HARDI signal of size [m n k N]
% %      - Voxelsize , edge length of voxels in (millimeter) as [xlen ylen zlen]
% %      - name      , string with e.g. patient name
% %
% % commandstr == 'loadHARDI'
% %     fiberGT_tool('loadHARDI',filename)
% %     loads a HARDI-mrstruct.
% % -----------------------------------------------------------------
% % commandstr == 'setDTD'
% %     fiberGT_tool('setDTD',EigStr,Voxelsize,name)
% %     loads diffusiontensor-data to the internal datastructure of the tracker.
% %      - EigStr    , structure containing EigStr.Eval of size [m n k 3] and 
% %                                         EigStr.EV1, EigStr.EV2, EigStr.EV3 containing the 
% %                                         corresponding eigenvectors, each of size [m n k 3]
% %      - bval      , b-value of measurement
% %      - Voxelsize , edge length of voxels in (millimeter) as [xlen ylen zlen]
% %      - name      , string with e.g. patient name
% %
% % commandstr == 'loadDTD'
% %     fiberGT_tool('loadDTD',filename,bvalue)
% %     loads a dtdstruct.
% % ----------------------------------------------------------------
% % commandstr == 'setMask'
% %     fiberGT_tool('setMask',mask)
% %      - mask  , White matter mask of size [m n k] or 'all' if the whole
% %                area should be taken
% %
% % commandstr == 'loadMask'
% %     fiberGT_tool('loadMask',filename,threshold)
% %     loads a mrstruct or a nifti as white matter mask
% %     (if a nifti is given a HARDI must be loaded before)
% % 
% % commandstr == 'loadMaskStruct'
% %     fiberGT_tool('loadMaskStruct',filename,name)
% %     loads the ROI <name> from a maskstruct as white matter mask
% %
% % commandstr == 'estimateMask'
% %      estimates white matter mask (more or less a brain mask).
% %      This function only works well for 'usual' whole brain images,
% %      because global histograms are involved to determine thresholds.
% %
% % commandstr == 'showMask'
% %     shows the white matter mask
% % ----------------------------------------------------------------
% % commandstr == 'loadFTR'
% %     fiberGT_tool('setFTR',filename)
% %      - filename  , filename of ftr-file (the fibertracks) formerly created with fiberGT_tool
% %
% % commandstr == 'send'
% %     send fibers to fiberviewer
% %
% % commandstr == 'setFTRname'
% %     fiberGT_tool('setFTR',filename)
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
% % commandstr == 'createCmat'
% %     cmat = fiberGT_tool('createCmat',roi_label_img,num_samples,num_its,Temp,Nsz,noavg)
% %
% %     computes a connectivity matrix (CM) accordings to the ROIs given in the
% %     labelimg <roi_label_img> (size [m n k]). The CM is averaged
% %     <num_samples> times, where between each sample <num_its> iterations are
% %     performed at a Temperature <Temp>. To compute the overlap between
% %     the fiber terminals and the ROIs the fiber terminal is assumed to have width <Nsz>
% %     cmat.cc contains the CM and cmat.ll the average fiber lengths
% %     if noavg=='noavg' then the indivdual cmat-samples are returned and
% %     no avgeraging is performed
% % ----------------------------------------------------------------
% % commandstr == 'suggestSparse'
% % commandstr == 'suggestDense'
% %     suggest parameter settings. Works well for 'usual' whole 
% %     brain images.
% %
% % commandstr == 'setparam'
% %     fiberGT_tool('setparam',parametername,value);
% %       possible parameternames:
% %       'width','length','weight','connlike','bfac','startTemp',
% %       'stopTemp','steps','itnum','ftrname','fiblength','chemPot1',
% %       'chemPot2',
% %       'inexbal' (-5..5 = in...ex emph.),
% %       'fiblength' (string of the form 'minlen;maxlen')
% %     
% % commandstr == 'setVoxelsize'
% %     set voxelsize in (millimeter) as [xlen ylen zlen]
% %
% % commandstr == 'showSignalHist'
% %     shows an intensity histogram of the FRT-signal
% %
% % commandstr == 'showParameters'
% %     show current parameters
% %
% % commandstr == 'smoothODF'
% %    fiberGT_tool('smoothODF',numit)
% % 
% % commandstr == 'getDS'
% %    return internal datastruct
% %
% % commandstr == 'reactivate'
% %    reactive uicontrols after crash
% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  A simple example for loading some hardi-data and start the tracking:
% %  Example: > fiberGT_tool('setHARDI',bTensor,Signal,[2 2 2],'Donald')
% %           > fiberGT_tool('estimateMask')
% %           > fiberGT_tool('suggestSparse')
% %           > fiberGT_tool('start')
% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ret = [];

hmain = findobj('Tag','fiberGT_main');
if isempty(hmain),
    hmain = figure('Name','Global Fiber Tracking','Tag','fiberGT_main','MenuBar','none','NumberTitle','off', ...
                   'Position',[0 0 800 550],'resize','off','CloseRequestFcn',@my_closereq); 
    
               
               
               
    %%%%%%%%%%%%%%%%%%%%%
                                                            
    uicontrol('Style','edit','String','0','Position',[0 0 100 20],'Tag','fiberGT_editchemPot1','BackGroundColor',[1 1 1],'Visible','off');
    uicontrol('Style','edit','String','0.5','Position',[0 0 100 20],'Tag','fiberGT_editconnlike','BackGroundColor',[1 1 1],'Visible','off');                         
    uicontrol('Style','edit','String','0','Position',[0 0 100 20],'Tag','fiberGT_editinexbal','BackGroundColor',[1 1 1],'Visible','off');                         
               
    %%%%%%%%%%%%%%%%%%%%%%%%%%               
    pos = [5 300 130 20];
    twid = [140 0 -50 0];
    thei = -[0 25 0 0];
    uicontrol('Style','text','String','Iteration Parameters','Position',pos + [0 -80 130 100],'HorizontalAlignment','left','FontWeight','bold');
    
    uicontrol('Style','text','String','start Temp.','Position',pos,'Tag','fiberGT_text');
    uicontrol('Style','edit','String','0.1','Position',pos + twid,'Tag','fiberGT_editstartTemp','BackGroundColor',[1 1 1]);
    pos = pos + thei;
    uicontrol('Style','text','String','stop Temp.','Position',pos,'Tag','fiberGT_text');
    uicontrol('Style','edit','String','0.001','Position',pos + twid,'Tag','fiberGT_editstopTemp','BackGroundColor',[1 1 1]);
    pos = pos + thei;
    uicontrol('Style','text','String','# Steps','Position',pos,'Tag','fiberGT_text');
    uicontrol('Style','edit','String','50','Position',pos + twid,'Tag','fiberGT_editsteps','BackGroundColor',[1 1 1]);
    pos = pos + thei;
    uicontrol('Style','text','String','# Iterations','Position',pos,'Tag','fiberGT_text');
    uicontrol('Style','edit','String','5*10^8','Position',pos + twid,'Tag','fiberGT_edititnum','BackGroundColor',[1 1 1]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pos = [5 170 130 20];
    twid = [140 0 -50 0];
    thei = -[0 25 0 0];
    uicontrol('Style','text','String','Cylinder Parameters','Position',pos + [0 -105 130 125],'HorizontalAlignment','left','FontWeight','bold');
    
    uicontrol('Style','text','String','Width (mm)','Position',pos,'Tag','fiberGT_text');
    uicontrol('Style','edit','String','???','Position',pos + twid,'Tag','fiberGT_editwidth','BackGroundColor',[1 1 1]);
    pos = pos + thei;
    uicontrol('Style','text','String','Length (mm)','Position',pos,'Tag','fiberGT_text');
    uicontrol('Style','edit','String','???','Position',pos + twid,'Tag','fiberGT_editlength','BackGroundColor',[1 1 1]);
    pos = pos + thei;
    uicontrol('Style','text','String','Weight (unitless)','Position',pos,'Tag','fiberGT_text');
    uicontrol('Style','edit','String','???','Position',pos + twid,'Tag','fiberGT_editweight','BackGroundColor',[1 1 1]);
    pos = pos + thei;
    uicontrol('Style','text','String','Dens. Penalty','Position',pos,'Tag','fiberGT_text');
    uicontrol('Style','edit','String','0.2','Position',pos + twid,'Tag','fiberGT_editchemPot2','BackGroundColor',[1 1 1]);
    pos = pos + thei;
    uicontrol('Style','text','String','b-value','Position',pos,'Tag','fiberGT_text');
    uicontrol('Style','edit','String','1.0','Position',pos + twid,'Tag','fiberGT_editbfac','BackGroundColor',[1 1 1]);
    pos = pos + thei;
    uicontrol('Style','text','String','Suggest Parameters','Position',pos + [0 -5 0 0],'Tag','fiberGT_text');
    uicontrol('Style','pushbutton','String','Dense','Position',pos + [140 -5 -70 0],'Tag','fiberGT_dense','Callback',{@suggest_Callback,gcbo,[],[],1});
    uicontrol('Style','pushbutton','String','Sparse','Position',pos + [200 -5 -70 0],'Tag','fiberGT_sparse','Callback',{@suggest_Callback,gcbo,[],[],2});
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pos = [15 410 100 18];
    uicontrol('Style','text','String','Fiber Data','Position',pos + [-10 -60 160 80],'HorizontalAlignment','left','FontWeight','bold');
 
    twid = [70 0 70 0];
    uicontrol('Style','text','String','FTR name','Position',pos,'Tag','fiberGT_text','HorizontalAlignment','left');
    uicontrol('Style','pushbutton','String','default_FTR','Position',pos + twid,'Tag','fiberGT_editftrname','BackGroundColor',[1 1 1],'Callback',{@setftrname_Callback,gcbo,[],[],1});
    pos = pos + [0 -25 100 0]; %%[15 395 200 18];
    twid = [150 0 -145 0];
    uicontrol('Style','text','String','Fiber length [min;max]','Position',pos,'Tag','fiberGT_text','HorizontalAlignment','left');
    uicontrol('Style','edit','String','10;inf','Position',pos + twid,'Tag','fiberGT_editfiblength','BackGroundColor',[1 1 1]);
    uicontrol('Style','text','String','segs','Position',pos + twid + [60 -4 -30 0],'Tag','fiberGT_text','HorizontalAlignment','left');

    pos = pos + [0 -30 -150 8];%%[15 365 50 25];
    uicontrol('Style','pushbutton','String','Load','Position',pos + [0 0 0 0],'Tag','fiberGT_savefibres','Callback',{@loadftr_Callback,gcbo,[],[],1});
    uicontrol('Style','pushbutton','String','Save','Position',pos + [50 0 0 0],'Tag','fiberGT_savefibres','Callback',{@saveftr_Callback,gcbo,[],[],1});
    uicontrol('Style','pushbutton','String','Send to FV','Position',pos + [100 0 30 0],'Tag','fiberGT_send','Callback',{@sendto_Callback,gcbo,[],[]});   
    uicontrol('Style','pushbutton','String','Save FD','Position',pos + [180 0 10 0],'Tag','fiberGT_saveFD','Callback',{@savefd_Callback,gcbo,[],[]});   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pos = [15 490 100 25];
  
    uicontrol('Style','text','String','Diffusion MR-Signal','Position',pos + [-10 -30 680 50],'HorizontalAlignment','left','FontWeight','bold');
    
    uicontrol('Style','text','String','<empty>','Position',pos + [310,-25,350,25],'Tag','fiberGT_info','HorizontalAlignment','left','BackGroundColor',0.8*[1 1 1]);
 
    uicontrol('Style','pushbutton','String','Load HARDI','Position',pos,'Tag','fiberGT_loadhardi','Callback',{@loadhardi_Callback,gcbo,[],[]},'ToolTip','Load mrStruct with HARDI data and gradient direction information');  
    uicontrol('Style','text','String','<empty>','Position',pos + [205 0 0 0],'Tag','fiberGT_hardistat','BackGroundColor',0.8*[1 1 1],'ButtonDownFcn',@buttonsignalinfo_Callback);
    uicontrol('Style','pushbutton','String','Load DTD','Position',pos + [100 0 0 0 ],'Tag','fiberGT_loaddti','Callback',{@loaddti_Callback,gcbo,[],[]},'ToolTip','Estimate Whitematter mask');
    pos = pos + [0 -25 0 0];
    uicontrol('Style','pushbutton','String','Load Mask','Position',pos,'Tag','fiberGT_loadmask','Callback',{@loadmask_Callback,gcbo,[],[]},'ToolTip','Load mrStruct with Whitematter mask');
    uicontrol('Style','text','String','<empty>','Position',pos + [205 0 0 0],'Tag','fiberGT_maskstat','Callback',{@loadmask_Callback,gcbo,[],[]},'BackGroundColor',0.8*[1 1 1],'ButtonDownFcn',@buttonmaskinfo_Callback);
    uicontrol('Style','pushbutton','String','Estimate Mask','Position',pos + [100 0 0 0],'Tag','fiberGT_estmask','Callback',{@estmask_Callback,gcbo,[],[]},'ToolTip','Estimate Whitematter mask');
   
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    uicontrol('Style','text','String','status: ready ...','Position',[5,5,300,20],'Tag','fiberGT_status','HorizontalAlignment','left');
     
    uicontrol('Style','pushbutton','String','Start Tracking','Position',[280,240,100,25],'Tag','fiberGT_startstop','Callback',{@startstop_Callback,gcbo,[],[]});
    uicontrol('Style','pushbutton','String','Reset State','Position',[680,240,100,25],'Tag','fiberGT_reset','Callback',{@reset_Callback,gcbo,[],[]});   
  
    uicontrol('Style','text','String','Tracking Log','Position',[280 425 90 20],'HorizontalAlignment','left');
    h = uicontrol('Style','listbox','String',{'Welcome to fiberGT!'},'Units','pixels','Position',[280 270 500 160],'Tag','fiberGT_log','BackGroundColor',[1 1 1],'ForeGroundColor',[1 1 1]*0);
  
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
    
    
%    sinterp = load('sinterp252struct.mat');
%    dataStruct.sphereInterpolation = sinterp;
    
    sinterp = load('sinterp128.mat');
    dataStruct.sphereInterpolation = sinterp.sinterp128;
                
    dataStruct.state = single([]);
    dataStruct.currentiteration = 1;
    dataStruct.conratio = [];   
    dataStruct.lengthhistogram = [];
    dataStruct.updownratio = [];
    dataStruct.mask_fname = '';
    dataStruct.edges = [];
    dataStruct.name = [];
    dataStruct.signal_fname = [];
    dataStruct.signal_type = [];
    dataStruct.offset = [];
    
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
        if strcmp(varargin{k},'setHARDI'),
            bTensor = varargin{k+1};
            signal = varargin{k+2};
            vox = varargin{k+3};
            patient = varargin{k+4};
            preprocessODF(bTensor,signal,patient,vox);
            addlog(sprintf('HARDI data was set manually'));
            reportstatus('ready ...');
            k = k + 4;
        elseif strcmp(varargin{k},'setDSI'),
            bTensor = varargin{k+1};
            signal = varargin{k+2};
            vox = varargin{k+3};
            patient = varargin{k+4};
            preprocessDSI(bTensor,signal,patient,vox);
            addlog(sprintf('DSI data was set manually'));
            reportstatus('ready ...');
            k = k + 4;

        elseif strcmp(varargin{k},'loadHARDI'),
            fn = varargin{k+1};
            loadhardi(fn);
            reportstatus('ready ...');
            k = k + 1;

        elseif strcmp(varargin{k},'setDTD'),
            EigStr = varargin{k+1};
            bval = varargin{k+2};
            vox = varargin{k+3};
            patient = varargin{k+4};
            sz = size(EigStr.Eval);
            [ten sig bval] = reconstructSigIntern(EigStr.Eval,EigStr.EV1,EigStr.EV2,EigStr.EV3,ones(sz(1:3)),bval,'V1.1'); 
            setdti(ten,sig,bval,patient,vox,patient,[]);            
            addlog(sprintf('DTD data was set manually'));
            reportstatus('ready ...');
            k = k + 4;
            
        elseif strcmp(varargin{k},'loadDTD'),
            fn = varargin{k+1};
            bval = varargin{k+2};
            loaddti(fn,bval);
            reportstatus('ready ...');
            k = k + 2;
            
        elseif strcmp(varargin{k},'loadMask'),
            fn = varargin{k+1};            
            threshold = varargin{k+2};
            loadmask(fn,threshold);
            reportstatus('ready ...');
            k = k + 2;
  
        elseif strcmp(varargin{k},'loadMaskStruct'),
            fn = varargin{k+1};
            name = varargin{k+2};
            mastr = maskstruct_read(fn);
            mask = maskstruct_query(mastr,'getMask',name);
            masknum = maskstruct_query(mastr,'getMaskId',name);
            if isempty(mask),
                display('mask not found!');
            else
                setWMmask(mask,sprintf('<loaded> \n\n %s %i',fn,masknum));
                addlog(sprintf('WM mask file loaded from MaskStruct'));     
                reportstatus('ready ...');
            end;    
            reportstatus('ready ...');
            k = k + 2;

        elseif strcmp(varargin{k},'setMask'),
            mask = varargin{k+1};
            if strcmp(mask,'all'),
                datastruct = get(MainHandle,'UserData');
                if isfield(datastruct,'b0avg'),
                    setWMmask(ones(size(datastruct.b0avg)),'<all>');
                    addlog(sprintf('WM mask was set manually'));     
                    reportstatus('ready ...');                    
                else
                    errordlg('Load data first!','Setting Mask');                    
                end;
            else
                setWMmask(mask>0,'<manually>');
                addlog(sprintf('WM mask was set manually'));     
                reportstatus('ready ...');                
            end;
            k = k + 1;

        elseif strcmp(varargin{k},'loadFTR'),
            fname = varargin{k+1};
            [ftr err] = ftrstruct_read(fname);
            if isempty(ftr),
                errordlg(sprintf('error: %s',err),'Loading FTR');
            else
                loadftr(ftr,fname);
            end;
            k = k + 1;
       
        elseif strcmp(varargin{k},'saveFD'),
            saveFD([],1);            

        elseif strcmp(varargin{k},'saveFDparam'),
            sel = varargin{k+1};
            szo = varargin{k+2};
            saveFD(sel,szo);            
            k = k + 2;
            
        elseif strcmp(varargin{k},'estimateMask'),
             estmask_Callback(gcbo,[],[],[],[]);
           
        elseif strcmp(varargin{k},'suggestDense'),     
            suggest_Callback(gcbo,[],[],[],[],1);
  
        elseif strcmp(varargin{k},'suggestSparse'),     
            suggest_Callback(gcbo,[],[],[],[],2);
               
          elseif strcmp(varargin{k},'setFTRname');
            setFTRname(varargin{k+1});
            k = k + 1;                
            
            
        elseif strcmp(varargin{k},'createCmat');
            roi_label_img = varargin{k+1};
            num_samples =  varargin{k+2};
            num_its =  varargin{k+3};
            Temp =  varargin{k+4};
            Nsz =  varargin{k+5};
            avgstr = varargin{k+6};
            ret = sampleCmat(roi_label_img,num_samples,num_its,Temp,Nsz,avgstr);
           
            k = k + 6;   
            
        elseif strcmp(varargin{k},'createEFTR');
            num_its =  varargin{k+1};
            num_samples =  varargin{k+2};
            Temp =  varargin{k+3};
            ret = sampleEftr(num_samples,num_its,Temp);           
            k = k + 3;                
            
        elseif strcmp(varargin{k},'genIso');
            bfac =  varargin{k+1};
            k = k + 1;           
            datastruct = get(MainHandle,'UserData');
            dirs = datastruct.sphereInterpolation.bDir;
            S = exp(-bfac*([1 0 0]*dirs).^2);
            S = S - mean(S);
            sz = size(datastruct.signal);
            for k = 1:prod(sz(2:end)),
                datastruct.signal(:,k) = S(:);
            end;            
            datastruct.spatialProbabilities = single(ones(size(datastruct.spatialProbabilities)));
            set(MainHandle,'UserData',datastruct);   
            
       elseif strcmp(varargin{k},'genCrossing');
            bfac =  varargin{k+1};
            ang = varargin{k+2};   ang2 = varargin{k+3};
            k = k + 3;           
            genCrossing(ang,ang2,bfac);
            
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

        elseif strcmp(varargin{k},'send');
            sendto_Callback(gcf,[],[]);

        elseif strcmp(varargin{k},'showMask');        
            showMask;

        elseif strcmp(varargin{k},'showSignalHist');      
            showSignalHistogram;

        elseif strcmp(varargin{k},'showParameters');      
            p = getParamsString;
            fields = fieldnames(p);
            for k = 1:length(fields);
                val = p.(fields{k});
                fprintf('%s : %s \n',fields{k},val);
            end;    
        elseif strcmp(varargin{k},'setVoxelsize');
            scaling = varargin{k+1};
            datastruct = get(MainHandle,'UserData');
            datastruct.offset(1,1) = scaling(1);
            datastruct.offset(2,2) = scaling(2);
            datastruct.offset(3,3) = scaling(3);
            set(MainHandle,'UserData',datastruct);
            infobox_update;
            k = k + 1;                
        elseif strcmp(varargin{k},'smoothODF');
            numit = varargin{k+1};
            datastruct = get(MainHandle,'UserData');
            busy;
            datastruct.signal = smoothODF(datastruct.signal,datastruct.sphereInterpolation.bDir,numit);                
            set(MainHandle,'UserData',datastruct);
            addlog(sprintf('ODF was smoothed with strength %i',numit));     
            ready;        
            k = k + 1;                
      elseif strcmp(varargin{k},'getDS');           
            ret = get(MainHandle,'UserData');
            
      elseif strcmp(varargin{k},'reactivate');           
            ready;
           
      elseif strcmp(varargin{k},'setparam');
            para = varargin{k+1};
            val = varargin{k+2};
            if isnumeric(val),
                val = num2str(val);
            end;
            h = findobj('Tag',sprintf('fiberGT_edit%s',para));
            if isempty(h),
                display('parameter not existing!');
            else
                set(h,'String',num2str(val));                        
            end;    
            k = k + 2;
        else
            fprintf('Command %s unknown!\n',['<' varargin{k} '>']);
            fprintf('Type "help fiberGT_tool" for usage information.\n');
            break;
        end;
        k = k + 1;
    end;
%catch ME1
%    display('Wrong usage! Type "help fiberGT_tool" for usage information.');
%end;
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
    set(MainHandle,'Pointer','watch');
    drawnow;
    
function ready
    h = findobj(MainHandle,'Style','pushbutton','-or','Style','edit');
    set(h,'Enable','on');
    set(MainHandle,'Pointer','arrow');
    reportstatus('ready ...');
    drawnow;
         
%%%%% overrides closefun of mainfig
function my_closereq(src,evnt)
% User-defined close request function 
% to display a question dialog box 
%    selection = questdlg('Are you sure? The tracking state will be lost!',...
%       'Close FiberGT',...
%       'Yes','No','Yes'); 
%    switch selection, 
%       case 'Yes',
%          delete(gcf)
%       case 'No'
%       return 
%    end
delete(MainHandle);

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
        sz = size(datastruct.signal);
        str = sprintf('Patient: %s\nDimensions: [%i %i %i], %i DE\nElement size [%.2f %.2f %.2f] mm',datastruct.name,sz(2:4),datastruct.original_size(4),datastruct.offset(1,1),datastruct.offset(2,2),datastruct.offset(3,3));
    end;
    set(h,'String',str);
    
    
function buttonsignalinfo_Callback(src,event)
busy;
    showSignalHistogram
ready;    

function showSignalHistogram
       ds =  get(MainHandle,'UserData');
        if isfield(ds,'signal'),
            figure;
            [hist bins] = computeSignalHistogram;
            plot(bins,hist);
        else
            display('No data has been loaded!');
        end;  
    
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
%    load/save/reset DTI/HARDI/Mask Data
%    
%    
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    
%%%%% callback for loading WM mask    
function varargout = loadmask_Callback(h,eventdata,handles,varargin)     
    loadmask('',[]);

function loadmask(fn,threshold,pmap)
        mainhandle = MainHandle; 
          ds =  get(MainHandle,'UserData');
        reportstatus('reading');
        if isempty(fn),
            [fn path] = uigetfile({'*.mat;*.nii','Accepted Files (*.mat,*.nii,*.hdr)'},'Load mrStruct/maskStruct');
            if fn == 0,
                return;
            end;
            fn = [path fn];
        end;          
        [fp fndum fext] = fileparts(fn);
        if strcmp(fext(1:4),'.nii') || strcmp(fext(1:4),'.hdr'),          
           mrdum = mrstruct_init;
           mrdum.dataAy = ds.b0avg;
           mrdum.edges = ds.edges;
           mrdum.dim3 = 'size_z';
           h = spm_vol(fn);
           if any(svd(ds.edges(1:3,1:3))*0.6 >  svd(h.mat(1:3,1:3))) %% high resolution mask
               mrdum.edges(1:3,1:3) = mrdum.edges(1:3,1:3)/2;   
               mrdum.dataAy = zeros(size(ds.b0avg)*2);               
           end;
                      [mm err] = nifti_to_mrstruct('volume',{fn},mrdum);
           if isempty(err),
                 if isempty(threshold),
                        threshold = chooseThreshold_stackview(ds.b0avg,mm.dataAy);
                        if isempty(threshold)
                           return;
                        end;
                 end; 
                 setWMmask(mm.dataAy>threshold,sprintf('<loaded> \n\n %s %f',fn,threshold));  
                 addlog(sprintf('WM mask file loaded'));    
                 reportstatus('ready ...');
                 return;
           else
                errordlg('Error during nifti_to_mrstruct','Loading Mask');
                reportstatus('ready ...');
                return;
           end;
            
    
            
        else
            if maskstruct_istype(load(fn),'Ver2')
                mastr = maskstruct_read(fn);        
                if isempty(threshold),
                    [res ok] = listdlg('ListString',mastr.maskNamesCell,'SelectionMode','single','Name','Select a Mask');
                else
                    ok = true;
                    res = threshold;
                end;
                if ok,
                    setWMmask(mastr.maskCell{res},sprintf('<loaded> \n\n %s %i',fn,res));
                    addlog(sprintf('WM mask file loaded from MaskStruct'));          
                    reportstatus('ready ...');                
                end;
            else            
                [mrstruct fname] = mrstruct_read(fn);
                if isempty(mrstruct),
                    errordlg('Error while reading','Loading Mask');
                    reportstatus('ready ...');
                else
                   if isfield(mrstruct.user,'probInfo'),
                      ok = 1;
                      if not(exist('pmap')),
                          pmap = [];
                      end;                      
                      if isempty(pmap),
                          maptype = probstruct_query(mrstruct,'mapTypes');
                          [res ok] = listdlg('ListString', maptype,'SelectionMode','single','Name','Select a Maptype');
                          pmap = maptype{res};                    
                      end;
                      if ok,
                          wm = probstruct_query(mrstruct,'getMap',pmap);   
                          if isempty(threshold),
                              threshold = chooseThreshold_stackview(ds.b0avg,wm);
                              if isempty(threshold)
                                 return;                          
                              end;                          
                          end; 
                          setWMmask(wm>threshold,sprintf('<loaded> \n\n %s %f %s',fn,threshold,pmap));
                          addlog(sprintf('WM mask file loaded from ProbMap'));        
                      end;
                      reportstatus('ready ...');
                   else
                       if not(islogical(mrstruct.dataAy)),
                            if isempty(threshold),
                                threshold = chooseThreshold_stackview(ds.b0avg,mrstruct.dataAy(:,:,:,1));
                                if isempty(threshold)
                                   return;                          
                                end;
                            end; 
                       else
                            threshold = 0.5;
                       end;
                       setWMmask(mrstruct.dataAy>threshold,sprintf('<loaded> \n\n %s %f',fn,threshold));
                       addlog(sprintf('WM mask file loaded'));          
                       reportstatus('ready ...');
                   end;
                end;
            end;
        end;
%%%%% callback for estimating WM mask    
function varargout = estmask_Callback(h,eventdata,handles,varargin)     
        reportstatus('estimating WM-mask');
        busy;
        datastruct = get(MainHandle,'UserData');
        if isfield(datastruct,'original_signal'),
            estimate = estimateWMmask(datastruct.original_signal,datastruct.original_bTensor);            
            setWMmask(estimate,'<estimate>');       
            addlog('WM mask was estimated');                       
        else
            errordlg('Load DW-data first!','Estimating Mask');       
        end;
        ready;
        reportstatus('ready ...');
   
%%%%% sets internal WM-mask to mask
function setWMmask(mask,fname)
       mainhandle = MainHandle;
       datastruct = get(mainhandle,'UserData');
       datastruct.spatialProbabilities = single(mask);         
       datastruct.mask_fname = fname;
       set(mainhandle,'UserData',datastruct);  
     
       hstat = findobj('Tag','fiberGT_maskstat');
       set(hstat,'String',fname);  

 
    
%%%%% callback for loading HARDIdata       
function varargout = loadhardi_Callback(h,eventdata,handles,varargin)  
   loadhardi('');

function loadhardi(fn)

        if isempty(fn),
            [files path] = uigetfile({'*_HARDI.mat;*.img;*.nii;*.txt', 'All accepted filetypes';'*.*','All files'},'Select HARDI data to open', 'MultiSelect', 'on');
        else
            files = fn;
            path = '';
        end;
        
        if isempty(files),
            return;
        end;
        
        reportstatus('reading HARDI mrstruct');
        busy;
        
        if iscell(files),
            files = cellfun(@(x) fullfile(path,x),files,'uniformoutput',false);
            mrstruct = nifti_to_hardi(files);
            fname = files{1};
        else
            files = fullfile(path,files);
            [p n ext] = fileparts(files);
            fname = files;
            if strcmp(ext,'.txt')
                mrstruct = nifti_to_hardi(files);
            else
                [mrstruct fname] = mrstruct_read(files);
            end;
        end;
           
        if not(isempty(mrstruct)),
            
           if false, %isfield(mrstruct.user,'mFOD'),
               preprocessFOD(mrstruct.user.bDir,mrstruct.user.sym,mrstruct.dataAy,mrstruct.patient,mrstruct.vox);              
               infobox_update ;  
               
               datastruct = get(MainHandle,'UserData');
               datastruct.signal_type = sprintf('FOD');
               datastruct.signal_fname = fname;      
               datastruct.edges = mrstruct.edges;
               sz = size(datastruct.signal);
               if isfield(datastruct,'spatialProbabilities'),
                   if any(sz(2:4) ~= size(datastruct.spatialProbabilities)),
                       datastruct = rmfield(datastruct,'spatialProbabilities');
                       hstat = findobj('Tag','fiberGT_maskstat');
                       set(hstat,'String',sprintf('<empty>'));                 
                   end;
               end;                      
               set(MainHandle,'UserData',datastruct);

               addlog(sprintf('FOD file %s loaded',fname));          
               hstat = findobj('Tag','fiberGT_hardistat');
               set(hstat,'String','FOD loaded');           
               [path fnamec] = fileparts(fname);
               setFTRname(fullfile(path,[strrep(fnamec,'FOD','') 'FTR']));
               ready;         
               
               setWMmask(mrstruct.user.mask>0.5,sprintf('<loaded> \n\n %s %f','fromFOD',0.5));
               
           else


               if ~isfield(mrstruct.user,'bTensor'),
                   ready;
                   errordlg('Gradient direction data missing!','Loading HARDI data');
                   return;
               end;


               preprocessODF(mrstruct.user.bTensor,mrstruct.dataAy,mrstruct.patient,mrstruct.vox);                    
               infobox_update ;  

               datastruct = get(MainHandle,'UserData');
               datastruct.signal_type = sprintf('HARDI');
               datastruct.signal_fname = fname;      
               datastruct.edges = mrstruct.edges;
               sz = size(datastruct.signal);
               if isfield(datastruct,'spatialProbabilities'),
                   if any(sz(2:4) ~= size(datastruct.spatialProbabilities)),
                       datastruct = rmfield(datastruct,'spatialProbabilities');
                       hstat = findobj('Tag','fiberGT_maskstat');
                       set(hstat,'String',sprintf('<empty>'));                 
                   end;
               end;                      
               set(MainHandle,'UserData',datastruct);

               addlog(sprintf('HARDI file %s loaded',fname));          
               hstat = findobj('Tag','fiberGT_hardistat');
               set(hstat,'String','HARDI loaded');           
               [path fnamec] = fileparts(fname);
               setFTRname(fullfile(path,[strrep(fnamec,'HARDI','') 'FTR']));
               ready;         
           end;
        end;


function [ten sig bval] = reconstructSignalFromTensor(mrstruct,bval,version)   
           ds = get(MainHandle,'UserData');
           reportstatus('Rebuilding signal from DTI');
           sz = size(mrstruct.b0_image_struc.dataAy);
           B0 = single(mrstruct.b0_image_struc.dataAy);
           Eval = single(mrstruct.eigenVal_struc.dataAy);
           EV1 = single(mrstruct.eigenVec_struc.dataAy(:,:,:,:,1));
           EV2 = single(mrstruct.eigenVec_struc.dataAy(:,:,:,:,2));
           EV3 = single(mrstruct.eigenVec_struc.dataAy(:,:,:,:,3));
           [ten sig bval] = reconstructSigIntern(Eval,EV1,EV2,EV3,B0,bval,version);
           
function [ten sig bval] = reconstructSigIntern(Eval,EV1,EV2,EV3,B0,bval,version)   
           ds = get(MainHandle,'UserData');
           reportstatus('Rebuilding signal from DTI');
           dirs = ds.sphereInterpolation.bDir';          
           N = size(dirs,1);
           sz = size(B0);           
           B0 = B0(:);
           Eval =  reshape(Eval,[prod(sz(1:3)) 3]);
           EV1 = reshape(EV1,[prod(sz(1:3)) 3]);
           EV2 = reshape(EV2,[prod(sz(1:3)) 3]);
           EV3 = reshape(EV3,[prod(sz(1:3)) 3]);  
           if strcmp(version,'V1.0'),
               R = [0 1 0; 1 0 0; 0 0 -1];
               EV1 = EV1*R;
               EV2 = EV2*R;
               EV3 = EV3*R;
           end;           
           for k = 1:N,
               dirs(k,:) = single(dirs(k,:) / norm(dirs(k,:)));
               ten(:,:,k) = 1000*dirs(k,:)'*dirs(k,:);
           end;
           sig = exp(-bval*((EV1*dirs').^2 .* repmat(Eval(:,1),[1 N]) + (EV2*dirs').^2 .* repmat(Eval(:,2),[1 N]) + (EV3*dirs').^2 .* repmat(Eval(:,3),[1 N])));           
           sig(:,N+1) = 1;
           sig = sig .* repmat(B0,[1 N+1]);
           
           ten(:,:,N+1) = zeros(3,3);
           sig = reshape(sig,[sz(1:3) (N+1)]);           
                
%%%%% callback for loading DTIdata                  
function varargout = loaddti_Callback(h,eventdata,handles,varargin)     
        bval = inputdlg('What b-value? (s/mm??)','Reconstruct Signal from Tensor',1,{'1000'});
        if ~isempty(bval)            
            bval = str2num(bval{1});
            loaddti([],bval);
        end;
        
function loaddti(fn,bval)
   
        reportstatus('reading DTI mrstruct');
        busy;
        if isempty(fn),
            [mrstruct err fname] = dtdstruct_read;
        else
            [mrstruct err fname] = dtdstruct_read(fn);
        end;
        if isempty(mrstruct),
            ready;      
            if ~isempty(fname),
                errordlg(sprintf('error: %s',err),'Loading DTI data');
            end;
            return;
        else
           version = 'V1.0';
           if isfield(mrstruct,'version'),
               version = mrstruct.version;
           end;                        
           [ten sig bval] = reconstructSignalFromTensor(mrstruct,bval,version);                         
           setdti(ten,sig,bval,mrstruct.b0_image_struc.patient,mrstruct.b0_image_struc.vox,fname,mrstruct.b0_image_struc.edges);
        end;      
        
function setdti(ten,sig,bval,name,vox,fname,edges)
           preprocessODF(ten,sig,name,vox);   
           infobox_update ;   

           datastruct = get(MainHandle,'UserData');           
           datastruct.signal_type = sprintf('DTI, b-value %.4f',bval);
           datastruct.signal_fname = fname;
           datastruct.edges = edges;
            
           sz = size(datastruct.signal);
           if isfield(datastruct,'spatialProbabilities'),
               if any(sz(2:4) ~= size(datastruct.spatialProbabilities)),
                   datastruct = rmfield(datastruct,'spatialProbabilities');
                   hstat = findobj('Tag','fiberGT_maskstat');
                   set(hstat,'String',sprintf('<empty>'));                 
               end;
           end;                
           set(MainHandle,'UserData',datastruct);
           
           addlog(sprintf('DTI file %s loaded',fname));
           hstat = findobj('Tag','fiberGT_hardistat');
           set(hstat,'String','DTI loaded');
           
           [path fnamec] = fileparts(fname);
           setFTRname(fullfile(path,[strrep(fnamec,'DTD','') 'FTR']));
           ready;
       
   
%%%%% callback for suggesting Parameters   
function varargout = suggest_Callback(h,eventdata,handles,varargin)  

busy;
reportstatus('suggesting parameters');
datastruct = get(MainHandle,'UserData');

hlen = findobj('Tag','fiberGT_editlength');
hwid = findobj('Tag','fiberGT_editwidth');
hwei = findobj('Tag','fiberGT_editweight');
hits = findobj('Tag','fiberGT_edititnum');

if isfield(datastruct,'offset') && isfield(datastruct,'spatialProbabilities'),
    [Hist Bins] = computeSignalHistogram;
    widthSignalDistrib = sqrt((Bins.^2)*Hist/sum(Hist));
    minelsz = min([datastruct.offset(1,1) datastruct.offset(2,2) datastruct.offset(3,3)]);
    if varargin{3} == 2,
        set(hits,'String','5*10^7');
        wmul = 1;
    else
        set(hits,'String','3*10^8'); 
        wmul = 1/4;
    end;
    set(hwid,'String',sprintf('%.3f',minelsz*0.5));
    set(hwei,'String',sprintf('%.3f',widthSignalDistrib*wmul));
    set(hlen,'String',sprintf('%.3f',minelsz*1.5));    
else
    errordlg('Load a DW-dataset and a mask first.','Parameter suggestion'); 
end;


ready;
      
function [h bins] = computeSignalHistogram    
    ds = get(MainHandle,'UserData');
    bins = -1:0.001:1;    
    
    ovsamp = max(floor(size(ds.spatialProbabilities)./size(ds.b0avg)));
    spPro = ds.spatialProbabilities;
    if ovsamp > 1,
        spPro = imresize3D_diff(spPro,size(spPro)/ovsamp,'nearest');
    end;
    
    mask = spPro(:) > 0; sigmasked = ds.signal(:,mask);
    h = histc(sigmasked(:),bins);

%%%%% set up everything st tracking can start (based on HARDI)
function preprocessODF(bTensor,signal,patient,vox)

       datastruct = get(MainHandle,'UserData');

       reportstatus('computing b0 average');       
       b0indicator = squeeze(squeeze(bTensor(1,1,:)+bTensor(2,2,:)+bTensor(3,3,:)));
       b0indicator = b0indicator/max(b0indicator);
       b0idx = find(b0indicator < 0.101);
       datastruct.b0avg = single(sum(signal(:,:,:,b0idx),4) /length(b0idx));

       reportstatus('preparing signal');
       didx = setdiff(1:size(bTensor,3),b0idx);
       datastruct.signal = single(zeros(length(didx),size(signal,1),size(signal,2),size(signal,3)));
       if isempty(b0idx),
            datastruct.signal = single(signal);
            datastruct.b0avg = single(ones(size(signal,2),size(signal,3),size(signal,4)));
       else
           for k = 1:length(didx),
                datastruct.signal(k,:,:,:) = single((signal(:,:,:,didx(k)) ./(datastruct.b0avg+0.001))); 
           end;
       end;
       datastruct.original_signal = single(signal);
       datastruct.original_bTensor = bTensor;

       
       datastruct.signal(datastruct.signal>1) = 1;
           
       reportstatus('correlating with model');          
       frt = FunkRadonTrans(bTensor(:,:,didx),datastruct.sphereInterpolation.bDir,40,0.0002);      
       frt = single(frt);
       sz = size(datastruct.signal);
       sz(1) = size(frt,1);
       datastruct.signal = single(reshape(frt*datastruct.signal(:,:),sz));
      
       datastruct.offset = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
       datastruct.offset(1,1) = vox(1);
       datastruct.offset(2,2) = vox(2);
       datastruct.offset(3,3) = vox(3);        
                   
       datastruct.name = patient;
       datastruct.original_size = size(signal);
 
       set(MainHandle,'UserData',datastruct);     

%%%%% set up everything st tracking can start (based on HARDI)
function preprocessFOD(bDir,sym,signal,patient,vox)

       datastruct = get(MainHandle,'UserData');

       datastruct.signal = permute(signal,[4 1 2 3]);
       datastruct.b0avg = squeeze(mean(datastruct.signal));
       datastruct.offset = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
       datastruct.offset(1,1) = vox(1);
       datastruct.offset(2,2) = vox(2);
       datastruct.offset(3,3) = vox(3);        
                   
       datastruct.name = patient;
       datastruct.original_size = size(signal);
 
       if sym == 1,
            sinterp = sphereInterpolLUT([bDir' ; -bDir']);
            datastruct.signal = cat(1,datastruct.signal,datastruct.signal);
       else,
            sinterp = sphereInterpolLUT([bDir']);
       end;
       
       datastruct.sphereInterpolation = sinterp;
       set(MainHandle,'UserData',datastruct);     

       
       
%%%%% set up everything st tracking can start (based on DSI)

function preprocessDSI(bTensor,signal,patient,vox)

       datastruct = get(MainHandle,'UserData');

       reportstatus('computing b0 average');      
       %% find the idx of the b=0 images
       b0indicator = squeeze(sum(sum(abs(bTensor),1),2));
       b0idx = find(b0indicator <= 0.001);
       %% compute b=0 average
       datastruct.b0avg = single(sum(signal(:,:,:,b0idx),4) /length(b0idx));

       reportstatus('preparing signal');
       %% select images b != 0
       didx = setdiff(1:size(bTensor,3),b0idx);
       datastruct.signal = single(zeros(length(didx),size(signal,1),size(signal,2),size(signal,3)));

       if isempty(b0idx), %% if there are no b=0 images just take the unnormalized DW-images
            datastruct.signal = single(signal);
            datastruct.b0avg = single(ones(size(signal,2),size(signal,3),size(signal,4)));
       else  %% normalize with b=0
           for k = 1:length(didx),
                datastruct.signal(k,:,:,:) = single((signal(:,:,:,didx(k)) ./(datastruct.b0avg+0.001))); 
           end;
       end;
       
       %% keep all the original data
       datastruct.original_signal = single(signal);
       datastruct.original_bTensor = bTensor;      
       
       %% 'unphysical' samples above 1 are set to 1
       datastruct.signal(datastruct.signal>1) = 1;
       
       %% select the directions which correspond to b!=0
       bTensor = bTensor(:,:,didx);
           
       
       reportstatus('correlating with model');          
       %% build model matrix
       Pstruc = getParams;
       bD = Pstruc.b_fac;
       bDir = datastruct.sphereInterpolation.bDir;
       for k = 1:size(bDir,2),
            for j = 1:size(bTensor,3),
                modelmat(k,j) = exp( -bD*bDir(:,k)'*bTensor(:,:,j)*bDir(:,k) /1000)*0.05;
            end;
       end;                         
       modelmat = single(modelmat);
       sz = size(datastruct.signal);
       sz(1) = size(modelmat,1);
       
       %% correlate data with model
       tmp = reshape(modelmat*datastruct.signal(:,:),sz);
       
       %% subtract mean
       tmp = tmp - repmat(mean(tmp),[size(tmp,1) 1]);
       

       %% ready and write to DS       
       datastruct.signal = single(tmp);
       
       datastruct.offset = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
       datastruct.offset(1,1) = vox(1);
       datastruct.offset(2,2) = vox(2);
       datastruct.offset(3,3) = vox(3);        
       datastruct.signal_type = 'DSI';
       datastruct.name = patient;
       datastruct.original_size = size(signal);
 
       set(MainHandle,'UserData',datastruct);     
       
       
       
       
%%%%% clears tracking state
function varargout = reset_Callback(h,eventdata,handles,varargin)  
    datastruct = get(MainHandle,'UserData');
    datastruct.state = single([]);
    datastruct.currentiteration = 1;
    datastruct.conratio = [];
    datastruct.updownratio = [];
    cla(datastruct.axeshandles(1));
    cla(datastruct.axeshandles(2));
    set(MainHandle,'UserData',datastruct);        
    addlog('Tracks deleted, state reseted');
      
    
      
    
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
info.params = getParamsString;
info.Log = getlog;
info.edges = datastruct.edges;
info.name = datastruct.name;
info.signal_fname = datastruct.signal_fname;
info.signal_type = datastruct.signal_type; 
info.mask_fname = datastruct.mask_fname;
info.conratio = datastruct.conratio;
info.lengthhistogram = datastruct.lengthhistogram;
info.updownratio = datastruct.updownratio;
info.currentiteration = datastruct.currentiteration;
    
%%%%% callback for sending fibers to fiberviewer
function varargout = sendto_Callback(h,eventdata,handles,varargin)  

datastruct = get(MainHandle,'UserData');

ftrname = getftrnameFromGUI;
    

fiblen = getFiblenFromGUI;

if ~isfield(datastruct,'state'),
    errordlg('No proper Tracking State to send!','Send fibres to Fibre Viewer');
    return;
end

if isempty(datastruct.state),
    errordlg('No fibers to send!','Send fibres to Fibre Viewer');
    return;
end;

info = gatherFTRinfo(datastruct);

reportstatus('building fibers');
fibidx = BuildFibres(datastruct.state,fiblen);
ftr = saveFTR(ftrname,datastruct.state,fibidx,datastruct.offset,info,false);

reportstatus('sending fibers to fiberviewer');
addlog(sprintf('%i fibers sent to fiberviewer',length(fibidx)));
try
    fiberviewer_model('set', 'ftrStruct',ftr);
catch
    errordlg('Problems while sending fibers to fiberviewer!','Send fibres to Fibre Viewer');
end;    

reportstatus('ready');

%%%%% callback for saving fibers
function varargout = saveftr_Callback(h,eventdata,handles,varargin)  
busy;
datastruct = get(MainHandle,'UserData');


ftrname = getftrnameFromGUI;    
fiblen = getFiblenFromGUI;

info = gatherFTRinfo(datastruct);

reportstatus('building fibers');
fibidx = BuildFibres(datastruct.state,fiblen); % reparamstep = 2

reportstatus('saving fibers');

[ftr lens] = saveFTR(ftrname,datastruct.state,fibidx,datastruct.offset,info,true);
datastruct.lengthhistogram = lenhist(lens);
updatePlots(datastruct);
set(MainHandle,'UserData',datastruct);    
addlog(sprintf('%i fibers saved',length(fibidx)));
ready;


%%%%% callback for loading fibers
function varargout = loadftr_Callback(h,eventdata,handles,varargin)  

[ftr err fname] = ftrstruct_read;
if isempty(ftr),
    %errordlg(sprintf('error: %s',err),'Loading FTR');
    ready;
    return;
end;
loadftr(ftr,fname);



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
    errordlg('Not a valid fiberGT-FTR file!','Loading FTR');
    ready;
    return;   
end;


signal_type = ftr.trackParam.signal_type;
signal_fname = ftr.trackParam.signal_fname;
if strcmp(signal_type(1:3),'HAR')
    if ~fileexists(signal_fname),
       if usejava('awt')
           uiwait(warndlg('Associated HARDI file does not exist!','Loading FTR','modal'));
       else
           error('fiberGT_tool:loadftr:noHARDI', 'Associated HARDI file ''%s'' does not exist!', signal_fname);
       end
    else
       loadhardi(signal_fname);
    end;
else
    bval = sscanf(signal_type,'DTI, b-value %f');
    loaddti(signal_fname,bval);
end;

mask_fname = ftr.trackParam.mask_fname;
if strcmp(mask_fname,'<estimate>'),
     estmask_Callback(MainHandle,[],[]);
elseif strcmp(mask_fname,'<manually>'),
     uiwait(warndlg('Manually selected mask is not inside ftr!','Loading FTR','modal'));
elseif strcmp(mask_fname,'<all>'),
        datastruct = get(MainHandle,'UserData');
        if isfield(datastruct,'b0avg'),
            setWMmask(ones(size(datastruct.b0avg)),'<all>');
            addlog(sprintf('WM mask was set manually'));     
        end
elseif strcmp(mask_fname(1:8),'<loaded>'),
   strip = mask_fname(10:end);
   [mask_fname strip] = strtok(strip);
   [threshold strip] = strtok(strip);
   [pmap strip] = strtok(strip);
   threshold = str2num(threshold);
   if ~fileexists(mask_fname),
       if usejava('awt')
           uiwait(warndlg('Associated Mask file does not exist!','Loading FTR','modal'));
           loadmask([],threshold);
       else
           error('fiberGT_tool:loadftr:nomask', 'Associated Mask file ''%s'' does not exist!', mask_fname);
       end
   else
       loadmask(mask_fname,threshold,pmap);
   end;
end;

[fpath ffname] = fileparts(fname);
setFTRname(ffname);
datastruct = get(MainHandle,'UserData');

datastruct.conratio = ftr.trackParam.conratio;
datastruct.lengthhistogram = ftr.trackParam.lengthhistogram;
datastruct.state = ftr.user;
datastruct.currentiteration = ftr.trackParam.currentiteration;
if isfield(ftr.trackParam,'updownratio'),
    datastruct.updownratio = ftr.trackParam.updownratio;
else,
    datastruct.updownratio = zeros(size(datastruct.conratio));
end;
setlog(ftr.trackParam.Log);
setParamsString(ftr.trackParam.params);
set(MainHandle,'UserData',datastruct);

updatePlots(datastruct);

ready;


function varargout = savefd_Callback(h,eventdata,handles,varargin)  
saveFD([],1);


function saveFD(sel,M)
busy;

datastruct = get(MainHandle,'UserData');
ftrname = getftrnameFromGUI;
    

fiblen = getFiblenFromGUI;

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
fibidx = BuildFibres(datastruct.state,fiblen);
ftr = saveFTR(ftrname,datastruct.state,fibidx,datastruct.offset,info,false);

reportstatus('building fiber densities');
sz = datastruct.original_size(1:3);
[fdRGB fd epdens] = ftr2FDmaps(ftr,sz,sel,M);

[fpath fname ext] = fileparts(ftrname);

reportstatus('saving fiber densities');
dummy = mrstruct_init('volume',[]);
dummy.vox = ftr.vox/M;
dummy.edges = ftr.hMatrix; dummy.edges(1:3,1:3) = dummy.edges(1:3,1:3)/M;
mrstruct_write(mrstruct_init('volume',fd,dummy),fullfile(fpath,[fname '_fd' ext]));
mrstruct_write(mrstruct_init('volumeEchos',fdRGB,dummy),fullfile(fpath,[fname '_fdrgb' ext]));
mrstruct_write(mrstruct_init('volume',epdens,dummy),fullfile(fpath,[fname '_epd' ext]));
addlog(sprintf('fiber densities were built (%i fibers)',length(fibidx)));

ready;    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%    
%    Tracking
%    
%    
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    

%%%%% callback for start tracking
function varargout = startstop_Callback(h,eventdata,handles,varargin)  

datastruct = get(MainHandle,'UserData');
if ~checkdata(datastruct),
    return;
end;

reportstatus('tracking');
Pstruc = getParams;

bfac = Pstruc.b_fac*bFacMultiplier ;


[besselexpansion meanval_sq] = computeFiberCorrelation(datastruct.sphereInterpolation.bDir,bfac);

if strcmp(datastruct.signal_type,'FOD');
    meanval_sq = 0;
end;


p_wei = Pstruc.p_weight;
params = [p_wei Pstruc.p_wid Pstruc.p_len Pstruc.c_likli Pstruc.p_chempot Pstruc.inex_balance Pstruc.p_chempot_2nd meanval_sq];
%params = [(Pstruc.p_weight*2^(3/2)) Pstruc.p_wid Pstruc.p_len Pstruc.c_likli Pstruc.p_chempot Pstruc.inex_balance Pstruc.p_chempot_2nd meanval_sq];


addlog('Tracking started');
busy;

info.params = Pstruc;
info.name = datastruct.name;

datastruct = trackit(datastruct,Pstruc.Tstart,Pstruc.Tend,Pstruc.numstep,Pstruc.numiterations,params,besselexpansion,info);
set(MainHandle,'UserData',datastruct);

ready;

addlog('Tracking stopped');
reportstatus('ready ...');

%%%%%%%%%%% get Params as Strings
function p = getParamsString

    p.Tstart = (get(findobj('Tag','fiberGT_editstartTemp'),'String'));
    p.Tend = (get(findobj('Tag','fiberGT_editstopTemp'),'String'));
    p.numstep = (get(findobj('Tag','fiberGT_editsteps'),'String'));
    p.numiterations = (get(findobj('Tag','fiberGT_edititnum'),'String'));
    p.p_weight = (get(findobj('Tag','fiberGT_editweight'),'String'));
    p.p_len = (get(findobj('Tag','fiberGT_editlength'),'String'));
    p.p_wid = (get(findobj('Tag','fiberGT_editwidth'),'String'));
    p.c_likli = (get(findobj('Tag','fiberGT_editconnlike'),'String'));
    p.b_fac = (get(findobj('Tag','fiberGT_editbfac'),'String'));
    p.p_chempot = (get(findobj('Tag','fiberGT_editchemPot1'),'String'));
    p.p_chempot_2nd = (get(findobj('Tag','fiberGT_editchemPot2'),'String'));
    p.inex_balance = (get(findobj('Tag','fiberGT_editinexbal'),'String'));
    p.fibrange = getFiblenFromGUI;

%%%%% sets tracking params in GUI
function setParamsString(p)

    set(findobj('Tag','fiberGT_editstartTemp'),'String',p.Tstart);
    set(findobj('Tag','fiberGT_editstopTemp'),'String',p.Tend);
    set(findobj('Tag','fiberGT_editsteps'),'String',p.numstep);
    set(findobj('Tag','fiberGT_edititnum'),'String',p.numiterations);
    set(findobj('Tag','fiberGT_editweight'),'String',p.p_weight);
    set(findobj('Tag','fiberGT_editlength'),'String',p.p_len);
    set(findobj('Tag','fiberGT_editwidth'),'String',p.p_wid);
    set(findobj('Tag','fiberGT_editconnlike'),'String',p.c_likli);
    set(findobj('Tag','fiberGT_editbfac'),'String',p.b_fac);
    set(findobj('Tag','fiberGT_editchemPot1'),'String',p.p_chempot);
    set(findobj('Tag','fiberGT_editchemPot2'),'String',p.p_chempot_2nd);
    set(findobj('Tag','fiberGT_editinexbal'),'String',p.inex_balance);
    if isfield(p,'fibrange'),
        set(findobj('Tag','fiberGT_editfiblength'),'String',sprintf('%i;%i',p.fibrange));
    end;
  
%%%%% gets tracking params from GUI as numbers
function p = getParams
    p = getParamsString;
    fields = fieldnames(p);
    for k = 1:length(fields);
        val = p.(fields{k});
        if ischar(val),
            val = str2num(val);
        end;
        p.(fields{k}) = val;
    end;
  
function fiblen = getFiblenFromGUI
    hfl = findobj('Tag','fiberGT_editfiblength');
    fiblen = eval(['[' get(hfl,'String') ']']);
    if length(fiblen) == 1,
        fiblen(2) = inf;
    end;
        
    


    
function m = bFacMultiplier    
m =1;
    
%%%%% checks data before tracking is started
function ret = checkdata(datastruct)

ret = false;
if ~isfield(datastruct,'spatialProbabilities'),
    errordlg('No White Matter Mask loaded','Data consistency');
    return
end;

if ~isfield(datastruct,'signal'),
    errordlg('No HARDI signal loaded','Data consistency');
    return
end;

if ~isfield(datastruct,'sphereInterpolation');
    errordlg('Sphere Interpolator not initialized (missing sinterp128.mat)','Data consistency');
    return
end;

szhardi = size(datastruct.signal);
szmask = size(datastruct.spatialProbabilities);
if ~isequal(mod(szmask,szhardi(2:4)),[0 0 0]),
    errordlg('dimensions of WM mask and HARDI signal not consistent','Data consistency');
    return
end;

if any(cellfun(@isempty,struct2cell(getParams))),
    errordlg('Parameters are not set correct','Data consistency');
    return
end;    

ret = true;
return;

%%%%% tracking
function datastruct = trackit(datastruct,Tstart,Tend,ni,its,params,besselexpansion,info)

if its > 1,
    ipr = round(its/ni);
else
    ipr = its;
end;

lens = [];
h = createStopButton;
sz = size(datastruct.spatialProbabilities);
alpha = log(Tend/Tstart);
vox = [datastruct.offset(1,1) datastruct.offset(2,2) datastruct.offset(3,3)];


ftrname = getftrnameFromGUI;
fiblen = getFiblenFromGUI;
    
doexit = false;
cnt = 1;
for k = datastruct.currentiteration:ni,    
    
      datastruct.currentiteration = k;
    
      tic

      T = Tstart * exp(alpha*(k-1)/ni);
      
      [datastruct.state stat] = pcRJMCMC(datastruct.state,datastruct.signal,datastruct.spatialProbabilities,double(vox),double([T ipr 0.35 params]),single(besselexpansion),datastruct.sphereInterpolation,h);    
      
      concnt = (sum(datastruct.state(9,:) ~= -1) + sum(datastruct.state(10,:) ~= -1))/2; 
      segcnt = size(datastruct.state,2);
     
      fibidx = BuildFibres(datastruct.state,fiblen);
     
      addlog(sprintf('%i/%i) #p %i, #c %i, #f %i, T %.4f, t per %.0e its :%.2f min\n',k,ni,segcnt,concnt,length(fibidx), T,stat.iterations, toc/60));
       
      datastruct.conratio(k) = concnt/segcnt;
      datastruct.updownratio(k) = stat.up/(stat.down+1);
 
      datastruct.lengthhistogram = lenhist(lens);

      updatePlots(datastruct);
    
      drawnow;
   

      if isfield(datastruct,'maxnumits'),
          if cnt >= datastruct.maxnumits,
              datastruct.currentiteration = k+1;
              doexit = true;
          end;
      end;      
      cnt = cnt + 1;
            
      
     if stopButtonPressed,
          addlog('interrupted by user');
          doexit = true;
      end;
      
      info = gatherFTRinfo(datastruct);
      [ftr lens] = saveFTR(ftrname,datastruct.state,fibidx,datastruct.offset,info);
      

      if doexit,
          break;
      end;
      
end;


if ~stopButtonPressed,
    close(h)
end;


%%%%%%%%%% cmat sampling




function ret = sampleEftr(num_samples,num_its,Temp);     

sampleparams.num_samples = num_samples;
sampleparams.num_its = num_its;
sampleparams.Temp = Temp;

ret = startSampler(sampleparams);


function ret = sampleCmat(rois,num_samps,num_its,Temp,Nsz,avgstr)

sampleparams.rois = rois;
sampleparams.Nsz = Nsz;
sampleparams.avgstr = avgstr;
sampleparams.num_samples = num_samps;
sampleparams.num_its = num_its;
sampleparams.Temp = Temp;

ret = startSampler(sampleparams);


function ret = startSampler(sampleparams)

datastruct = get(MainHandle,'UserData');
if ~checkdata(datastruct),
    return;
end;

reportstatus('sampling Cmat');
Pstruc = getParams;

bfac = Pstruc.b_fac*bFacMultiplier ;

[besselexpansion meanval_sq] = computeFiberCorrelation(datastruct.sphereInterpolation.bDir,bfac);


if strcmp(datastruct.signal_type,'FOD');
    meanval_sq = 0;
end;

if isempty(sampleparams.Temp),
    alpha = log(Pstruc.Tend/Pstruc.Tstart);
    sampleparams.Temp = Pstruc.Tstart * exp(alpha*(datastruct.currentiteration-1)/Pstruc.numstep);
end;


p_wei = Pstruc.p_weight;
params = [p_wei Pstruc.p_wid Pstruc.p_len Pstruc.c_likli Pstruc.p_chempot Pstruc.inex_balance Pstruc.p_chempot_2nd meanval_sq];
%params = [(Pstruc.p_weight*2^(3/2)) Pstruc.p_wid Pstruc.p_len Pstruc.c_likli Pstruc.p_chempot Pstruc.inex_balance Pstruc.p_chempot_2nd meanval_sq];


addlog('Sampling started');
busy;

info.params = Pstruc;
info.name = datastruct.name;

ret = sampleit(datastruct,params,besselexpansion,info,sampleparams);

ready;

addlog('Sampling stopped');
reportstatus('ready ...');




%%%%% tracking
function ret = sampleit(datastruct,params,besselexpansion,info,sampleparams)

lens = [];
h = createStopButton;
sz = size(datastruct.spatialProbabilities);
vox = [datastruct.offset(1,1) datastruct.offset(2,2) datastruct.offset(3,3)];


ftrname = getftrnameFromGUI;
fiblen = getFiblenFromGUI;
    
mstruc.cc = 0;
mstruc.lens = 0;

cnt = 1;

ftr_all = [];

for k = 1:sampleparams.num_samples,
        
      tic
      
      [datastruct.state stat] = pcRJMCMC(datastruct.state,datastruct.signal,datastruct.spatialProbabilities,double(vox),...
          double([sampleparams.Temp sampleparams.num_its 0.35 params]),single(besselexpansion),datastruct.sphereInterpolation,h);    
  
      concnt = (sum(datastruct.state(9,:) ~= -1) + sum(datastruct.state(10,:) ~= -1))/2; 
      segcnt = size(datastruct.state,2);
     
      fibidx = BuildFibres(datastruct.state,fiblen);
     
      addlog(sprintf('%i/%i) #p %i, #c %i, #f %i, T %.4f, t per %.0e its :%.2f min\n',k,sampleparams.num_samples,segcnt,concnt,length(fibidx), sampleparams.Temp,round(sampleparams.num_its), toc/60));

      info = gatherFTRinfo(datastruct);
      
      [ftr lens] = saveFTR(ftrname,datastruct.state,fibidx,datastruct.offset,info,false);
 
      
      if isfield(sampleparams,'rois'),

          newmstruc = createCM_GT(ftr,sampleparams.rois,sampleparams.Nsz);

          if strcmp(sampleparams.avgstr,'noavg'),
              mstruc_noavg(k) = newmstruc;          
          else
              mstruc.cc = mstruc.cc + newmstruc.cc;
              mstruc.lens = mstruc.lens + newmstruc.lens;
          end;
      end;
      ftr_all = mergeFTR(ftr,ftr_all);          
      
      cnt = cnt + 1;
      
      datastruct.conratio(k) = concnt/segcnt;
      datastruct.updownratio(k) = stat.up/(stat.down+1);
 
      datastruct.lengthhistogram = lenhist(lens);

      updatePlots(datastruct);
    
      drawnow;
   
     if stopButtonPressed,
         addlog('interrupted by user');
          break;
      end;
end;

if isfield(sampleparams,'rois'),
    if strcmp(sampleparams.avgstr,'noavg'),
        mstruc.cc = cat(3,mstruc_noavg(:).cc);
        mstruc.lens = cat(3,mstruc_noavg(:).lens);    
    else,    
        mstruc.cc = mstruc.cc/cnt;
        mstruc.lens = mstruc.lens/cnt;
    end;
    ret = mstruc;
    ret.ftr = ftr_all;
else
    ret = ftr_all;
end;


if ~stopButtonPressed,
    close(h)
end;
 




function ftr1 = mergeFTR(ftr1,ftr2)

if isempty(ftr2),
    return;
end;

ftr1.curveSegCell = {ftr1.curveSegCell{:} ftr2.curveSegCell{:}};
ftr1.fiber = {};
ftr1.user = [];
ftr1.connectCell = arrayfun(@(x) x, 1:length(ftr1.curveSegCell),'uniformoutput',false)';
ftr1.curveSegCell = ftr1.curveSegCell';
 







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
        plot(datastruct.axeshandles(1),1:n,datastruct.conratio,'b',1:n,datastruct.updownratio,'r'); 
        
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
  

    
    
function genCrossing(angle,ang2,bfac)    
    
            datastruct = get(MainHandle,'UserData');
            dirs = datastruct.sphereInterpolation.bDir;
            N1 = [cos(ang2) sin(ang2) 0];
            S1 = exp(-bfac*(N1*dirs).^2);
           
            N2 = [sin(angle) cos(angle) 0];
            S2 = exp(-bfac*(N2*dirs).^2);
           
            
            d = 2^2;
            sig = zeros(length(dirs),64,64);
            [X Y] = ndgrid(-31:32);
            R2 = X.^2 + Y.^2; 
            vx =  (R2-(N1(1)*X+N1(2)*Y).^2 < d); vx = find(vx(:));
            for k = 1:length(vx),
                sig(:,vx(k)) = sig(:,vx(k)) + S1(:);
            end;
            vx =  (R2-(N2(1)*X+N2(2)*Y).^2 < d); vx = find(vx(:));
            for k = 1:length(vx),
                sig(:,vx(k)) = sig(:,vx(k)) + S2(:);
            end;
            mask = abs(squeeze(sig(1,:,:)))>0;
            
            
            sig = permute(sig,[2 3 1]);
            for k = 1:3,
                signal(:,:,k,:) = reshape(sig,[64 64 1 length(dirs)]);
                datastruct.spatialProbabilities(:,:,k) = mask;
            end;            
            
            datastruct.spatialProbabilities = single(datastruct.spatialProbabilities);
            set(MainHandle,'UserData',datastruct);    
            for k = 1:length(dirs),
                tens(:,:,k) = dirs(:,k)*(dirs(:,k))';
            end;
            
            tens(:,:,end+1) = 0;
            signal(:,:,:,end+1) = 1;
            preprocessODF(tens,signal,'crossing',[3 3 3]);
    
    
    
    
