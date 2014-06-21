function ea_export_server(hobj,ev,options)

disp('Exporting data to LEAD Server...');
if ~isfield(options.prefs.ls,'dir')
    % configure server output directory for the first time.
    
    serverdir=uigetdir([],'Please choose output directory for LEAD-server.');
    
    if ~serverdir % user pressed cancel.
        return
    end
    serverdir=[serverdir,filesep]; % append /
    % store directory in ea_prefs
    fid = fopen([fileparts(which('lead')),filesep,'ea_prefs.m'],'a');
    fwrite(fid,['prefs.ls.dir=','''',serverdir,'''',';']);
    fclose(fid);
    options.prefs.ls.dir=serverdir;
    
    if ~exist([serverdir,'server.js'],'file')
        initializeserver(serverdir,options);
    end
    
    % try to start server using system events
    try
        system(['cd ',serverdir]);
        system(['node server.js']);
        disp('LEAD-Server seems to have been established successfully.');
        disp('Navigate to http://localhost:5000/ in your browser to see results.');
    end
    
end

if ~exist(options.prefs.ls.dir,'file');
    mkdir(options.prefs.ls.dir);
end
if ~exist([options.prefs.ls.dir,'data'],'file');
    mkdir([options.prefs.ls.dir,'data']);
end
if ~exist([options.prefs.ls.dir,'data',filesep,options.patientname],'file');
    
    mkdir([options.prefs.ls.dir,'data',filesep,options.patientname]);
end
% export model
bbstruct=getappdata(gcf,'bbstruct');
if ~isempty(bbstruct)
    ea_savejson('',bbstruct,'FileName',[options.prefs.ls.dir,'data',filesep,options.patientname,filesep,'bb_scene.json'],'ArrayToStruct',0);
    % export html
    copyfile([options.earoot,'ls',filesep,'index.html'],[options.prefs.ls.dir,'data',filesep,options.patientname,filesep,'index.html']);
    
else
    warning('JSON-file not set ? set electrode rendering to standard electrodes display and re-run LEAD for export in webserver.');
    return
end
% export html
copyfile([options.earoot,'ls',filesep,'index.html'],[options.prefs.ls.dir,'data',filesep,options.patientname,filesep,'index.html']);


% export ftracking and vat results if present
PL=getappdata(gcf,'PL');
if ~isempty(PL)
    nowstr=[date,'_',num2str(now)];
    nowstr=nowstr(1:end-5);
    mkdir([options.prefs.ls.dir,'data',filesep,options.patientname,filesep,nowstr])
    vatstruct=ea_viz2brainbrowser(PL.vatfv);
    
    ea_savejson('',vatstruct,'FileName',[options.prefs.ls.dir,'data',filesep,options.patientname,filesep,nowstr,filesep,'bb_vat.json'],'ArrayToStruct',0);
    % export html
    copyfile([options.earoot,'ls',filesep,'index_vat.html'],[options.prefs.ls.dir,'data',filesep,options.patientname,filesep,nowstr,filesep,'index.html']);
    
    try
        
        fibstruct=ea_viz2brainbrowser(PL.bbfibfv,'line');
        %fibstruct=rmfield(fibstruct,'normals');
        
        ea_savejson('',fibstruct,'FileName',[options.prefs.ls.dir,'data',filesep,options.patientname,filesep,nowstr,filesep,'bb_fibs.json'],'ArrayToStruct',0);
        % export html
        copyfile([options.earoot,'ls',filesep,'index_fibs.html'],[options.prefs.ls.dir,'data',filesep,options.patientname,filesep,nowstr,filesep,'index.html']);
        
    end

    
    
    
    
end

% append to index
ea_export_ls_index(options.prefs.ls.dir,options);


function initializeserver(serverdir,options)

% initialize server inside that folder...

copyfile([options.earoot,'ls',filesep,'ls_bb_extract.zip'],[serverdir,'ls_bb_extract.zip']);
% extract zipfile
unzip([serverdir,'ls_bb_extract.zip'],serverdir);
delete([serverdir,'ls_bb_extract.zip']);