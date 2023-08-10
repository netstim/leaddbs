function ea_export_server(hobj,ev,options)

disp('Exporting data to LEAD Server...');
%keyboard
if ~isfield(options.prefs.ls,'dir')
    % configure server output directory for the first time.
    serverdir=uigetdir([],'Please choose output directory for LEAD-server.');
    if ~serverdir % user pressed cancel.
        return
    end
    serverdir=[serverdir,filesep]; % append /

    % store directory in ea_prefs
    fid = fopen([ea_getearoot,'ea_prefs.m'],'a');
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

if ~exist(options.prefs.ls.dir,'file')
    mkdir(options.prefs.ls.dir);
end

if ~exist([options.prefs.ls.dir,'data'],'file')
    mkdir([options.prefs.ls.dir,'data']);
end

if ~exist([options.prefs.ls.dir,'data',filesep,options.patientname],'file')
    mkdir([options.prefs.ls.dir,'data',filesep,options.patientname]);
end

% export model
bbstruct=getappdata(gcf,'bbstruct');
if ~isempty(bbstruct)
    savejson('',bbstruct,'FileName',[options.prefs.ls.dir,'data',filesep,options.patientname,filesep,'bb_scene.json'],'ArrayToStruct',0);
    % export html
    copyfile([options.earoot,'ls',filesep,'index.html'],[options.prefs.ls.dir,'data',filesep,options.patientname,filesep,'index.html']);
else
    warning('JSON-file not set ? set electrode rendering to standard electrodes display and re-run LEAD for export in webserver.');
    return
end

% export html
copyfile([options.earoot,'ls',filesep,'index.html'],[options.prefs.ls.dir,'data',filesep,options.patientname,filesep,'index.html']);
fid = fopen([options.prefs.ls.dir,'data',filesep,options.patientname,filesep,'index.html'], 'a');
fwrite(fid,'</body> ');
fwrite(fid,'</html>');

% export ftracking and vat results if present
PL=getappdata(gcf,'PL');
stimparams=getappdata(gcf,'stimparams');
if ~isempty(PL)
    stimparamsstr{1}='Stimulation parameters:';
    stimparamsstr{2}='Right Hemisphere:';
    stimparamsstr{3}=['Voltage: ',num2str(stimparams(1).U),'.'];
    stimparamsstr{4}=['Impedance: ',num2str(stimparams(1).Im),'.'];
    stimparamsstr{5}='Left Hemisphere:';
    stimparamsstr{6}=['Voltage: ',num2str(stimparams(2).U),'.'];
    stimparamsstr{7}=['Impedance: ',num2str(stimparams(2).Im),'.'];

    nowstr=[date,'_'];
    c=clock;
    nowstr=[nowstr,num2str(c(4)),'_',num2str(c(5))];
    clear c
    mkdir([options.prefs.ls.dir,'data',filesep,options.patientname,filesep,nowstr])
    save([options.prefs.ls.dir,'data',filesep,options.patientname,filesep,nowstr,filesep,'stimparams.mat'],'stimparams','stimparamsstr');

    vatstruct=ea_viz2brainbrowser(PL.vatfv);
    savejson('',vatstruct,'FileName',[options.prefs.ls.dir,'data',filesep,options.patientname,filesep,nowstr,filesep,'bb_vat.json'],'ArrayToStruct',0);

    % export html
    copyfile([options.earoot,'ls',filesep,'index_vat.html'],[options.prefs.ls.dir,'data',filesep,options.patientname,filesep,nowstr,filesep,'index.html']);
    fid = fopen([options.prefs.ls.dir,'data',filesep,options.patientname,filesep,nowstr,filesep,'index.html'], 'a');
    fwrite(fid,'<div><p>');
    fwrite(fid,'<span style="font-family:Arial,Arial">');
    fwrite(fid,'<span style="font-weight:bold">');
    fwrite(fid,['<br>',stimparamsstr{1}]);
    fwrite(fid,'</span>');
    for s=2:length(stimparamsstr)
        fwrite(fid,['<br>',stimparamsstr{s}]);
    end

    fwrite(fid,'</span></p></div>');
    fwrite(fid,'</body> ');
    fwrite(fid,'</html>');

    try
        fibstruct=ea_viz2brainbrowser(PL.bbfibfv,'line');
        %fibstruct=rmfield(fibstruct,'normals');
        savejson('',fibstruct,'FileName',[options.prefs.ls.dir,'data',filesep,options.patientname,filesep,nowstr,filesep,'bb_fibs.json'],'ArrayToStruct',0);

        % export html
        copyfile([options.earoot,'ls',filesep,'index_fibs.html'],[options.prefs.ls.dir,'data',filesep,options.patientname,filesep,nowstr,filesep,'index.html']);
        fid = fopen([options.prefs.ls.dir,'data',filesep,options.patientname,filesep,nowstr,filesep,'index.html'], 'a');
        fwrite(fid,'<div><p>');
        for s=1:length(stimparamsstr)
            fwrite(fid,['<br>',stimparamsstr{s}]);
        end
        fwrite(fid,'</div></p>');
        fwrite(fid,'</body> ');
        fwrite(fid,'</html>');
    end
end

mkdir([options.prefs.ls.dir,'data',filesep,options.patientname,filesep,'slices']);
orients={'axial','coronal','saggital'};
oricon=zeros(3,16); % matrix that shows which images are available
for con=0:15
    for or=1:length(orients)
        if exist([options.root,options.patientname,filesep,'K',num2str(con),'_',orients{or},'.png'],'file')
            copyfile([options.root,options.patientname,filesep,'K',num2str(con),'_',orients{or},'.png'],[options.prefs.ls.dir,'data',filesep,options.patientname,filesep,'slices',filesep,'K',num2str(con),'_',orients{or},'.png']);
            oricon(or,con+1)=1;
        end
    end
end

% export slices html:
if any(oricon(:))
    fid = fopen([options.prefs.ls.dir,'data',filesep,options.patientname,filesep,'slices',filesep,'index.html'], 'w+');

    fwrite(fid,'<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">');
    fwrite(fid,'<html xmlns="http://www.w3.org/1999/xhtml">');
    fwrite(fid,'<head>');
    fwrite(fid,'<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />');
    fwrite(fid,['<title>',options.patientname,': Slice views','</title>']);
    fwrite(fid,'<style type="text/css">');
    fwrite(fid,'<!--');
    fwrite(fid,' body {');
    fwrite(fid,'  font-family: Arial;');
    fwrite(fid,'  text-align:left;');
    fwrite(fid,' }');
    fwrite(fid,'//-->');
    fwrite(fid,'</style>');
    fwrite(fid,'</head>');
    fwrite(fid,'');
    fwrite(fid,'<body>');
    fwrite(fid,['<p style="text-align:left"><h2>Slice cut views for ',options.patientname,'</h2>']);
    fwrite(fid,['<p style="text-align:left">Contacts K0-K7 refer to the right hemisphere, K8-K15 to the left hemisphere.']);
    fwrite(fid,['</p>']);

    allorients=[1:3];
    for or=allorients(logical(sum(oricon,2)))
        fwrite(fid,['<p style="text-align:left"><h2>',orients{or},' views</h2></p>']);
        allcontacts=[0:15];
        for con=allcontacts(logical(oricon(or,:))')
            fwrite(fid,['<p style="text-align:left"><h3>K ',num2str(con),':</h3></p>']);
            fwrite(fid,['<p style="text-align:left"><a href="K',num2str(con),'_',orients{or},'.png"><img src="K',num2str(con),'_',orients{or},'.png" alt="L',num2str(con),'_',orients{or},'" width="800"></a></p>']);
        end
    end

    fwrite(fid,'</body>');
    fwrite(fid,'</html>');
    fclose(fid);
else
    rmdir([options.prefs.ls.dir,'data',filesep,options.patientname,filesep,'slices'],'s');
end

% append to index
ea_export_ls_index(options.prefs.ls.dir,options);


function initializeserver(serverdir,options)

% initialize server inside that folder...
copyfile([options.earoot,'ls',filesep,'ls_bb_extract.zip'],[serverdir,'ls_bb_extract.zip']);
% extract zipfile
unzip([serverdir,'ls_bb_extract.zip'],serverdir);
delete([serverdir,'ls_bb_extract.zip']);
