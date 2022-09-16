function success=ea_hastoolbox(cmd)
success=0;
switch cmd
    case 'freesurfer'
        
        if ispc
           errordlg('FreeSurfer is not available for Windows.');
           return
        end
        
        checkfsstandard
        
        options.prefs=ea_prefs;
        while ~isfield(options.prefs,'fspath')
            errordlg('FreeSurfer installation not set properly, please select FreeSurfer base folder in the next step.');
            uiwait
            pth=uigetdir('','Please select FreeSurfer installation folder');
            if ~ischar(pth) % user pressed cancel
                return
            end
            if exist(fullfile(pth,'bin','recon-all'),'file')
                success=1;
                ea_injectprefstring('fspath',[pth,filesep]);
            end
            options.prefs=ea_prefs;
        end
        setenv('FREESURFER_HOME',options.prefs.fs.dir);
        system(['source ',options.prefs.fs.dir,filesep,'SetUpFreeSurfer.sh']);
        setenv('PATH', [getenv('PATH') ':',options.prefs.fs.dir,'bin']);
        setenv('PATH', [getenv('PATH') ':',options.prefs.fs.dir,'mni/bin']);
        if ismac
            try
                copyfile([ea_gethome,'.bash_profile'],[ea_getearoot,'bp']);
                system(['chmod +x ',ea_getearoot,'bp']);
                system([ea_getearoot,'bp']);
                delete([ea_getearoot,'bp']);
            end
        end
        success=1;
        
    case 'fsl'
        if ispc
            errordlg('FSL is not available for Windows.');
            return
        end
        
        
        options.prefs=ea_prefs;
        while ~isfield(options.prefs,'fsldir')
            errordlg('FSL installation not set properly, please select FSL base folder in the next step.');
            uiwait
            pth=uigetdir('','Please select FSL installation folder');
            if ~ischar(pth) % user pressed cancel
                return
            end
            if exist(fullfile(pth,'bin',filesep,'fsl'),'file')
                success=1;
                ea_injectprefstring('fsldir',[pth,filesep]);
            end
            options.prefs=ea_prefs;
        end
        success=1;
        
        system(['FSLDIR=',options.prefs.fsldir]);
        system(['PATH=${FSLDIR}/bin:${PATH}']);
        system([options.prefs.fsldir,'etc/fslconf/fsl.sh']);
        system(['export FSLDIR PATH']);
        setenv('PATH', [getenv('PATH') ':',options.prefs.fsldir,'bin']);
    case 'cat'
        success=exist([spm('dir'),filesep,'toolbox',filesep,'cat12'],'dir');
    case 'slicer'
        checkslicerstandard;
        options.prefs=ea_prefs;
        while ~isfield(options.prefs,'slicer') || isempty(options.prefs.slicer.dir)
            errordlg('3DSlicer installation not set properly, please select 3DSlicer in the next step.');
            uiwait
            [fn,pth]=uigetfile('*.*','Please select 3DSlicer application');
            if ~ischar(fn) % user pressed cancel
                return
            end
                success=1;
                ea_injectprefstring('slicer','dir',fullfile(pth,fn));
            options.prefs=ea_prefs;
        end
        success=1;
end






function checkfsstandard % check if standard installations exist.

if ismac
    if exist(fullfile('/Applications','freesurfer','bin','recon-all'),'file')
       ea_injectprefstring('fspath',['/Applications/freesurfer',filesep]);
    end
end
        
        
      

function checkslicerstandard % check if standard installations exist.

if ismac
    if exist(fullfile('/Applications','Slicer.app'),'file')
       ea_injectprefstring('slicer','dir',fullfile('/Applications','Slicer.app'));
    end
end
        
        
        
        
