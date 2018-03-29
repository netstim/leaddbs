function success=ea_hastoolbox(cmd)
success=0;
switch cmd
    case 'freesurfer'
        
        if ispc
           msgbox('Freesurfer is not available for Windows.'); 
           return
        end
        
        checkfsstandard
        
        options.prefs=ea_prefs;
        while ~isfield(options.prefs,'fspath')
            msgbox('FreeSurfer installation not set properly, please select FreeSurfer base folder in the next step.');
            pth=uigetdir('','Please select freesurfer installation folder');
            if ~ischar(pth) % user pressed cancel
                return
            end
            if exist(fullfile(pth,'bin','recon-all'),'file')
                success=1;
                ea_injectprefstring(['prefs.fspath=''',pth,filesep,''';']);
            end
            options.prefs=ea_prefs;
        end
        setenv('FREESURFER_HOME',options.prefs.fspath);
        system(['source ',options.prefs.fspath,filesep,'SetUpFreeSurfer.sh']);
        setenv('PATH', [getenv('PATH') ':',options.prefs.fspath,'bin']);
        setenv('PATH', [getenv('PATH') ':',options.prefs.fspath,'mni/bin']);
        if ismac
            copyfile([ea_gethome,'.bash_profile'],[ea_getearoot,'bp']);
            system(['chmod +x ',ea_getearoot,'bp']);
            system([ea_getearoot,'bp']);
            delete([ea_getearoot,'bp']);
        end
        success=1;
        
    case 'fsl'
        if ispc
            msgbox('FSL is not available for Windows.');
            return
        end
        
        
        options.prefs=ea_prefs;
        while ~isfield(options.prefs,'fsldir')
            msgbox('FSL installation not set properly, please select FSL base folder in the next step.');
            pth=uigetdir('','Please select FSL installation folder');
            if ~ischar(pth) % user pressed cancel
                return
            end
            if exist(fullfile(pth,'bin',filesep,'fsl'),'file')
                success=1;
                ea_injectprefstring(['prefs.fsldir=''',pth,filesep,''';']);
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
end






function checkfsstandard % check if standard installations exist.

if ismac
    if exist(fullfile('/Applications','freesurfer','bin','recon-all'),'file')
       ea_injectprefstring(['prefs.fspath=''','/Applications/freesurfer',filesep,''';']);
    end
end
        
        
        
        
        
        
        
        