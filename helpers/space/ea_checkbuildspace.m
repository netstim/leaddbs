function ea_checkbuildspace


if exist([ea_getearoot,'templates',filesep,'space',filesep,ea_getspace,filesep,'need_build'],'file')
    answ=questdlg('The template space you selected needs to be unpacked/built before it may be used. This process can take several hours and will take up additional space on your harddrive. Do you wish to proceed unpacking the data?',...
        'Unpack template space','Sure','Cancel','Cancel');
    
    switch answ
        case 'Sure'
            ea_unpackspace
        case 'Cancel'
            msgbox('It is recommended to switch the space back to a different space.');
    end
elseif exist([ea_getearoot,'templates',filesep,'space',filesep,ea_getspace,filesep,'need_install'],'file')
    answ=questdlg('The template space you selected needs to be installed and unpacked/built before it may be used. This process needs an internet connection and can take several hours and will take up additional space on your harddrive. Do you wish to proceed unpacking the data?',...
        'Unpack template space','Sure','Cancel','Cancel');
    
    switch answ
        case 'Sure'
            ea_installspace
        case 'Cancel'
            msgbox('It is recommended to switch the space back to a different space.');
    end
    
end