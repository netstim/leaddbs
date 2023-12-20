function ea_checkbuildspace

if isfile([ea_space,'need_build'])
    answ=questdlg('The template space you selected needs to be unpacked/built before it may be used. This process can take several hours and will take up additional space on your harddrive. Do you wish to proceed unpacking the data?',...
        'Unpack template space','Sure','Cancel','Cancel');

    switch answ
        case 'Sure'
            ea_unpackspace
        case 'Cancel'
            msgbox('It is recommended to switch the space back to a different space.');
    end
elseif isfile([ea_space,'need_install'])
    % in user environment ('need_install' will be deleted after space installed) OR
    % in dev environment('need_install' is always there) && space not installed
    if ~(isfolder([ea_getearoot,'.git']) && isfile([ea_space,'wires.mat']))
        answ=questdlg('The template space you selected needs to be installed and unpacked/built before it may be used. This process needs an internet connection and can take several hours and will take up additional space on your harddrive. Do you wish to proceed unpacking the data?',...
            'Unpack template space','Sure','Cancel','Cancel');

        switch answ
            case 'Sure'
                ea_installspace
            case 'Cancel'
                msgbox('It is recommended to switch the space back to a different space.');
        end
    end
end
