function ea_flattenfiducialhelpers(~,~,handles)

uipatdir=getappdata(handles.leadfigure,'uipatdir');
cnt=1;
for pt=1:length(uipatdir)
    fids=dir([uipatdir{pt},filesep,'fiducials',filesep,'native',filesep,'*.nii.gz']);
    clear fidcell tfidcell
    for fi=1:length(fids)
        % pt folder
        gunzip([uipatdir{pt},filesep,'fiducials',filesep,'native',filesep,fids(fi).name]);
        delete([uipatdir{pt},filesep,'fiducials',filesep,'native',filesep,fids(fi).name]);
        fidcell{fi}=[uipatdir{pt},filesep,'fiducials',filesep,'native',filesep,fids(fi).name(1:end-3)];
        
        % templates
        if exist([uipatdir{pt},filesep,'fiducials',filesep,ea_getspace,filesep,fids(fi).name],'file')
            gunzip([uipatdir{pt},filesep,'fiducials',filesep,ea_getspace,filesep,fids(fi).name]);
            delete([uipatdir{pt},filesep,'fiducials',filesep,ea_getspace,filesep,fids(fi).name]);
            tfidcell{fi}=[uipatdir{pt},filesep,'fiducials',filesep,ea_getspace,filesep,fids(fi).name(1:end-3)];
        end
    end
    
    % flatten ROI:
    fguid=ea_generate_uuid;
    for tc=1:length(tfidcell)
        nii=ea_load_nii(tfidcell{tc});
        nii.img=nii.img./max(nii.img(:));
        if ~exist('AllX','var')
            AllX=nii.img;
        else
            try
                AllX=[AllX+nii.img];
            catch
                ea_conformspaceto(tfidcell{1},tfidcell{tc},1);
                nii=ea_load_nii(tfidcell{tc});
                nii.img=nii.img./max(nii.img(:));
                AllX=[AllX+nii.img];
            end
        end
    end
    nii.img=AllX;
    clear AllX
    nii.img=nii.img./max(nii.img(:));
    nii.fname=fullfile([uipatdir{pt},filesep,'fiducials',filesep,ea_getspace],[fguid,'.nii']);
    nii.dt(1) = 16;
    ea_write_nii(nii);
    gziplocal([uipatdir{pt},filesep,'fiducials',filesep,ea_getspace],[fguid,'.nii']);
    ea_delete(tfidcell);
    
    for tc=1:length(fidcell)
        try
            nii=ea_load_nii(fidcell{tc});
        catch
            keyboard
        end
        nii.img=nii.img./max(nii.img(:));
        if ~exist('AllX','var')
            AllX=nii.img;
        else
            try
                AllX=[AllX+nii.img];
            catch
                ea_conformspaceto(fidcell{1},fidcell{tc},1);
                nii=ea_load_nii(fidcell{tc});
                nii.img=nii.img./max(nii.img(:));
                AllX=[AllX+nii.img];
            end
        end
    end
    nii.img=AllX;
    clear AllX
    nii.img=nii.img./max(nii.img(:));
    nii.fname=fullfile([uipatdir{pt},filesep,'fiducials',filesep,'native'],[fguid,'.nii']);
    nii.dt(1) = 16;
    ea_write_nii(nii);
    
    gziplocal([uipatdir{pt},filesep,'fiducials',filesep,'native'],[fguid,'.nii']);
    ea_delete(fidcell);
end



function gziplocal(pathn,filen)

gzip(fullfile(pathn,filen));
delete(fullfile(pathn,filen));


