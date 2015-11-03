function ea_exportisovolume(elstruct,options)
% this function exports an isovolume to a nifti file. It is largely based
% on ea_showisovolume
% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn


disp('*** Exporting isovolume to nifti files.');

if size(options.d3.isomatrix{1},2)==4-1 % 3 contact pairs
    shifthalfup=1;
elseif size(options.d3.isomatrix{1},2)==4 % 4 contacts
    shifthalfup=0;
else
    ea_error('Isomatrix has wrong size. Please specify a correct matrix.')
end
for side=options.sides
    
    cnt=1;
    for sub=1:length(elstruct)
        for cont=1:size(options.d3.isomatrix{1},2)
            if ~isnan(options.d3.isomatrix{side}(sub,cont));
                if ~shifthalfup
                    X{side}(cnt)=elstruct(sub).coords_mm{side}(cont,1);
                    Y{side}(cnt)=elstruct(sub).coords_mm{side}(cont,2);
                    Z{side}(cnt)=elstruct(sub).coords_mm{side}(cont,3);
                else % using pairs of electrode contacts (i.e. 3 pairs if there are 4 contacts)
                    try
                        X{side}(cnt)=mean([elstruct(sub).coords_mm{side}(cont,1),elstruct(sub).coords_mm{side}(cont+1,1)]);
                        Y{side}(cnt)=mean([elstruct(sub).coords_mm{side}(cont,2),elstruct(sub).coords_mm{side}(cont+1,2)]);
                        Z{side}(cnt)=mean([elstruct(sub).coords_mm{side}(cont,3),elstruct(sub).coords_mm{side}(cont+1,3)]);
                    catch
                        ea_error(['Please check localization for subject no. ',num2str(sub),'.']);
                    end
                end
                V{side}(cnt)=options.d3.isomatrix{side}(sub,cont);
                
                cnt=cnt+1;
            end
        end
    end
    
    X{side}=X{side}(:);        Y{side}=Y{side}(:);        Z{side}=Z{side}(:); V{side}=V{side}(:);
    
    Vol=spm_vol([options.earoot,'templates',filesep,'bb.nii']);
    nii{side}=spm_read_vols(Vol);
    nii{side}(:)=nan;
    XYZ=[X{side},Y{side},Z{side},ones(length(X{side}),1)]';
    XYZ=Vol.mat\XYZ; % to voxel space.
    XYZ=round(XYZ(1:3,:)');
    % repeat the above but in voxel space..
    clear bb
    bb(1,:)=[min(XYZ(:,1)),max(XYZ(:,1))];
    bb(2,:)=[min(XYZ(:,2)),max(XYZ(:,2))];
    bb(3,:)=[min(XYZ(:,3)),max(XYZ(:,3))];
    clear XI YI ZI
    [XI,YI,ZI]=meshgrid([bb(1,1):bb(1,2)],[bb(2,1):bb(2,2)],[bb(3,1):bb(3,2)]);
    warning('off')
    F = scatteredInterpolant(XYZ(:,1),XYZ(:,2),XYZ(:,3),double(V{side}));
    warning('on')
    F.ExtrapolationMethod='none';
    
    %    [p,idx]=ea_isosignificance([XYZ,double([V{side}])],1,0.5);
    
    xix{side}=bb(1,1):bb(1,2); yix{side}=bb(2,1):bb(2,2); zix{side}=bb(3,1):bb(3,2);
    
    nii{side}(xix{side},yix{side},zix{side})=F({xix{side},yix{side},zix{side}});
    
    
    
    
    
    
    switch side
        case 1
            lr='right';
        case 2
            lr='left';
    end
    %% the following commented code can be used to write out singe-hemisphere volumes
    %         Vol.fname=[options.root,options.patientname,filesep,options.d3.isomatrix_name,'_',lr,'.nii'];
    %         Vol.dtype=[32 1];
    %         spm_write_vol(Vol,nii{side});
    %
    %         matlabbatch{1}.spm.spatial.smooth.data = {[options.root,options.patientname,filesep,options.d3.isomatrix_name,'_',lr,'.nii,1']};
    %         matlabbatch{1}.spm.spatial.smooth.fwhm = [2 2 2];
    %         matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    %         matlabbatch{1}.spm.spatial.smooth.im = 1;
    %         matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    %         jobs{1}=matlabbatch;
    %         cfg_util('run',jobs);
    %         clear jobs matlabbatch
    
    
    
    if side==2; % write out combined volume with separate interpolations for each side.
        %% old part
        Vol.fname=[options.root,options.patientname,filesep,options.d3.isomatrix_name,'_lr.nii'];
        niic=ea_nanmean(cat(4,nii{1},nii{2}),4);
        spm_write_vol(Vol,niic);
        
        %ea_crop_nii([options.root,options.patientname,filesep,options.d3.isomatrix_name,'_lr.nii'],'','nn');
        
        % smooth image.
        
        
        matlabbatch{1}.spm.spatial.smooth.data = {[options.root,options.patientname,filesep,options.d3.isomatrix_name,'_lr.nii,1']};
        matlabbatch{1}.spm.spatial.smooth.fwhm = [0.7 0.7 0.7];
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 1;
        matlabbatch{1}.spm.spatial.smooth.prefix = 's';
        jobs{1}=matlabbatch;
        cfg_util('run',jobs);
        clear jobs matlabbatch
        
        %% Also write out volume with combined information on both sides (symmetric image).
        
        Vol=spm_vol([options.earoot,'templates',filesep,'bb.nii']);
        niic=spm_read_vols(Vol);
        niic(:)=nan;
        niicsig=niic;
        for inside=1:2
            switch inside
                case 1 % flip infos from right to left
                    XYZ=[[-X{1};X{2}],[Y{1};Y{2}],[Z{1};Z{2}],ones(length([X{1};X{2}]),1)]';
                case 2 % flip infos from left to right
                    XYZ=[[X{1};-X{2}],[Y{1};Y{2}],[Z{1};Z{2}],ones(length([X{1};X{2}]),1)]';
            end
            XYZ=Vol.mat\XYZ; % to voxel space.
            XYZ=round(XYZ(1:3,:)');
            % repeat the above but in voxel space..
            clear bb
            bb(1,:)=[min(XYZ(:,1)),max(XYZ(:,1))];
            bb(2,:)=[min(XYZ(:,2)),max(XYZ(:,2))];
            bb(3,:)=[min(XYZ(:,3)),max(XYZ(:,3))];
            clear XI YI ZI
            [XI,YI,ZI]=meshgrid([bb(1,1):bb(1,2)],[bb(2,1):bb(2,2)],[bb(3,1):bb(3,2)]);
            
            
            
            warning('off');
            F = scatteredInterpolant(XYZ(:,1),XYZ(:,2),XYZ(:,3),double([V{1};V{2}]));
            
            F.ExtrapolationMethod='none';
            warning('on');
            
            
            
            
            xixc=bb(1,1):bb(1,2); yixc=bb(2,1):bb(2,2); zixc=bb(3,1):bb(3,2);
            
            niic(xixc,yixc,zixc)=F({xixc,yixc,zixc});
            
            
            %% write out significant volume:
            if options.d2.write % only needs to be done once..
            XYZV=[XYZ,[V{1};V{2}]];
            if inside==1
            [ixes]=ea_eigvcentrality_significance(XYZV);
            end
            
            if sum(ixes)>3
            XYZV=XYZV(ixes,:); % only significant entries..
            XYZV(:,4)=1;
            warning('off');
            Fsig = scatteredInterpolant(XYZV(:,1),XYZV(:,2),XYZV(:,3),XYZV(:,4));
            
            Fsig.ExtrapolationMethod='none';
            warning('on');
            
            
            niicsig(xixc,yixc,zixc)=Fsig({xixc,yixc,zixc});
            end
            end
        end
        
        Vol.fname=[options.root,options.patientname,filesep,options.d3.isomatrix_name,'_combined.nii'];
        spm_write_vol(Vol,niic);
        %ea_crop_nii([options.root,options.patientname,filesep,options.d3.isomatrix_name,'_combined.nii'],'','nn');
        % smooth image.
        matlabbatch{1}.spm.spatial.smooth.data = {[options.root,options.patientname,filesep,options.d3.isomatrix_name,'_combined.nii,1']};
        matlabbatch{1}.spm.spatial.smooth.fwhm = [1 1 1];
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 1;
        matlabbatch{1}.spm.spatial.smooth.prefix = 's';
        jobs{1}=matlabbatch;
        cfg_util('run',jobs);
        clear jobs matlabbatch
        
        %% write out significant volume:
        
        
        Vol.fname=[options.root,options.patientname,filesep,options.d3.isomatrix_name,'_combined_p05.nii'];
        spm_write_vol(Vol,niicsig);

        
        
        
    end
    
    
    
    
    
end

disp('*** Done exporting isovolume to nifti files.');