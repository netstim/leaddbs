function ea_exportisovolume(elstruct,options)
% this function exports an isovolume to a nifti file. It is largely based
% on ea_showisovolume
% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn


disp('*** Exporting isovolume to nifti files.');

if size(options.d3.isomatrix{1},2)==get_maxNumContacts(elstruct)-1 % number of contact pairs
    shifthalfup=1;
elseif size(options.d3.isomatrix{1},2)==get_maxNumContacts(elstruct) % number of contacts
    shifthalfup=0;
else
    ea_cprintf('CmdWinErrors', 'Be careful! Isomatrix might have wrong size, or numbers of contacts are not consistent across patients.\n');
end

for iside=1:length(options.sides)
    side=options.sides(iside);

    cnt=1;
    for sub=1:length(elstruct)
        for cont=1:size(options.d3.isomatrix{1},2)
            if ~isnan(options.d3.isomatrix{side}(sub,cont))
                 if ea_arenopoints4side(elstruct(sub).coords_mm,side)
                    %if there are no coordinates, don't parse anything, and skip to next
                    warning_printf=@(str_in) fprintf(['ATTENTION!! : ' str_in '\n']);
                    if side==1
                        warning_printf(['no isovolume will be exported for the right side of subj #' num2str(sub) ' as there is no lead in it.']);
                    elseif side==2
                        warning_printf(['no isovolume will be exported for the right side of subj #' num2str(sub) ' as there is no lead in it.']);
                    else
                        warning_printf(['no isovolume will be exported for side=' num2str(side) ' of subj #' num2str(sub) ' as there is no lead in it.']);
                    end
                 else
                    %if there are coordinates, parse them, otherwise skip to next
                    try
                        if ~shifthalfup
                            X{side}(cnt)=elstruct(sub).coords_mm{side}(cont,1);
                            Y{side}(cnt)=elstruct(sub).coords_mm{side}(cont,2);
                            Z{side}(cnt)=elstruct(sub).coords_mm{side}(cont,3);
                        else % using pairs of electrode contacts (i.e. 3 pairs if there are 4 contacts)
                            X{side}(cnt)=mean([elstruct(sub).coords_mm{side}(cont,1),elstruct(sub).coords_mm{side}(cont+1,1)]);
                            Y{side}(cnt)=mean([elstruct(sub).coords_mm{side}(cont,2),elstruct(sub).coords_mm{side}(cont+1,2)]);
                            Z{side}(cnt)=mean([elstruct(sub).coords_mm{side}(cont,3),elstruct(sub).coords_mm{side}(cont+1,3)]);
                        end
                        V{side}(cnt)=options.d3.isomatrix{side}(sub,cont);

                        PT{side}(cnt)=sub;
                        cnt=cnt+1;
                    catch
                        warning(['Please check localization for subject no. ',num2str(sub),'.']);
                    end
                 end
            end
        end
    end

    if cnt==1
        %no elements/leads were present in this side, skip it
        continue
    end

    X{side}=X{side}(:);
    Y{side}=Y{side}(:);
    Z{side}=Z{side}(:);
    V{side}=V{side}(:);
    PT{side}=PT{side}(:);

    Vol=spm_vol([ea_space,'bb.nii']);
    nii{side}=spm_read_vols(Vol);
    nii{side}(:)=nan;
    
    XYZ=[X{side},Y{side},Z{side},ones(length(X{side}),1)]';
    
    XYZ=Vol.mat\XYZ; % to voxel space.
    XYZ=(XYZ(1:3,:)');
    % repeat the above but in voxel space..
    clear bb
    bb(1,:)=[round(min(XYZ(:,1))),round(max(XYZ(:,1)))];
    bb(2,:)=[round(min(XYZ(:,2))),round(max(XYZ(:,2)))];
    bb(3,:)=[round(min(XYZ(:,3))),round(max(XYZ(:,3)))];
    clear XI YI ZI
    [XI,YI,ZI]=meshgrid([bb(1,1):bb(1,2)],[bb(2,1):bb(2,2)],[bb(3,1):bb(3,2)]);
    warning('off')

    nanix=~isnan(V{side});
    F = scatteredInterpolant(XYZ(nanix,1),XYZ(nanix,2),XYZ(nanix,3),double(V{side}(nanix)),'natural');
    warning('on')
    F.ExtrapolationMethod='none';

    %    [p,idx]=ea_isosignificance([XYZ,double([V{side}])],1,0.5);

    xix{side}=bb(1,1):bb(1,2); yix{side}=bb(2,1):bb(2,2); zix{side}=bb(3,1):bb(3,2);

    nii{side}(xix{side},yix{side},zix{side})=F({xix{side},yix{side},zix{side}});
    nii2{side}=nii{side};

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
    %         spm_jobman('run',jobs);
    %         clear jobs matlabbatch

    if side==2 % write out combined volume with separate interpolations for each side.
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
        
        spm_jobman('run',jobs);
        clear jobs matlabbatch
        %% Also write out volume with combined information on both sides (symmetric image).

        niic=ea_load_nii([ea_space,'bb.nii']);
        niic.dt(1) = 16;

        %niic=spm_read_vols(Vol);
        niic.img(:)=nan;
        niicsig=niic.img;
        for inside=1:2
            switch inside
                case 1 % flip infos from right to left
                    
                    XYZtf=[X{1},Y{1},Z{1}];
                    XYZflip=ea_flip_lr_nonlinear(XYZtf);
                    
                    XYZ=[[XYZflip(:,1);X{2}],[Y{1};Y{2}],[Z{1};Z{2}],ones(length([X{1};X{2}]),1)]';
                case 2 % flip infos from left to right
                    XYZtf=[X{2},Y{2},Z{2}];
                    XYZflip=ea_flip_lr_nonlinear(XYZtf);
                    
                    XYZ=[[X{1};XYZflip(:,1)],[Y{1};Y{2}],[Z{1};Z{2}],ones(length([X{1};X{2}]),1)]';
            end
            XYZ=niic.mat\XYZ; % to voxel space.
            XYZ=(XYZ(1:3,:)');
            % repeat the above but in voxel space..
            clear bb
            bb(1,:)=[round(min(XYZ(:,1))),round(max(XYZ(:,1)))];
            bb(2,:)=[round(min(XYZ(:,2))),round(max(XYZ(:,2)))];
            bb(3,:)=[round(min(XYZ(:,3))),round(max(XYZ(:,3)))];
            clear XI YI ZI
            [XI,YI,ZI]=meshgrid([bb(1,1):bb(1,2)],[bb(2,1):bb(2,2)],[bb(3,1):bb(3,2)]);

            warning('off');

            nanix=[(~isnan(V{1}));(~isnan(V{2}))];
            AllV=[V{1};V{2}];
            F = scatteredInterpolant(XYZ(nanix,1),XYZ(nanix,2),XYZ(nanix,3),double(AllV(nanix)),'natural');

            F.ExtrapolationMethod='none';
            warning('on');

            xixc=bb(1,1):bb(1,2); yixc=bb(2,1):bb(2,2); zixc=bb(3,1):bb(3,2);

            niic.img(xixc,yixc,zixc)=F({xixc,yixc,zixc});

            %% write out significant volume:
            if options.d2.write % only needs to be done once..
                XYZV=[XYZ,[V{1};V{2}]];
                PTb=[PT{1};PT{2}];
                if inside==1
                    significancemethod=0;
                    switch significancemethod
                        case 1 % estimate significance by estimating centrality
                            [ixes]=ea_centrality_significance(XYZV);
                            if sum(ixes)>3
                                XYZV=XYZV(ixes,:); % only significant entries..
                                XYZV(:,4)=1;
                                warning('off');
                                Fsig = scatteredInterpolant(XYZV(:,1),XYZV(:,2),XYZV(:,3),XYZV(:,4));

                                Fsig.ExtrapolationMethod='none';
                                warning('on');


                                niicsig(xixc,yixc,zixc)=Fsig({xixc,yixc,zixc});
                            end
                        case 2 % estimate significance by smoothing and using SPM
                            niicsig=ea_smooth_significance(XYZV,PTb,niic,niic.img,options);

                        case 3 % estimate significance by applying leave-one-out permutations on scattered interpolant
                            [ixes,R,p]=ea_leoo_significance(XYZV);
                            ea_dumpsigtxt([options.root,options.patientname,filesep,options.d3.isomatrix_name,'_leoo_sig_scatteredinterpolant.txt'],R,p);

                            if sum(ixes)>3
                                XYZV=XYZV(ixes,:); % only significant entries..
                                XYZV(:,4)=1;
                                warning('off');
                                Fsig = scatteredInterpolant(XYZV(:,1),XYZV(:,2),XYZV(:,3),XYZV(:,4));

                                Fsig.ExtrapolationMethod='none';
                                warning('on');


                                niicsig(xixc,yixc,zixc)=Fsig({xixc,yixc,zixc});
                            end
                        case 4 % estimate significance by applying leave-one-out permutations on weighted average (by distance) from rest of data
                            [ixes,R,p]=ea_leoo_significance_weightedave(XYZV);
                            ea_dumpsigtxt([options.root,options.patientname,filesep,options.d3.isomatrix_name,'_leoo_sig_weightedave.txt'],R,p);
                            if sum(ixes)>3
                                XYZV=XYZV(ixes,:); % only significant entries..
                                XYZV(:,4)=1;
                                warning('off');
                                Fsig = scatteredInterpolant(XYZV(:,1),XYZV(:,2),XYZV(:,3),XYZV(:,4));

                                Fsig.ExtrapolationMethod='none';
                                warning('on');
                                niicsig(xixc,yixc,zixc)=Fsig({xixc,yixc,zixc});
                            end

                        case 5 % estimate significance by applying leave-one-out permutations on weighted distance (by value) from rest of data
                            [ixes,R,p]=ea_leoo_significance_weighteddist(XYZV);
                            ea_dumpsigtxt([options.root,options.patientname,filesep,options.d3.isomatrix_name,'_leoo_sig_weighteddist.txt'],R,p);

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
            end
        end

        niic.fname=[options.root,options.patientname,filesep,options.d3.isomatrix_name,'_combined.nii'];
        ea_write_nii(niic);

        %ea_crop_nii([options.root,options.patientname,filesep,options.d3.isomatrix_name,'_combined.nii'],'','nn');
        % smooth image.

        clear jobs matlabbatch

        %% write out significant volume:

        try
            if any(niicsig(:))
                niic.fname=[options.root,options.patientname,filesep,options.d3.isomatrix_name,'_combined_p05.nii'];
                niic.img=niicsig;
                ea_write_nii(niic);
                ea_crop_nii([options.root,options.patientname,filesep,'s',options.d3.isomatrix_name,'_combined_p05.nii']);
            end
        end

        ea_crop_nii([options.root,options.patientname,filesep,options.d3.isomatrix_name,'_lr.nii'],'','nz',1,1);
        ea_crop_nii([options.root,options.patientname,filesep,'s',options.d3.isomatrix_name,'_lr.nii'],'','nz',1,1);
        ea_crop_nii([options.root,options.patientname,filesep,options.d3.isomatrix_name,'_combined.nii'],'','nz',1,1);

        matlabbatch{1}.spm.spatial.smooth.data = {[options.root,options.patientname,filesep,options.d3.isomatrix_name,'_combined.nii,1']};
        matlabbatch{1}.spm.spatial.smooth.fwhm = [1 1 1];
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 1;
        matlabbatch{1}.spm.spatial.smooth.prefix = 's';
        jobs{1}=matlabbatch;
        spm_jobman('run',jobs);
    end
end

disp('*** Done exporting isovolume to nifti files.');


function ea_dumpsigtxt(filepath,R,p)
fid=fopen(filepath,'w');

if p<0.05 && R>0
    fprintf(fid,['Found significant positive relationship in data without permutations (R=',num2str(R),', p=',num2str(p),').']);
else
    fprintf(fid,['No significant positive relationship in data found (R=',num2str(R),', p=',num2str(p),').']);
end
fclose(fid);


function maxNumContacts = get_maxNumContacts(elstruct)
coords = {elstruct.coords_mm};
coords = horzcat(coords{:})';
maxNumContacts = max(cellfun(@(x) size(x,1), coords));
