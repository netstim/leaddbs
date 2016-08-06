function varargout=ea_normalize_maget(options)

if ischar(options) % return name of method.
    varargout{1}='MAGeT Brain-like Segmentation/Normalization DISTAL atlas (Chakravarty 2013, Ewert 2016)';
    varargout{2}={'SPM8','SPM12'};
    return
end


atlastouse='DISTAL_manual';

peerfolders={'/Volumes/Neuro_Charite/ACPC/mni-hisub25/s01'
'/Volumes/Neuro_Charite/ACPC/mni-hisub25/s02'
'/Volumes/Neuro_Charite/ACPC/mni-hisub25/s03'
'/Volumes/Neuro_Charite/ACPC/mni-hisub25/s04'
'/Volumes/Neuro_Charite/ACPC/mni-hisub25/s05'
'/Volumes/Neuro_Charite/ACPC/mni-hisub25/s06'
'/Volumes/Neuro_Charite/ACPC/mni-hisub25/s07'
'/Volumes/Neuro_Charite/ACPC/mni-hisub25/s08'
'/Volumes/Neuro_Charite/ACPC/mni-hisub25/s09'
'/Volumes/Neuro_Charite/ACPC/mni-hisub25/s10'
'/Volumes/Neuro_Charite/ACPC/mni-hisub25/s11'
'/Volumes/Neuro_Charite/ACPC/mni-hisub25/s12'
'/Volumes/Neuro_Charite/ACPC/mni-hisub25/s13'
'/Volumes/Neuro_Charite/ACPC/mni-hisub25/s14'
'/Volumes/Neuro_Charite/ACPC/mni-hisub25/s15'
'/Volumes/Neuro_Charite/ACPC/mni-hisub25/s16'
'/Volumes/Neuro_Charite/ACPC/mni-hisub25/s17'
'/Volumes/Neuro_Charite/ACPC/mni-hisub25/s18'
'/Volumes/Neuro_Charite/ACPC/mni-hisub25/s19'
'/Volumes/Neuro_Charite/ACPC/mni-hisub25/s20'
'/Volumes/Neuro_Charite/ACPC/mni-hisub25/s21'
'/Volumes/Neuro_Charite/ACPC/mni-hisub25/s22'
'/Volumes/Neuro_Charite/ACPC/mni-hisub25/s23'
'/Volumes/Neuro_Charite/ACPC/mni-hisub25/s24'
'/Volumes/Neuro_Charite/ACPC/mni-hisub25/s25'};

%% step 0: check if all subjects have been processed with ANTs multimodal or ANTs
for peer=1:length(peerfolders)
   if ~ismember(ea_whichnormmethod([peerfolders{peer},filesep]),ea_getantsnormfuns)
      ea_error('Please make sure that all peers selected have been normalized using ANTs.')
   end
end

subdirec=[options.root,options.patientname,filesep];
if ~ismember(ea_whichnormmethod([peerfolders{peer},filesep]),ea_getantsnormfuns)
   ea_normalize_ants_multimodal(options); 
end

%% step 1, warp DISTAL back to each peer brain
earoot=ea_getearoot;
atlasbase=[earoot,'atlases',filesep,atlastouse,filesep];
for peer=1:length(peerfolders)

    peerdirec=[peerfolders{peer},filesep];
    poptions=options;
    [poptions.root,poptions.patientname]=fileparts(peerfolders{peer});
    poptions.root=[poptions.root,filesep];
    poptions=ea_assignpretra(poptions);

        
        if ~exist([peerdirec,'MAGeT'],'file') % set up folder structure
            mkdir([peerdirec,'MAGeT']);
            mkdir([peerdirec,'MAGeT',filesep,'atlases']);
            mkdir([peerdirec,'MAGeT',filesep,'atlases',filesep,atlastouse]);
            mkdir([peerdirec,'MAGeT',filesep,'atlases',filesep,atlastouse,filesep,'lh']);
            mkdir([peerdirec,'MAGeT',filesep,'atlases',filesep,atlastouse,filesep,'rh']);
            mkdir([peerdirec,'MAGeT',filesep,'atlases',filesep,atlastouse,filesep,'midline']);
            mkdir([peerdirec,'MAGeT',filesep,'atlases',filesep,atlastouse,filesep,'mixed']);
            mkdir([peerdirec,'MAGeT']);
            mkdir([peerdirec,'MAGeT',filesep,'atlasreceives']);
            mkdir([peerdirec,'MAGeT',filesep,'atlasreceives',filesep,poptions.patientname]);
            mkdir([peerdirec,'MAGeT',filesep,'atlasreceives',filesep,poptions.patientname,filesep,atlastouse]);
            mkdir([peerdirec,'MAGeT',filesep,'atlasreceives',filesep,poptions.patientname,filesep,atlastouse,filesep,'lh']);
            mkdir([peerdirec,'MAGeT',filesep,'atlasreceives',filesep,poptions.patientname,filesep,atlastouse,filesep,'rh']);
            mkdir([peerdirec,'MAGeT',filesep,'atlasreceives',filesep,poptions.patientname,filesep,atlastouse,filesep,'midline']);
            mkdir([peerdirec,'MAGeT',filesep,'atlasreceives',filesep,poptions.patientname,filesep,atlastouse,filesep,'mixed']);
        end
        oatlasbase=[peerdirec,'MAGeT',filesep,'atlases',filesep,atlastouse,filesep];
        satlasbase=[subdirec,'MAGeT',filesep,'atlasreceives',filesep,poptions.patientname,filesep,atlastouse,filesep];

        load([atlasbase,'atlas_index.mat']);
        cnt=1;
        for atlas=1:length(atlases.names)
            % warp atlas to peer

            switch atlases.types(atlas)
                case 1 % LH only
                    froms{cnt}=[atlasbase,'lh',filesep,atlases.names{atlas}];
                    tos{cnt}=[oatlasbase,'lh',filesep,atlases.names{atlas}];
                    sub_tos{cnt}=[satlasbase,'lh',filesep,atlases.names{atlas}];
                    
                    cnt=cnt+1;
                case 2 % RH
                    froms{cnt}=[atlasbase,'rh',filesep,atlases.names{atlas}];
                    tos{cnt}=[oatlasbase,'rh',filesep,atlases.names{atlas}];
                    sub_tos{cnt}=[satlasbase,'rh',filesep,atlases.names{atlas}];

                    cnt=cnt+1;
                case 3 % both RH / LH present
                    froms{cnt}=[atlasbase,'lh',filesep,atlases.names{atlas}];
                    tos{cnt}=[oatlasbase,'lh',filesep,atlases.names{atlas}];
                    sub_tos{cnt}=[satlasbase,'lh',filesep,atlases.names{atlas}];

                    froms{cnt+1}=[atlasbase,'rh',filesep,atlases.names{atlas}];
                    tos{cnt+1}=[oatlasbase,'rh',filesep,atlases.names{atlas}];
                    sub_tos{cnt+1}=[satlasbase,'rh',filesep,atlases.names{atlas}];

                    cnt=cnt+2;
                case 4 % Mixed
                    froms{cnt}=[atlasbase,'mixed',filesep,atlases.names{atlas}];
                    tos{cnt}=[oatlasbase,'mixed',filesep,atlases.names{atlas}];
                    sub_tos{cnt}=[satlasbase,'mixed',filesep,atlases.names{atlas}];

                    cnt=cnt+1;
                case 5 % Midline
                    froms{cnt}=[atlasbase,'midline',filesep,atlases.names{atlas}];
                    tos{cnt}=[oatlasbase,'midline',filesep,atlases.names{atlas}];
                    sub_tos{cnt}=[satlasbase,'midline',filesep,atlases.names{atlas}];

                    cnt=cnt+1;
            end
            
        end
        
        % check for .gz files:
        for fi=1:length(froms)
            if ~exist(froms{fi},'file') && exist([froms{fi},'.gz'],'file')
                froms{fi}=[froms{fi},'.gz'];
                tos{fi}=[tos{fi},'.gz']; % stay consistent here
                sub_tos{fi}=[sub_tos{fi},'.gz'];
            else
                ea_error(['Atlas file not found or dual .nii.gz/.nii version present: ',froms{fi}]);
            end
        end
        
        if ~exist([peerdirec,'MAGeT',filesep,'atlases',filesep,atlastouse],'file') % assume the work has been done already
            
            ea_ants_applytransforms(poptions,froms,tos,1,[peerdirec,poptions.prefs.prenii_unnormalized]);
            
        end
    
    
    
    %% step 2, generate warps from peers to the selected patient brain

    
    if ~exist([subdirec,'MAGeT',filesep,'warpreceives',filesep,poptions.patientname,'Composite.h5'],'file')
    
    
        % build up tos and froms

        [~,peerpresentfiles]=ea_assignpretra(poptions);
        [~,subpresentfiles]=ea_assignpretra(options);
        

        presentinboth=ismember(subpresentfiles,peerpresentfiles);
        peerpresentfiles=peerpresentfiles(presentinboth);
        if ~isequal(subpresentfiles,peerpresentfiles) % then I made something wrong.
            keyboard
        else
           presentfiles=subpresentfiles;
           clear subpresentfiles peerpresentfiles
        end
        for anatfi=1:length(presentfiles)
            spfroms{anatfi}=[peerdirec,presentfiles{anatfi}];
            sptos{anatfi}=[subdirec,presentfiles{anatfi}];
            metrics{anatfi}='MI';
        end
        weights=repmat(1.5,length(presentfiles),1);
       
        % add FA if present ? add to beginning since will use last entry to
        % converge
        if exist([subdirec,options.prefs.fa2anat],'file') && exist([peerdirec,options.prefs.fa2anat],'file')
           spfroms=[{[peerdirec,options.prefs.fa2anat]};spfroms];
           sptos=[{[subdirec,options.prefs.fa2anat]};spfroms];
           weights=[1,weights];
           metrics=[{'CC'},metrics];
        end
        
        defoutput=[subdirec,'MAGeT',filesep,'Warpreceives',filesep,poptions.patientname];
        if ~exist([subdirec,'MAGeT',filesep,'Warpreceives',filesep],'file')
            mkdir([subdirec,'MAGeT',filesep,'Warpreceives',filesep]);
        end
        
        
        ea_ants_nonlinear(sptos,spfroms,[subdirec,'MAGeT',filesep,'Warpreceives',filesep,poptions.patientname,'.nii'],weights,metrics,options);
        delete([subdirec,'MAGeT',filesep,'Warpreceives',filesep,poptions.patientname,'.nii']); % we only need the warp
        delete([subdirec,'MAGeT',filesep,'Warpreceives',filesep,poptions.patientname,'InverseComposite.h5']); % we dont need the inverse warp
    end
    
    
    %% step 3: warp all atlas nuclei from peers to sub
    
    %if ~exist([peerdirec,'MAGeT',filesep,'atlases',filesep,atlastouse],'file') % assume the work has been done already
        % caution: old tos is new from here! this is correct.
       
       transformfile=[subdirec,'MAGeT',filesep,'Warpreceives',filesep,poptions.patientname,filesep,'Composite.h5'];
        ea_ants_applytransforms(poptions,tos,sub_tos,1,[peerdirec,poptions.prefs.prenii_unnormalized],transformfile);
        
    %end
    
    
end



