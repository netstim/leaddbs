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

directory=[options.root,options.patientname,filesep];
if ~ismember(ea_whichnormmethod([peerfolders{peer},filesep]),ea_getantsnormfuns)
   ea_normalize_ants_multimodal(options); 
end

%% step 1, warp DISTAL back to each peer brain
earoot=ea_getearoot;
atlasbase=[earoot,'atlases',filesep,atlastouse,filesep];
for peer=1:length(peerfolders)

    peerdirec=[peerfolders{peer},filesep];
    
    if ~exist([peerdirec,'MAGeT'],'file')
        mkdir([peerdirec,'MAGeT']);
        mkdir([peerdirec,'MAGeT',filesep,'atlases']);
        mkdir([peerdirec,'MAGeT',filesep,'atlases',filesep,atlastouse]);
        mkdir([peerdirec,'MAGeT',filesep,'atlases',filesep,atlastouse,filesep,'lh']);
        mkdir([peerdirec,'MAGeT',filesep,'atlases',filesep,atlastouse,filesep,'rh']);
        mkdir([peerdirec,'MAGeT',filesep,'atlases',filesep,atlastouse,filesep,'midline']);
        mkdir([peerdirec,'MAGeT',filesep,'atlases',filesep,atlastouse,filesep,'mixed']);
        
    end
    oatlasbase=[peerdirec,'MAGeT',filesep,'atlases',filesep,atlastouse,filesep];
    load([atlasbase,'atlas_index.mat']);
    cnt=1;
    for atlas=1:length(atlases.names)
        % warp atlas to peer
        poptions=options;
        [poptions.root,poptions.patientname]=fileparts(peerfolders{peer});
        poptions.root=[poptions.root,filesep];
        switch atlases.types(atlas)
            case 1 % LH only
                froms{cnt}=[atlasbase,'lh',filesep,atlases.names{atlas}];
                tos{cnt}=[oatlasbase,'lh',filesep,atlases.names{atlas}];
                cnt=cnt+1;
            case 2 % RH
                froms{cnt}=[atlasbase,'rh',filesep,atlases.names{atlas}];
                tos{cnt}=[oatlasbase,'rh',filesep,atlases.names{atlas}];
                cnt=cnt+1;
            case 3 % both RH / LH present
                froms{cnt}=[atlasbase,'lh',filesep,atlases.names{atlas}];
                tos{cnt}=[oatlasbase,'lh',filesep,atlases.names{atlas}];
                froms{cnt+1}=[atlasbase,'rh',filesep,atlases.names{atlas}];
                tos{cnt+1}=[oatlasbase,'rh',filesep,atlases.names{atlas}];
                cnt=cnt+2;
            case 4 % Mixed
                froms{cnt}=[atlasbase,'mixed',filesep,atlases.names{atlas}];
                tos{cnt}=[oatlasbase,'mixed',filesep,atlases.names{atlas}];
                cnt=cnt+1;
            case 5 % Midline
                froms{cnt}=[atlasbase,'midline',filesep,atlases.names{atlas}];
                tos{cnt}=[oatlasbase,'midline',filesep,atlases.names{atlas}];
                cnt=cnt+1;
        end
        poptions=ea_assignpretra(poptions);
    ea_ants_applytransforms(poptions,froms,tos,1,[peerdirec,poptions.prefs.prenii_unnormalized]);

    keyboard
    end
end



