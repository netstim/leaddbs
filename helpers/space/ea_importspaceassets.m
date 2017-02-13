function ea_importspaceassets(~,~,fromspace,what,infname,outfname)

ea_genwarp2space(fromspace);
norm_method_applied{1}='ea_normalize_spmdartel';
save([ea_space,fromspace,filesep,'ea_normmethod_applied.mat'],'norm_method_applied');

switch what
    case 'both'
        ea_warplabelassets(fromspace);
        ea_warpatlasassets(fromspace);
    case 'atlases'
        ea_warpatlasassets(fromspace);
        
    case 'labeling'
        ea_warplabelassets(fromspace);
        
    case 'custom'
        ea_warpfilefromspace(fromspace,infname,outfname);

end
if ~strcmp(what,'custom')
rmdir([ea_space,fromspace],'s'); % cleanup all.
end

function ea_warplabelassets(fromspace)
labelinginotherspace=[ea_getearoot,'templates',filesep,'space',filesep,fromspace,filesep,'labeling',filesep];
if exist(labelinginotherspace,'dir')
   copyfile(fileparts(labelinginotherspace),[ea_space,fromspace,filesep,'labeling']); 
else
    disp('No whole-brain parcellations found.')
    return
end

wlabeling=[ea_space,fromspace,filesep,'labeling',filesep];
ll=dir([wlabeling,'*.nii']);
if isempty(ll)
    disp('No whole-brain parcellations found.')
    rmdir(wlabeling,'s');
    return
end
parcellation={};
for lab=1:length(ll)
    [~,n]=fileparts(ll(lab).name);
    parcellation{lab}=n;
end

for p=1:length(parcellation)
   from{1}=[wlabeling,parcellation{p},'.nii'];
   to{1}=[ea_space([],'labeling'),parcellation{p},' [imported from ',fromspace,']','.nii'];
       if exist(to{1},'file')
          disp(['A whole-brain parcellation with the same name (',to{1},') already exists! Skipping!']);
          continue
       else
           directory=[ea_space,fromspace,filesep];
            options.prefs=ea_prefs('');
            ea_apply_normalization_tofile(options,from,to,directory,0,0);
       end
    movefile([wlabeling,parcellation{p},'.txt'],[ea_space([],'labeling'),parcellation{p},' [imported from ',fromspace,']','.txt']);

    
end


function ea_warpfilefromspace(fromspace,infname,outfname)


options.prefs=ea_prefs('');
directory=[ea_space,fromspace,filesep];
ea_apply_normalization_tofile(options,{infname},{outfname},directory,0,0);






function ea_warpatlasassets(fromspace)

options.prefs=ea_prefs('');
% list atlases:
as=dir([ea_getearoot,'templates',filesep,'space',filesep,fromspace,filesep,'atlases']);

asc=cell(0);
cnt=1;
for i=1:length(as)
    if as(i).isdir && ~strcmp(as(i).name(1),'.')
        asc{cnt}=as(i).name;
        cnt=cnt+1;
    end
end


for atlasset=1:length(asc)
    if exist([ea_space([],'atlases'),asc{atlasset},' [imported from ',fromspace,']'],'dir')
        disp([ea_space([],'atlases'),asc{atlasset},' already exists! Skipping...'])
        continue
    else
        copyfile([ea_getearoot,'templates',filesep,'space',filesep,fromspace,filesep,'atlases',filesep,asc{atlasset}],[ea_space,fromspace,filesep,asc{atlasset}]);
    end
    if ~exist([ea_space,fromspace,filesep,asc{atlasset},filesep,'atlas_index.mat'],'file')
        disp([asc{atlasset},' needs to be visualized in original space before import! Skipping.']);
        rmdir([ea_space,fromspace,filesep,asc{atlasset}],'s');
        continue
    end
    
    atlroot=[ea_space,fromspace,filesep,asc{atlasset},filesep];
    
    load([atlroot,'atlas_index.mat'])
    for atlas=1:length(atlases.names)
        [~,~,ext]=fileparts(atlases.names{atlas});
        if strcmp(ext,'.gz')
            wasgzip=1;
        else
            wasgzip=0;
        end
        
        switch atlases.types(atlas)
            case 1 % right hemispheric atlas.
                nii{1}=[atlroot,'rh',filesep,atlases.names{atlas}];
            case 2 % left hemispheric atlas.
                nii{1}=[atlroot,'lh',filesep,atlases.names{atlas}];
            case 3 % both-sides atlas composed of 2 files.
                nii{1}=[atlroot,'lh',filesep,atlases.names{atlas}];
                nii{2}=[atlroot,'rh',filesep,atlases.names{atlas}];
            case 4 % mixed atlas (one file with both sides information).
                nii{1}=[atlroot,'mixed',filesep,atlases.names{atlas}];
                
            case 5 % midline atlas (one file with both sides information.
                nii{1}=[atlroot,'midline',filesep,atlases.names{atlas}];
        end
        
        for n=1:length(nii)
            if wasgzip
                gunzip(nii{n});
                nii{n}=strrep(nii{n},'.gz','');
            end
            from{1}=nii{n}; to{1}=nii{n};
            directory=[ea_space,fromspace,filesep];
            options.prefs=ea_prefs('');
            ea_apply_normalization_tofile(options,from,to,directory,0,1);
            if wasgzip
                gzip(nii{n});
                delete(nii{n});
            end
        end

    end
    % finally rebuild index:
    % and move to atlases dir
    delete([atlroot,'atlas_index.mat']);
    
    movefile(atlroot,[ea_space([],'atlases'),asc{atlasset},' [imported from ',fromspace,']']);
    options.atlasset=[asc{atlasset},' [imported from ',fromspace,']'];
    ea_genatlastable([],fileparts([ea_space([],'atlases')]),options);
    
    
end
