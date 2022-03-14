function ea_importspaceassets(~,~,fromspace,what,infname,outfname)


if strcmp(what,'custom')
    if ~exist('infname','var')
        [infi,inpth]=uigetfile('*.nii','Choose file to warp...');
        infname=fullfile(inpth,infi);
    end
end

ea_genwarp2space(fromspace);


switch what
    case 'both'
        ea_warplabelassets(fromspace);
        ea_warpatlasassets(fromspace);
    case 'atlases'
        ea_warpatlasassets(fromspace);

    case 'labeling'
        ea_warplabelassets(fromspace);

    case 'custom'
        ea_warpfilefromspace(fromspace,infname);

end
% if ~strcmp(what,'custom')
% rmdir([ea_space,fromspace],'s'); % cleanup all.
% end

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
   to{1}=[ea_space([],'labeling'),parcellation{p},'.nii'];
   if ~exist(ea_space([],'labeling'),'dir')
       mkdir(ea_space([],'labeling'));
   end
       if exist(to{1},'file')
          disp(['A whole-brain parcellation with the same name (',to{1},') already exists! Skipping!']);
          continue
       else
           directory=[ea_space,fromspace,filesep];
            options.prefs=ea_prefs('');
            ea_apply_normalization_tofile(options,from,to,0,0);
       end

       movefile(to{1},[ea_space([],'labeling'),parcellation{p},' - imported from ',fromspace,'.nii']);
       try
           movefile([wlabeling,parcellation{p},'.txt'],[ea_space([],'labeling'),parcellation{p},' - imported from ',fromspace,'.txt']);
       catch
           warning([wlabeling,parcellation{p},'.txt', ' not present. Labeling folder in source space seems inconsistent. Moving on.']);
       end

end


function ea_warpfilefromspace(fromspace,infname)

[pth,inf,ext]=fileparts(infname);
outfname=fullfile(pth,['w',inf,ext]);

directory=[ea_space,fromspace,filesep];
options=ea_getptopts(directory);
options.prefs=ea_prefs('');
options=ea_assignpretra(options);
ea_apply_normalization_tofile(options,{infname},{outfname},0,4,[ea_space,options.primarytemplate,'.nii']);






function ea_warpatlasassets(fromspace)
directory=[ea_space,fromspace,filesep];
options=ea_getptopts(directory);
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
    if exist([ea_space([],'atlases'),asc{atlasset},' - imported from ',fromspace],'dir')
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

    % for .mat case
    src = fullfile(ea_space,fromspace,'anat_t1.nii');
    invt = fullfile(ea_space,fromspace,'inverseTransform');
    srcV = ea_open_vol(src);

    atlroot=[ea_space,fromspace,filesep,asc{atlasset},filesep];
    load([atlroot,'atlas_index.mat'])

    for atlas=1:length(atlases.names)

        [~,~,ext] = fileparts(atlases.names{atlas});

        switch ext
            case {'.nii', '.gz'}
                wasgzip = strcmp(ext,'.gz');

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
                        if ~strcmp(nii{n}(end-2:end),'.gz')
                            gunzip([nii{n},'.gz']);
                        else
                            gunzip([nii{n}]);
                        end
                        nii{n}=strrep(nii{n},'.gz','');
                    end
                    from{1}=nii{n}; to{1}=nii{n};
                    directory=[ea_space,fromspace,filesep];
                    options.prefs=ea_prefs('');
                    ea_apply_normalization_tofile(options,from,to,0,1);
                    if wasgzip
                        gzip(nii{n});
                        delete(nii{n});
                    end
                end

            case '.mat'
                atlas_full_dir = dir(fullfile(atlroot,'*',atlases.names{atlas}));
                for n = 1:length(atlas_full_dir)
                    atl_name = fullfile(atlas_full_dir(n).folder,atlas_full_dir(n).name);
                    atl_load = load(atl_name);
                    fibs_mm = [atl_load.fibers(:,1:3)'; ones(1,size(atl_load.fibers,1))];
                    fibs_vox = srcV.mat \ fibs_mm;
                    XYZ_dest_mm = ea_map_coords(fibs_vox, src, invt);
                    atl_load.fibers(:,1:3) = XYZ_dest_mm';
                    save(atl_name, '-struct', 'atl_load');
                end

        end

    end
    % finally rebuild index:
    % and move to atlases dir
    delete([atlroot,'atlas_index.mat']);
if ~exist(ea_space([],'atlases'),'dir')
    mkdir(ea_space([],'atlases'));
end
    movefile(atlroot,[ea_space([],'atlases'),asc{atlasset},' - imported from ',fromspace]);
    options.atlasset=[asc{atlasset},' - imported from ',fromspace];
    ea_genatlastable([],ea_space([],'atlases'),options);


end
