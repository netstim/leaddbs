function atlases=ea_genatlastable(varargin)
% This function reads in atlases in the Lead-dbs/atlases directory and
% generates a table of all available atlas files.
% Atlastypes:   1 - RH
%               2 - LH
%               3 - both hemispheres (2 files present both in lhs and rhs
%               folder
%               4 - mixed (one file with one cluster on each hemisphere)
%               5 - midline (one file with one cluster in total)
%
%
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

atlases=varargin{1};
root=varargin{2};
if strcmp(root(end), filesep)
    root = fileparts(root);
end

options=varargin{3};
if nargin==4
    mifix=varargin{4};
else
    mifix='';
end

if isempty(atlases) || ~isfield(atlases,'roi') % old format
    disp('Generating Atlas table. This may take a while...');

    lhcell=cell(0); rhcell=cell(0); mixedcell=cell(0); midlinecell=cell(0);

    delete([root,filesep,mifix,options.atlasset,filesep,'lh',filesep,'*_temp.ni*']);
    lhatlases=dir([root,filesep,mifix,options.atlasset,filesep,'lh',filesep,'*.ni*']);
    lhtrks=dir([root,filesep,mifix,options.atlasset,filesep,'lh',filesep,'*.tr*']);
    lhmats=dir([root,filesep,mifix,options.atlasset,filesep,'lh',filesep,'*.mat']);
    lhatlases=[lhatlases;lhtrks;lhmats];
    for i=1:length(lhatlases)
        lhcell{i}=lhatlases(i).name;
    end

    delete([root,filesep,mifix,options.atlasset,filesep,'rh',filesep,'*_temp.ni*']);
    rhatlases=dir([root,filesep,mifix,options.atlasset,filesep,'rh',filesep,'*.ni*']);
    rhtrks=dir([root,filesep,mifix,options.atlasset,filesep,'rh',filesep,'*.tr*']);
    rhmats=dir([root,filesep,mifix,options.atlasset,filesep,'rh',filesep,'*.mat']);
    rhatlases=[rhatlases;rhtrks;rhmats];
    for i=1:length(rhatlases)
        rhcell{i}=rhatlases(i).name;
    end

    delete([root,filesep,mifix,options.atlasset,filesep,'mixed',filesep,'*_temp.ni*']);
    mixedatlases=dir([root,filesep,mifix,options.atlasset,filesep,'mixed',filesep,'*.ni*']);
    mixedtrks=dir([root,filesep,mifix,options.atlasset,filesep,'mixed',filesep,'*.tr*']);
    mixedmats=dir([root,filesep,mifix,options.atlasset,filesep,'mixed',filesep,'*.mat']);
    mixedatlases=[mixedatlases;mixedtrks;mixedmats];
    for i=1:length(mixedatlases)
        mixedcell{i}=mixedatlases(i).name;
    end

    delete([root,filesep,mifix,options.atlasset,filesep,'midline',filesep,'*_temp.ni*']);
    midlineatlases=dir([root,filesep,mifix,options.atlasset,filesep,'midline',filesep,'*.ni*']);
    midlinetrks=dir([root,filesep,mifix,options.atlasset,filesep,'midline',filesep,'*.tr*']);
    midlinemats=dir([root,filesep,mifix,options.atlasset,filesep,'midline',filesep,'*.mat']);
    midlineatlases=[midlineatlases;midlinetrks;midlinemats];
    for i=1:length(midlineatlases)
        midlinecell{i}=midlineatlases(i).name;
    end

    % concatenate lh and rh
    todeletelh=[];
    todeleterh=[];
    for i=1:length(lhcell)
        [ism, loc]=ismember(lhcell{i},rhcell);
        if ism
            todeletelh=[todeletelh,i];
            todeleterh=[todeleterh,loc];
        end
    end

    bothcell=lhcell(todeletelh);
    lhcell(todeletelh)=[];
    rhcell(todeleterh)=[];

    allcell=[rhcell,lhcell,bothcell,mixedcell,midlinecell];

    typecell=[ones(1,length(rhcell)),2*ones(1,length(lhcell)),3*ones(1,length(bothcell)),4*ones(1,length(mixedcell)),5*ones(1,length(midlinecell))];
    atlases.names=allcell;
    atlases.types=typecell;
    if ~isfield(atlases, 'threshold') || isempty(atlases.threshold)
        atlases.threshold.type='relative_intensity';
        atlases.threshold.value=0.5;
    end
end

if ~isfield(atlases,'tissuetypes')
    atlases.tissuetypes=ones(1,length(atlases.names));
end

if checkrebuild(atlases,options,root,mifix)
    %% build iXYZ tables:
    maxcolor=64; % change to 45 to avoid red / 64 to use all colors

    nm=1:2; % native and mni
    try
        nmind=[options.atl.can,options.atl.ptnative]; % which shall be performed?
    catch
        nmind=[1 0];
    end
    nm=nm(logical(nmind)); % select which shall be performed.
    if ~isfield(atlases,'colormap')
        atlases.colormap=ea_color_wes('all');
    end
    for nativemni=nm % switch between native and mni space atlases.
        switch nativemni
            case 1
                root=fileparts(ea_space([],'atlases'));
            case 2
                root=[options.root,options.patientname,filesep,'atlases'];
        end

        atlascnt=1;

        % iterate through atlases, visualize them and write out stats.
        disp('Building atlas table...');
        for atlas=1:length(atlases.names)
            %ea_dispercent(atlas/length(atlases.names));
            switch atlases.types(atlas)
                case 1 % right hemispheric atlas.
                    structure=load_structure([root,filesep,mifix,options.atlasset,filesep,'rh',filesep,atlases.names{atlas}]);
                case 2 % left hemispheric atlas.
                    structure=load_structure([root,filesep,mifix,options.atlasset,filesep,'lh',filesep,atlases.names{atlas}]);
                case 3 % both-sides atlas composed of 2 files.
                    lstructure=load_structure([root,filesep,mifix,options.atlasset,filesep,'lh',filesep,atlases.names{atlas}]);
                    rstructure=load_structure([root,filesep,mifix,options.atlasset,filesep,'rh',filesep,atlases.names{atlas}]);
                case 4 % mixed atlas (one file with both sides information).
                    structure=load_structure([root,filesep,mifix,options.atlasset,filesep,'mixed',filesep,atlases.names{atlas}]);
                case 5 % midline atlas (one file with both sides information).
                    structure=load_structure([root,filesep,mifix,options.atlasset,filesep,'midline',filesep,atlases.names{atlas}]);
            end

            for side=detsides(atlases.types(atlas))
                if atlases.types(atlas)==3 % both-sides atlas composed of 2 files.
                    if side==1
                        structure=rstructure;
                    elseif side==2
                        structure=lstructure;
                    end
                end

                colornames='bgcmywkbgcmywkbgcmywkbgcmywkbgcmywkbgcmywkbgcmywkbgcmywkbgcmywk'; % red is reserved for the VAT.

                colorc=colornames(1);
                colorc=rgb(colorc);
                if isfield(structure, 'img') % volumetric atlas
                        pobj.nii=structure;
                        % set cdata
                        try % check if explicit color info for this atlas is available.
                            cdat=atlases.colormap(atlases.colors(atlas),:);
                        catch
                            cdat=atlases.colormap(round(atlas*(maxcolor/length(atlases.names))),:);
                        end
                        pobj.color=cdat;
                     
                     roi=ea_roi([],pobj);

                    [xx,yy,zz]=ind2sub(size(structure.img),find(structure.img>0)); % find 3D-points that have correct value.
                    vv=structure.img(structure.img(:)>0);

                    if ~isempty(xx)
                        XYZ.vx=[xx,yy,zz]; % concatenate points to one matrix.
                        XYZ.val=vv;
                        XYZ.mm=map_coords_proxy(XYZ.vx,structure); % map to mm-space
                        XYZ.dims=structure.voxsize;
                    else
                        XYZ.vx=[];
                        XYZ.val=[];
                        XYZ.mm=[];
                        XYZ.dims=structure.voxsize;
                    end

                    % surface(xx(1:10)',yy(1:10)',zz(1:10)',ones(10,1)');
                    % hold on

                    if atlases.types(atlas)==4 && side==2 % restore from backup
                        structure=bstructure;
                        XYZ.mm=bXYZ.mm;
                        XYZ.val=bXYZ.val;
                        XYZ.vx=bXYZ.vx;
                    end

                    try
                        bb=[0,0,0;size(structure.img)];
                    catch
                        keyboard
                    end

                    bb=map_coords_proxy(bb,structure);
                    gv=cell(3,1);
                    for dim=1:3
                        gv{dim}=linspace(bb(1,dim),bb(2,dim),size(structure.img,dim));
                    end

                    if atlases.types(atlas)==4 % mixed atlas, divide
                        if side==1
                            bstructure=structure;
                            bXYZ=XYZ;
                            if ~any(gv{1}>0)
                                ea_error('Mixed atlas does not show positive voxels on the right side');
                            end
                            structure.img=structure.img(gv{1}>0,:,:);
                            gv{1}=gv{1}(gv{1}>0);
                            try
                                XYZ.vx=XYZ.vx(XYZ.mm(:,1)>0,:,:);
                            catch
                                ea_error(['Atlas is empty: ',atlases.names{atlas},'.']);
                            end
                            XYZ.val=XYZ.val(XYZ.mm(:,1)>0,:,:);
                            XYZ.mm=XYZ.mm(XYZ.mm(:,1)>0,:,:);

                            structure.dim=[length(gv{1}),length(gv{2}),length(gv{3})];
                        elseif side==2
                            if ~any(gv{1}<0)
                                ea_error('Mixed atlas does not show positive voxels on the left side');
                            end
                            structure.img=structure.img(gv{1}<0,:,:);

                            gv{1}=gv{1}(gv{1}<0);
                            XYZ.vx=XYZ.vx(XYZ.mm(:,1)<0,:,:);
                            XYZ.val=XYZ.val(XYZ.mm(:,1)<0,:,:);
                            XYZ.mm=XYZ.mm(XYZ.mm(:,1)<0,:,:);

                            structure.dim=[length(gv{1}),length(gv{2}),length(gv{3})];
                        end
                    end

                    [X,Y,Z]=meshgrid(gv{1},gv{2},gv{3});

                    % Xvx=linspace(1,length(gv{1}),3*length(gv{1}));
                    % Yvx=linspace(1,length(gv{2}),3*length(gv{2}));
                    % Zvx=linspace(1,length(gv{3}),3*length(gv{3}));
                    % [Xq,Yq,Zq]=meshgrid(interp1(gv{1},Xvx),...
                    %     interp1(gv{2},Yvx),...
                    %     interp1(gv{3},Zvx));

                    thresh=ea_detthresh(atlases,atlas,structure.img);
                    ea_addnii2lf(atlases,atlas,thresh,options,root,mifix)
                    try
                        fv=isosurface(X,Y,Z,permute(structure.img,[2,1,3]),thresh);
                        %fv=isosurface(Xq,Yq,Zq,permute(interp3(nii.img,Xvx,Yvx,Zvx),[2,1,3]),thresh);
                    catch
                        keyboard
                    end

                    fvc=isocaps(X,Y,Z,permute(structure.img,[2,1,3]),thresh);
                    fv.faces=[fv.faces;fvc.faces+size(fv.vertices,1)];
                    fv.vertices=[fv.vertices;fvc.vertices];
                    if options.prefs.hullsmooth
                        try
                            fv=ea_smoothpatch(fv,[],ceil(options.prefs.hullsmooth/2));
                        catch
                            keyboard
                        end
                    end
                    
                    if ischar(options.prefs.hullsimplify)   % for 'auto' hullsimplify
                        % get to 700 faces
                        simplify=700/length(fv.faces);
                        if simplify < 1 % skip volumes with fewer than 700 faces
                            fv=reducepatch(fv,simplify);
                        end
                    else
                        if options.prefs.hullsimplify<1 && options.prefs.hullsimplify>0
                            fv=reducepatch(fv,options.prefs.hullsimplify);
                        elseif options.prefs.hullsimplify>1
                            simplify=options.prefs.hullsimplify/length(fv.faces);
                            fv=reducepatch(fv,simplify);
                        end
                    end
                    
                    
                    try % works > ML 2015:
                        tr=triangulation(fv.faces,fv.vertices);
                        normals{atlas,side} = vertexNormal(tr);
                    catch % workaround for older versions:
                        % temporally plot atlas to get vertex normals..
                        tmp=figure('visible','off');
                        tp=patch(fv,'VertexNormalsMode','manual');
                        set(tp,'VertexNormalsMode','manual')
                        normals{atlas,side}=get(tp,'VertexNormals')*-1;
                        delete(tmp);
                    end

            

                    iroi{atlas,side}=roi; % later stored
                    icdat{atlas,side}=cdat; % later stored
                    try
                        iXYZ{atlas,side}=XYZ; % later stored
                    catch
                        keyboard
                    end

                    ipixdim{atlas,side}=structure.voxsize(1:3); % later stored
                    icolorc{atlas,side}=colorc; % later stored
                    pixdim=ipixdim{atlas,side};
                    atlascnt=atlascnt+1;

                elseif isfield(structure, 'fibers') % fibertract
                    % concat fibers to one patch object
                    addobjr=ea_showfiber(structure.fibers,structure.idx,colorc);

                    fv.vertices=addobjr.Vertices;
                    fv.faces=addobjr.Faces;
                    delete(addobjr);

                    structure.mm=structure.fibers;
                    structure=rmfield(structure,'fibers');
                    iXYZ{atlas,side} = structure;

                    ifv{atlas,side}=fv;
                    icdat{atlas,side}=[];
                    ipixdim{atlas,side}='fibers';
                    icolorc{atlas,side}=colorc;
                    normals{atlas,side}=[];
                    try
                        atlases.colors(atlas); % check if predefined color exists
                    catch
                        atlases.colors(atlas)=atlas*(maxcolor/length(atlases.names));
                    end
                elseif isfield(structure, 'isdiscfibers') % discriminative fibers
                    ifv{atlas,side} = [];
                    icdat{atlas,side} = [];
                    iXYZ{atlas,side} = [];
                    ipixdim{atlas,side} = 'discfibers';
                    icolorc{atlas,side} = [];
                    normals{atlas,side} = [];
                    atlases.colors(atlas) = NaN;
                end
            end
        end

        %ea_dispercent(1,'end');

        % finish gm_mask file for leadfield computation.
        try % fibertract only atlases dont have a gm_mask
            V=spm_vol([root,filesep,mifix,options.atlasset,filesep,'gm_mask.nii']);
            X=spm_read_vols(V);

            X(X<1.5)=0;
            X(X>1.5)=1;

            V.dt=[16,0];
            delete([root,filesep,mifix,options.atlasset,filesep,'gm_mask.nii']);
            spm_write_vol(V,X);
            ea_crop_nii([root,filesep,mifix,options.atlasset,filesep,'gm_mask.nii']);

            clear X V
            
            ea_reslice_nii([root,filesep,mifix,options.atlasset,filesep,'gm_mask.nii'],[root,filesep,mifix,options.atlasset,filesep,'gm_mask.nii'],[0.3,0.3,0.3],0,0,1,[],[],3);
            spm_smooth([root,filesep,mifix,options.atlasset,filesep,'gm_mask.nii'],[root,filesep,mifix,options.atlasset,filesep,'gm_mask.nii'],[1,1,1]);
            gzip([root,filesep,mifix,options.atlasset,filesep,'gm_mask.nii']);
            delete([root,filesep,mifix,options.atlasset,filesep,'gm_mask.nii']);
        end

        % save table information that has been generated from nii files (on first run with this atlas set).
        atlases.fv=ifv;
        atlases.cdat=icdat;
        atlases.XYZ=iXYZ;
        atlases.pixdim=ipixdim;
        atlases.colorc=icolorc;
        atlases.normals=normals;

        atlases.rebuild=0; % always reset rebuild flag.
        save([root,filesep,mifix,options.atlasset,filesep,'atlas_index.mat'],'atlases','-v7.3');
    end
end


function structure=load_structure(fname)

if strcmp(fname(end-2:end),'.gz')
    wasgzip=1;
    gunzip(fname);
    fname=fname(1:end-3);
else
    wasgzip=0;
end

if strcmp(fname(end-3:end),'.nii') % volumetric
    warning('off');
    ea_crop_nii(fname);
    structure=ea_load_nii(fname);
    warning('on');

    if ~all(abs(structure.voxsize)<=0.8)
        ea_reslice_nii(fname,fname,[0.4,0.4,0.4],0,0,0,[],[],1);
        structure=ea_load_nii(fname);
    end

    if wasgzip
        gzip(fname);
        delete(fname); % since gunzip makes a copy of the zipped file.
    end
elseif strcmp(fname(end-3:end),'.trk') % tracts in trk format
    [fibers,idx]=ea_loadfibertracts(fname,1);
    structure.fibers=fibers;
    structure.idx=idx;

    if wasgzip
        if strcmp(fname(end-3:end),'.trk') % also delete converted .mat file
            [pth,fn]=fileparts(fname);
            delete(fullfile(pth,[fn,'.mat']));
        end
        delete(fname); % since gunzip makes a copy of the zipped file.
    end
elseif strcmp(fname(end-3:end),'.mat')
    warning('off', 'MATLAB:load:variableNotFound');
    if ~isempty(fieldnames(load(fname, 'ea_fibformat'))) % tracts in trk format
        [fibers,idx]=ea_loadfibertracts(fname,1);
        structure.fibers=fibers;
        structure.idx=idx;
    elseif ~isempty(fieldnames(load(fname, 'vals'))) % discriminative fibers
        structure.isdiscfibers = 1;
        % Set default color (blue and red) if not found in mat.
        if isempty(fieldnames(load(fname, 'fibcolor')))
            fibcolor = [0 0 1;1 0 0];
            save(fname, 'fibcolor', '-append');
        end
    end
    if wasgzip
        delete(fname); % since gunzip makes a copy of the zipped file.
    end
end


function sides=detsides(opt)

switch opt
    case 1 % left hemispheric atlas
        sides=1;
    case 2 % right hemispheric atlas
        sides=2;
    case 3
        sides=1:2;
    case 4
        sides=1:2;
    case 5
        sides=1; % midline
end


function coords=map_coords_proxy(XYZ,V)

XYZ=[XYZ';ones(1,size(XYZ,1))];

coords=V.mat*XYZ;
coords=coords(1:3,:)';


function C=rgb(C) % returns rgb values for the colors.

C = rem(floor((strfind('kbgcrmyw', C) - 1) * [0.25 0.5 1]), 2);


function reb=checkrebuild(atlases,options,root,mifix)

reb=1;

if isfield(atlases,'roi')
    reb=0;
end
% if ~exist([root,filesep,mifix,options.atlasset,filesep,'gm_mask.nii'],'file');
%     reb=1;
% end
try
    if atlases.rebuild
        reb=1;
    end
end


function ea_addnii2lf(atlases,atlas,thresh,options,root,mifix)

switch atlases.types(atlas)
    case 1 % right hemispheric atlas.
        atlnames{1}=[root,filesep,mifix,options.atlasset,filesep,'rh',filesep,atlases.names{atlas}];
    case 2 % left hemispheric atlas.
        atlnames{1}=[root,filesep,mifix,options.atlasset,filesep,'lh',filesep,atlases.names{atlas}];
    case 3 % both-sides atlas composed of 2 files.
        atlnames{1}=[root,filesep,mifix,options.atlasset,filesep,'lh',filesep,atlases.names{atlas}];
        atlnames{2}=[root,filesep,mifix,options.atlasset,filesep,'rh',filesep,atlases.names{atlas}];
    case 4 % mixed atlas (one file with both sides information).
        atlnames{1}=[root,filesep,mifix,options.atlasset,filesep,'mixed',filesep,atlases.names{atlas}];
    case 5 % midline atlas (one file with both sides information.
        atlnames{1}=[root,filesep,mifix,options.atlasset,filesep,'midline',filesep,atlases.names{atlas}];
end
[options] = ea_assignpretra(options);
for atl=1:length(atlnames)
    atlname=atlnames{atl};

    if strcmp(atlname(end-2:end),'.gz')
        wasgzip=1;
        gunzip(atlname);
        delete(atlname);
        atlname=atlname(1:end-3);
    else
        wasgzip=0;
    end

    if ~exist([root,filesep,mifix,options.atlasset,filesep,'gm_mask.nii'],'file') % first atlas, generate empty hdtemplate in atlas dir...
        if (~exist([root,filesep,mifix,options.atlasset,filesep,'gm_mask.nii.gz'],'file'))
            if ~options.native
                load([ea_space,'ea_space_def.mat'])
                copyfile([ea_space,spacedef.templates{1},'.nii'],[root,filesep,mifix,options.atlasset,filesep,'gm_mask.nii']);
                ea_reslice_nii([root,filesep,mifix,options.atlasset,filesep,'gm_mask.nii'],[root,filesep,mifix,options.atlasset,filesep,'gm_mask.nii'],...
                    [0.3,0.3,0.3],0,0,1,[],[],1);
            else
                copyfile([options.root,options.patientname,filesep,options.prefs.prenii_unnormalized],[root,filesep,mifix,options.atlasset,filesep,'gm_mask.nii']);
            end
            V=spm_vol([root,filesep,mifix,options.atlasset,filesep,'gm_mask.nii']);
            X=spm_read_vols(V);
            V.dt=[16,0];
            V.pinfo=[1;0;352];
            X(:)=0;
            spm_write_vol(V,X);
            clear V X
        else
            gunzip([root,filesep,mifix,options.atlasset,filesep,'gm_mask.nii.gz']);
            delete([root,filesep,mifix,options.atlasset,filesep,'gm_mask.nii.gz']);
        end
    end

    % add atlas file to hdtemplate in atlas dir
    matlabbatch{1}.spm.util.imcalc.input = {
        [root,filesep,mifix,options.atlasset,filesep,'gm_mask.nii,1'];
        [atlname,',1']
        };
    matlabbatch{1}.spm.util.imcalc.output = 'gm_mask.nii';
    matlabbatch{1}.spm.util.imcalc.outdir = {[root,filesep,mifix,options.atlasset,filesep]};
    matlabbatch{1}.spm.util.imcalc.expression = ['i1+(i2>',num2str(thresh),')'];
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 16;
    if atlases.tissuetypes(atlas)==1
        spm_jobman('run',{matlabbatch});
    end
    clear matlabbatch

    if wasgzip
        gzip(atlname); 
        delete(atlname); % since gunzip makes a copy of the zipped file.
    end
end

%figure, patch(afv)
