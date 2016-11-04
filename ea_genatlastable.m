function atlases=ea_genatlastable(varargin)
% This function reads in atlases in the Lead-dbs/atlases directory and
% generates a table of all available atlas files.
% Atlastypes:   1 ? LH
%               2 ??RH
%               3 ? both hemispheres (2 files present both in lhs and rhs
%               folder
%               4 ? mixed (one file with one cluster on each hemisphere)
%               5 ? midline (one file with one cluster in total)
%
%
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn


atlases=varargin{1};
root=varargin{2};
options=varargin{3};
if nargin==4
    mifix=varargin{4};
else
    mifix='';
end

if isempty(atlases) % create from scratch - if not empty, rebuild flag has been set.
    disp('Generating Atlas table (first run with new atlas only). This may take a while...');
    lhcell=cell(0); rhcell=cell(0); mixedcell=cell(0); midlinecell=cell(0);
    delete([root,'atlases',filesep,mifix,options.atlasset,filesep,'lh',filesep,'*_temp.ni*']);
    lhatlases=dir([root,'atlases',filesep,mifix,options.atlasset,filesep,'lh',filesep,'*.ni*']);
    lhtrks=dir([root,'atlases',filesep,mifix,options.atlasset,filesep,'lh',filesep,'*.tr*']);
    lhmats=dir([root,'atlases',filesep,mifix,options.atlasset,filesep,'lh',filesep,'*.mat']);
    lhatlases=[lhatlases;lhtrks;lhmats];


    for i=1:length(lhatlases);
        lhcell{i}=lhatlases(i).name;
    end
    delete([root,'atlases',filesep,mifix,options.atlasset,filesep,'rh',filesep,'*_temp.ni*']);
    rhatlases=dir([root,'atlases',filesep,mifix,options.atlasset,filesep,'rh',filesep,'*.ni*']);
    rhtrks=dir([root,'atlases',filesep,mifix,options.atlasset,filesep,'rh',filesep,'*.tr*']);
    rhmats=dir([root,'atlases',filesep,mifix,options.atlasset,filesep,'rh',filesep,'*.mat']);
    rhatlases=[rhatlases;rhtrks;rhmats];

    for i=1:length(rhatlases);
        rhcell{i}=rhatlases(i).name;
    end
    delete([root,'atlases',filesep,mifix,options.atlasset,filesep,'mixed',filesep,'*_temp.ni*']);
    mixedatlases=dir([root,'atlases',filesep,mifix,options.atlasset,filesep,'mixed',filesep,'*.ni*']);
    mixedtrks=dir([root,'atlases',filesep,mifix,options.atlasset,filesep,'mixed',filesep,'*.tr*']);
    mixedmats=dir([root,'atlases',filesep,mifix,options.atlasset,filesep,'mixed',filesep,'*.mat']);
    mixedatlases=[mixedatlases;mixedtrks;mixedmats];

    for i=1:length(mixedatlases);
        mixedcell{i}=mixedatlases(i).name;
    end
    delete([root,'atlases',filesep,mifix,options.atlasset,filesep,'midline',filesep,'*_temp.ni*']);

    midlineatlases=dir([root,'atlases',filesep,mifix,options.atlasset,filesep,'midline',filesep,'*.ni*']);
    midlinetrks=dir([root,'atlases',filesep,mifix,options.atlasset,filesep,'midline',filesep,'*.tr*']);
    midlinemats=dir([root,'atlases',filesep,mifix,options.atlasset,filesep,'midline',filesep,'*.mat']);
    midlineatlases=[midlineatlases;midlinetrks;midlinemats];

    for i=1:length(midlineatlases);
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
    typecell=[repmat(1,1,length(rhcell)),repmat(2,1,length(lhcell)),repmat(3,1,length(bothcell)),repmat(4,1,length(mixedcell)),repmat(5,1,length(midlinecell))];
    atlases.names=allcell;
    atlases.types=typecell;
    atlases.rebuild=0;
    atlases.threshold.type='relative_intensity';
    atlases.threshold.value=0.5;

end

mcr=ea_checkmacaque(options);

if checkrebuild(atlases,options,root,mifix)


    %% build iXYZ tables:

    maxcolor=64; % change to 45 to avoid red / 64 to use all colors


    nm=[0:2]; % native and mni
    try
        nmind=[options.atl.pt,options.atl.can,options.atl.ptnative]; % which shall be performed?
    catch
        nmind=[0 1 0];
    end
    nm=nm(logical(nmind)); % select which shall be performed.



    for nativemni=nm % switch between native and mni space atlases.

        switch nativemni
            case 0
                root=[options.root,options.patientname,filesep];
            case 1
                root=[options.earoot,mcr];
            case 2
                root=[options.root,options.patientname,filesep];
        end

        atlascnt=1;


        % iterate through atlases, visualize them and write out stats.
        disp('Building atlas table...');
        for atlas=1:length(atlases.names)
            %ea_dispercent(atlas/length(atlases.names));
            switch atlases.types(atlas)
                case 1 % right hemispheric atlas.
                    nii=load_nii_crop([root,'atlases',filesep,mifix,options.atlasset,filesep,'rh',filesep,atlases.names{atlas}],options);
                case 2 % left hemispheric atlas.
                    nii=load_nii_crop([root,'atlases',filesep,mifix,options.atlasset,filesep,'lh',filesep,atlases.names{atlas}],options);
                case 3 % both-sides atlas composed of 2 files.
                    lnii=load_nii_crop([root,'atlases',filesep,mifix,options.atlasset,filesep,'lh',filesep,atlases.names{atlas}],options);
                    rnii=load_nii_crop([root,'atlases',filesep,mifix,options.atlasset,filesep,'rh',filesep,atlases.names{atlas}],options);
                case 4 % mixed atlas (one file with both sides information).
                    nii=load_nii_crop([root,'atlases',filesep,mifix,options.atlasset,filesep,'mixed',filesep,atlases.names{atlas}],options);
                case 5 % midline atlas (one file with both sides information.
                    nii=load_nii_crop([root,'atlases',filesep,mifix,options.atlasset,filesep,'midline',filesep,atlases.names{atlas}],options);
            end



            for side=detsides(atlases.types(atlas));
                if atlases.types(atlas)==3 % both-sides atlas composed of 2 files.
                    if side==1
                        nii=rnii;
                    elseif side==2
                        nii=lnii;
                    end
                end






                colornames='bgcmywkbgcmywkbgcmywkbgcmywkbgcmywkbgcmywkbgcmywkbgcmywkbgcmywk'; % red is reserved for the VAT.

                colorc=colornames(1);
                colorc=rgb(colorc);
                if isfield(nii,'img') % volumetric atlas
                    
                    if options.prefs.hullsmooth
                        nii.img = smooth3(nii.img,'gaussian',options.prefs.hullsmooth);
                  
                    end
                    
                    [xx,yy,zz]=ind2sub(size(nii.img),find(nii.img>0)); % find 3D-points that have correct value.
                    vv=nii.img(nii.img(:)>0);

                    if ~isempty(xx)

                        XYZ.vx=[xx,yy,zz]; % concatenate points to one matrix.
                        XYZ.val=vv;
                        XYZ.mm=map_coords_proxy(XYZ.vx,nii); % map to mm-space
                        XYZ.dims=nii.hdr.dime.pixdim;


                    else
                        XYZ.vx=[];
                        XYZ.val=[];
                        XYZ.mm=[];
                        XYZ.dims=nii.hdr.dime.pixdim;

                    end

                    %surface(xx(1:10)',yy(1:10)',zz(1:10)',ones(10,1)');
                    %             hold on








                    if atlases.types(atlas)==4 && side==2 % restore from backup
                        nii=bnii;
                        XYZ.mm=bXYZ.mm;
                        XYZ.val=bXYZ.val;
                        XYZ.vx=bXYZ.vx;
                    end
                    try
                        bb=[0,0,0;size(nii.img)];
                    catch
                        keyboard
                    end
                    bb=map_coords_proxy(bb,nii);
                    gv=cell(3,1);
                    for dim=1:3
                        gv{dim}=linspace(bb(1,dim),bb(2,dim),size(nii.img,dim));
                    end



                    if atlases.types(atlas)==4 % mixed atlas, divide
                        if side==1
                            bnii=nii;
                            bXYZ=XYZ;
                            if ~any(gv{1}>0)
                                ea_error('Mixed atlas does not show positive voxels on the right side');
                            end
                            nii.img=nii.img(gv{1}>0,:,:);
                            gv{1}=gv{1}(gv{1}>0);

                            XYZ.vx=XYZ.vx(XYZ.mm(:,1)>0,:,:);
                            XYZ.val=XYZ.val(XYZ.mm(:,1)>0,:,:);
                            XYZ.mm=XYZ.mm(XYZ.mm(:,1)>0,:,:);


                            nii.dim=[length(gv{1}),length(gv{2}),length(gv{3})];
                        elseif side==2
                            if ~any(gv{1}<0)
                                ea_error('Mixed atlas does not show positive voxels on the left side');
                            end
                            nii.img=nii.img(gv{1}<0,:,:);
                            
                            gv{1}=gv{1}(gv{1}<0);
                            XYZ.vx=XYZ.vx(XYZ.mm(:,1)<0,:,:);
                            XYZ.val=XYZ.val(XYZ.mm(:,1)<0,:,:);
                            XYZ.mm=XYZ.mm(XYZ.mm(:,1)<0,:,:);

                            nii.dim=[length(gv{1}),length(gv{2}),length(gv{3})];
                        end
                    end


                    [X,Y,Z]=meshgrid(gv{1},gv{2},gv{3});
                    
                    
%                     Xvx=linspace(1,length(gv{1}),3*length(gv{1}));
%                     Yvx=linspace(1,length(gv{2}),3*length(gv{2}));
%                     Zvx=linspace(1,length(gv{3}),3*length(gv{3}));
%                     [Xq,Yq,Zq]=meshgrid(interp1(gv{1},Xvx),...
%                         interp1(gv{2},Yvx),...
%                         interp1(gv{3},Zvx));

  
                    thresh=ea_detthresh(atlases,atlas,nii.img);
                    ea_addnii2lf(atlases,atlas,thresh,options,root,mifix)
try
                    fv=isosurface(X,Y,Z,permute(nii.img,[2,1,3]),thresh);
                    %fv=isosurface(Xq,Yq,Zq,permute(interp3(nii.img,Xvx,Yvx,Zvx),[2,1,3]),thresh);
                    
catch
    keyboard
end
                    fvc=isocaps(X,Y,Z,permute(nii.img,[2,1,3]),thresh);
                    fv.faces=[fv.faces;fvc.faces+size(fv.vertices,1)];
                    fv.vertices=[fv.vertices;fvc.vertices];




                    try % works only in ML 2015:
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

                    % set cdata

                try % check if explicit color info for this atlas is available.
                    cdat=abs(repmat(atlases.colors(atlas),length(fv.vertices),1));
                catch
                    cdat=abs(repmat(atlas*(maxcolor/length(atlases.names)),length(fv.vertices),1));
                    atlases.colors(atlas)=atlas*(maxcolor/length(atlases.names));
                end

                ifv{atlas,side}=fv; % later stored
                icdat{atlas,side}=cdat; % later stored
                try
                    iXYZ{atlas,side}=XYZ; % later stored
                catch
                    keyboard
                end
                ipixdim{atlas,side}=nii.hdr.dime.pixdim(1:3); % later stored

                icolorc{atlas,side}=colorc; % later stored

                pixdim=ipixdim{atlas,side};




                atlascnt=atlascnt+1;

                elseif isfield(nii,'fibers') % fibertract
                  %% concat fibers to one patch object
                 [~,alnm]=fileparts(atlases.names{atlas});

                    ea_dispercent(0,['Concatenating ',alnm]);
                fibmax=length(nii.idx);
                fcnt=1;

                for fib=1:fibmax
                    ea_dispercent(fib/fibmax);
                    thisfib=nii.fibers(nii.fibers(:,4)==fib,:);
                    if size(thisfib,1)>5 % neglect very small fibertracts
                        % set 4th dim color



                        %atlassurfs(atlascnt,fib)=ea_plot3(thisfib(:,1),thisfib(:,2),thisfib(:,3),'-','color',[rr,gg,bb]);
                        %atlassurfs(atlascnt,fcnt)=ea_plot3t(thisfib(:,1),thisfib(:,2),thisfib(:,3),0.1,[rr,gg,bb],6,0);
                        
                            [~,thisfv]=ea_plot3t(thisfib(:,1),thisfib(:,2),thisfib(:,3),0.1,'r',6,0);
                            % need to manually shuffle the results in case
                            % of mixture fibertracking/nuclei atlases.
                            fv(fcnt).vertices=thisfv.vertices;
                            fv(fcnt).faces=thisfv.faces;
                            fv(fcnt).facevertexdata=thisfv.facevertexcdata;

                        fcnt=fcnt+1;
                    end
                end
                                    ea_dispercent(1,'end');


                fv=ea_concatfv(fv);
                if length(fv.vertices)>200000
                   simpl=200000/length(fv.vertices);
                   fv=reducepatch(fv,simpl);
                end


                    nii.mm=nii.fibers;
                    nii=rmfield(nii,'fibers');
                    iXYZ{atlas,side}=nii;
                    icolorc{atlas,side}=colorc;
                    ifv{atlas,side}=fv;
                    ipixdim{atlas,side}='fibers';
                                        icdat{atlas,side}=[];
                                        normals{atlas,side}=[];
                    try
                        atlases.colors(atlas); % check if predefined color exists
                    catch
                        atlases.colors(atlas)=atlas*(maxcolor/length(atlases.names));
                    end
                end

            end
        end


        %ea_dispercent(1,'end');

        % finish gm_mask file for leadfield computation.
        try % fibertract only atlases dont have a gm_mask
        V=spm_vol([root,'atlases',filesep,mifix,options.atlasset,filesep,'gm_mask.nii']);
        X=spm_read_vols(V);

        X(X<1.5)=0;
        X(X>1.5)=1;

        V.dt=[4,0];

        spm_write_vol(V,X);

        clear X V
        ea_crop_nii([root,'atlases',filesep,mifix,options.atlasset,filesep,'gm_mask.nii']);
        end

        % save table information that has been generated from nii files (on first run with this atlas set).

        atlases.fv=ifv;
        atlases.cdat=icdat;
        atlases.XYZ=iXYZ;
        atlases.pixdim=ipixdim;
        atlases.colorc=icolorc;

        atlases.normals=normals;


        atlases.rebuild=0; % always reset rebuild flag.
        save([root,'atlases',filesep,mifix,options.atlasset,filesep,'atlas_index.mat'],'atlases','-v7.3');


    end

end



function nii=load_nii_crop(fname,options)

if strcmp(fname(end-2:end),'.gz')
    wasgzip=1;
    gunzip(fname);
    delete(fname);
    fname=fname(1:end-3);
else
    wasgzip=0;
end
%try

if strcmp(fname(end-3:end),'.nii') % volumetric

    ea_crop_nii(fname);
    nii=spm_vol(fname);

    nii.img=spm_read_vols(nii);

    %catch

    %end


    nii.hdr.dime.pixdim=nii.mat(logical(eye(4)));
    if ~all(abs(nii.hdr.dime.pixdim(1:3))<=0.7)
        ea_reslice_nii(fname,fname,[0.4,0.4,0.4],3);

        nii=spm_vol(fname);
        nii.img=spm_read_vols(nii);
        nii.hdr.dime.pixdim=nii.mat(logical(eye(4)));

    end
    if wasgzip
        gzip(fname); % since gunzip makes a copy of the zipped file.
        delete(fname);
    end

elseif strcmp(fname(end-3:end),'.trk') || strcmp(fname(end-3:end),'.mat') % tracts in mat format % tracts in trk format

    [fibers,idx]=ea_loadfibertracts(fname);
    nii.fibers=fibers;
    nii.idx=idx;


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

if isfield(atlases,'fv')
    reb=0;
    if ~isfield(atlases.XYZ{1,1},'mm')
        reb=1;
    end
end
% if ~exist([root,'atlases',filesep,mifix,options.atlasset,filesep,'gm_mask.nii'],'file');
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
        atlnames{1}=[root,'atlases',filesep,mifix,options.atlasset,filesep,'rh',filesep,atlases.names{atlas}];
    case 2 % left hemispheric atlas.
        atlnames{1}=[root,'atlases',filesep,mifix,options.atlasset,filesep,'lh',filesep,atlases.names{atlas}];
    case 3 % both-sides atlas composed of 2 files.
        atlnames{1}=[root,'atlases',filesep,mifix,options.atlasset,filesep,'lh',filesep,atlases.names{atlas}];
        atlnames{2}=[root,'atlases',filesep,mifix,options.atlasset,filesep,'rh',filesep,atlases.names{atlas}];
    case 4 % mixed atlas (one file with both sides information).
        atlnames{1}=[root,'atlases',filesep,mifix,options.atlasset,filesep,'mixed',filesep,atlases.names{atlas}];
    case 5 % midline atlas (one file with both sides information.
        atlnames{1}=[root,'atlases',filesep,mifix,options.atlasset,filesep,'midline',filesep,atlases.names{atlas}];
end

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

    if ~exist([root,'atlases',filesep,mifix,options.atlasset,filesep,'gm_mask.nii'],'file') % first atlas, generate empty hdtemplate in atlas dir...
        copyfile([options.earoot,'templates',filesep,'mni_hires.nii'],[root,'atlases',filesep,mifix,options.atlasset,filesep,'gm_mask.nii']);
        V=spm_vol([root,'atlases',filesep,mifix,options.atlasset,filesep,'gm_mask.nii']);
        X=spm_read_vols(V);
        X(:)=0;
        spm_write_vol(V,X);
        clear V X
    end


    % add atlas file to hdtemplate in atlas dir
    matlabbatch{1}.spm.util.imcalc.input = {
        [root,'atlases',filesep,mifix,options.atlasset,filesep,'gm_mask.nii,1'];
        [atlname,',1']
        };
    matlabbatch{1}.spm.util.imcalc.output = 'gm_mask.nii';
    matlabbatch{1}.spm.util.imcalc.outdir = {[root,'atlases',filesep,mifix,options.atlasset,filesep]};
    matlabbatch{1}.spm.util.imcalc.expression = ['i1+(i2>',num2str(thresh),')'];
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    spm_jobman('run',{matlabbatch});
    clear matlabbatch

    if wasgzip
        gzip(atlname); % since gunzip makes a copy of the zipped file.
        delete(atlname);
    end
end



%figure, patch(afv)
