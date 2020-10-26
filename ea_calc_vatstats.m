function [PL]=ea_calc_vatstats(resultfig,options,hmchanged)
% Calculate VAT stats (atlas intersection and volume)
% Stimulation information w.r.t the VTA will also be stored in stats.
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

PL.ht=uitoolbar(resultfig);
set(0,'CurrentFigure',resultfig)

ea_dispt('Visualizing VTA...');

hold on
% get app data

stimparams=getappdata(resultfig,'stimparams');
S=getappdata(resultfig,'curS');
for side=1:length(stimparams)
    VAT{side}=stimparams(side).VAT;
end

if ~exist('hmchanged','var')
    hmchanged=1;
end

% clean downstreamfiles if necessary
if hmchanged
    ea_delete(fullfile(options.root,options.patientname,'stimulations',ea_nt(options),S.label,'bihem_vat_left.nii'));
    ea_delete(fullfile(options.root,options.patientname,'stimulations',ea_nt(options),S.label,'bihem_vat_right.nii'));
    ea_delete(fullfile(options.root,options.patientname,'stimulations',ea_nt(options),S.label,'fl_vat_left.nii'));
    ea_delete(fullfile(options.root,options.patientname,'stimulations',ea_nt(options),S.label,'fl_vat_right.nii'));

    ea_delete(fullfile(options.root,options.patientname,'stimulations',ea_nt(options),S.label,'bihem_vat_efield_left.nii'));
    ea_delete(fullfile(options.root,options.patientname,'stimulations',ea_nt(options),S.label,'bihem_vat_efield_right.nii'));
    ea_delete(fullfile(options.root,options.patientname,'stimulations',ea_nt(options),S.label,'fl_vat_efield_left.nii'));
    ea_delete(fullfile(options.root,options.patientname,'stimulations',ea_nt(options),S.label,'fl_vat_efield_right.nii'));

    ea_delete(fullfile(options.root,options.patientname,'stimulations',ea_nt(options),S.label,'bihem_vat_efield_gauss_left.nii'));
    ea_delete(fullfile(options.root,options.patientname,'stimulations',ea_nt(options),S.label,'bihem_vat_efield_gauss_right.nii'));
    ea_delete(fullfile(options.root,options.patientname,'stimulations',ea_nt(options),S.label,'fl_vat_efield_gauss_left.nii'));
    ea_delete(fullfile(options.root,options.patientname,'stimulations',ea_nt(options),S.label,'fl_vat_efield_gauss_right.nii'));

    ea_delete(fullfile(options.root,options.patientname,'stimulations',ea_nt(options),S.label,'vat_seed_compound_dMRI.nii'));
    ea_delete(fullfile(options.root,options.patientname,'stimulations',ea_nt(options),S.label,'vat_seed_compound_dMRI_l.nii'));
    ea_delete(fullfile(options.root,options.patientname,'stimulations',ea_nt(options),S.label,'vat_seed_compound_dMRI_r.nii'));

    ea_delete(fullfile(options.root,options.patientname,'stimulations',ea_nt(options),S.label,'vat_seed_compound_fMRI.nii'));
    ea_delete(fullfile(options.root,options.patientname,'stimulations',ea_nt(options),S.label,'vat_seed_compound_fMRI_l.nii'));
    ea_delete(fullfile(options.root,options.patientname,'stimulations',ea_nt(options),S.label,'vat_seed_compound_fMRI_r.nii'));
end

%prepare statvat exports once if needed.
if options.expstatvat.do % export statvat nifti images.
    tV=spm_vol([ea_space(options),'bb.nii']);
    tnii=spm_read_vols(tV);
    tnii(:)=0;
    %generate mesh of hires MNI
    [x,y,z]=ind2sub(size(tnii),1:numel(tnii));
    templatecoords=[x;y;z;ones(1,length(x))]; clear x y z
    templatecoords=tV.mat*templatecoords;
    templatecoords=templatecoords(1:3,:)';
end

togglenames={'vaton','quivon'};
for but=1:length(togglenames)

    eval([togglenames{but},'=getappdata(resultfig,''',togglenames{but},''');']);
    expand=1;
    if isempty(eval(togglenames{but}))
        %eval([togglenames{but},'=repmat(1,expand,length(options.sides));']);
        %changed to max, as to include for sure the array as large as the maximum side used, 
        %as this code was intended for the bilateral cases
        %maybe will have to change it to minimum have two elements, to always include at least R and L sides
        %For example, before if the side was only Left, the multiplier would have been only 1
        eval([togglenames{but},'=repmat(1,expand,max(options.sides));']);
    end

    setappdata(resultfig,togglenames{but},eval(togglenames{but}));

end
clear expand

load([options.root,options.patientname,filesep,'ea_stats']);

% assign the place where to write stim stats into struct

if isfield(options,'groupmode')
    if options.groupmode
        S.label=['gs_',options.groupid]; % hard reset label to group stim label
    end
end

[ea_stats,thisstim]=ea_assignstimcnt(ea_stats,S);

if (isfield(VAT{1},'VAT') && isstruct(VAT{1}.VAT)) || ((length(VAT)>1) && isfield(VAT{2},'VAT') && isstruct(VAT{2}.VAT)) % e.g. simbio model used
    vat=1;
    for iside=1:length(options.sides)
        side=options.sides(iside);
        try
            nVAT{side}.VAT{vat}=VAT{side}.VAT.vertices;
            K(side).K{vat}=VAT{side}.VAT.faces;
        catch
            nVAT{side}.VAT{vat}=[];
            K(side).K{vat}=[];
        end
    end
    VAT=nVAT;
end

for iside=1:length(options.sides)
    side=options.sides(iside);
    switch side
        case 1
            sidec='right';
        case 2
            sidec='left';
    end
    for vat=1:length(VAT{side}.VAT)
        if ~exist('K','var') % e.g. maedler model used
            try
                K(side).K{vat}=convhulln(VAT{side}.VAT{vat}+randn(size(VAT{side}.VAT{vat}))*0.000001); % create triangulation.
            catch
                keyboard
                if isnan(VAT{side}.VAT) % empty VTA
                    continue
                end
            end
        else % still maedler model used
            try
                K(side).K{vat}; % not defined
            catch
                K(side).K{vat}=convhulln(VAT{side}.VAT{vat}+randn(size(VAT{side}.VAT{vat}))*0.000001); % create triangulation.
            end
        end

        % show vat
        if ~isempty(K(side).K{vat})

            PL.vatsurfs(side,vat)=trisurf(K(side).K{vat},VAT{side}.VAT{vat}(:,1),VAT{side}.VAT{vat}(:,2),VAT{side}.VAT{vat}(:,3),...
                abs(repmat(60,length(VAT{side}.VAT{vat}),1)...
                +randn(length(VAT{side}.VAT{vat}),1)*2)');

            % the following is some code required for
            % Web/JSON/BrainBrowser-Export.
            PL.vatfv(side,vat).vertices=[VAT{side}.VAT{vat}(:,1),VAT{side}.VAT{vat}(:,2),VAT{side}.VAT{vat}(:,3)];
            PL.vatfv(side,vat).faces=K(side).K{vat};
            PL.vatfv(side,vat).normals=get(PL.vatsurfs(side,vat),'Vertexnormals');
            PL.vatfv(side,vat).colors=repmat([1,0,0,0.7],size(PL.vatfv(side,vat).vertices,1),1);

            ea_spec_atlas(PL.vatsurfs(side,vat),'vat',jet,1);

            vatgrad=getappdata(resultfig,'vatgrad');
            if ~isempty(vatgrad)
                try % only one hemisphere could be defined.
                    if stimparams(side).volume
                        reduc=ceil(length(vatgrad(side).x)/100000);

                        PL.quiv(side)=quiver3(vatgrad(side).x(1:reduc:end),vatgrad(side).y(1:reduc:end),vatgrad(side).z(1:reduc:end),vatgrad(side).qx(1:reduc:end),vatgrad(side).qy(1:reduc:end),vatgrad(side).qz(1:reduc:end),0,'w-','LineWidth',1);
                    end
                end
            end

            if options.writeoutstats
                ea_dispt('Writing out stats...');
                load([options.root,options.patientname,filesep,'ea_stats']);
                ea_stats.stimulation(thisstim).vat(side,vat).amp=S.amplitude{side};
                ea_stats.stimulation(thisstim).vat(side,vat).label=S.label;
                ea_stats.stimulation(thisstim).vat(side,vat).contact=vat;
                ea_stats.stimulation(thisstim).vat(side,vat).side=side;
                ea_stats.stimulation(thisstim).label=S.label;

                atlasName = options.atlasset;
                load([ea_space(options,'atlases'),atlasName,filesep,'atlas_index.mat']);

                for atlas=1:length(atlases.names)
                    if any(S.amplitude{side}) % stimulation on
                        switch atlases.types(atlas)
                            case 1 % right hemispheric atlas.
                                atlasfile = [ea_space([],'atlases'),options.atlasset,filesep,'rh',filesep,atlases.names{atlas}];
                            case 2 % left hemispheric atlas.
                                atlasfile = [ea_space([],'atlases'),options.atlasset,filesep,'lh',filesep,atlases.names{atlas}];
                            case 3 % both-sides atlas composed of 2 files.
                                switch sidec
                                    case 'right'
                                        atlasfile = [ea_space([],'atlases'),options.atlasset,filesep,'rh',filesep,atlases.names{atlas}];
                                    case 'left'
                                        atlasfile = [ea_space([],'atlases'),options.atlasset,filesep,'lh',filesep,atlases.names{atlas}];
                                end
                            case 4 % mixed atlas (one file with both sides information).
                                atlasfile = [ea_space([],'atlases'),options.atlasset,filesep,'mixed',filesep,atlases.names{atlas}];
                            case 5 % midline atlas (one file with both sides information.
                                atlasfile = [ea_space([],'atlases'),options.atlasset,filesep,'midline',filesep,atlases.names{atlas}];
                        end

                        if endsWith(atlasfile, {'.nii','.nii.gz'})
                            atlasfile = ea_niigz(atlasfile);
                        else % Skip fiber atlas
                            ea_stats.stimulation(thisstim).vat(side,vat).AtlasIntersection(atlas)=0;
                            ea_stats.stimulation(thisstim).vat(side,vat).nAtlasIntersection(atlas)=0;
                            continue;
                        end

                        vatfile = ea_niigz([options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options),S.label,filesep,'vat_',sidec]);
                        voxsize = prod(ea_detvoxsize(vatfile));

                        ea_stats.stimulation(thisstim).vat(side,vat).AtlasIntersection(atlas)=ea_vta_overlap(vatfile,atlasfile,sidec).*voxsize;
                        ea_stats.stimulation(thisstim).vat(side,vat).nAtlasIntersection(atlas)=ea_stats.stimulation(thisstim).vat(side,vat).AtlasIntersection(atlas)/stimparams(1,side).volume(vat);
                        ea_stats.stimulation(thisstim).vat(side,vat).volume=stimparams(1,side).volume(vat);

                        % now also add efield overlap:
                        if exist(ea_niigz([options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options)...
                                S.label,filesep,'vat_efield_',sidec]),'file')
                            vefieldfile=ea_niigz([options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options),S.label,filesep,'vat_efield_',sidec]);
                            Vefield=load_untouch_nii(vefieldfile);
                            efiedthreshold = options.prefs.machine.vatsettings.horn_ethresh*10^3;
                            Vefield.img(Vefield.img<=efiedthreshold) = 0;
                            ea_stats.stimulation(thisstim).efield(side,vat).AtlasIntersection(atlas)=ea_vta_overlap(vefieldfile,atlasfile,sidec);
                            ea_stats.stimulation(thisstim).efield(side,vat).nAtlasIntersection(atlas)=...
                                ea_stats.stimulation(thisstim).efield(side,vat).AtlasIntersection(atlas)./sum(Vefield.img(:));
                            ea_stats.stimulation(thisstim).efield(side,vat).volume=sum(Vefield.img(:));
                        end
                    else % no stimulation, simply set vi to zero.
                        ea_stats.stimulation(thisstim).vat(side,vat).AtlasIntersection(atlas)=0;
                        ea_stats.stimulation(thisstim).vat(side,vat).nAtlasIntersection(atlas)=0;
                    end
                end

                save([options.root,options.patientname,filesep,'ea_stats'],'ea_stats','-v7.3');
            end
        end
    end

    try
        vatbutton(side)=uitoggletool(PL.ht,'CData',ea_get_icn('vat'),'TooltipString','Volume of activated tissue','OnCallback',{@objvisible,PL.vatsurfs(side,:),resultfig,'vaton',[],side,1},'OffCallback',{@objvisible,PL.vatsurfs(side,:),resultfig,'vaton',[],side,0},'State',getstate(vaton(side)));
        quivbutton(side)=uitoggletool(PL.ht,'CData',ea_get_icn('quiver'),'TooltipString','E-field','OnCallback',{@objvisible,PL.quiv(side),resultfig,'quivon',[],side,1},'OffCallback',{@objvisible,PL.quiv(side),resultfig,'quivon',[],side,0},'State',getstate(quivon(side)));
    end
end

% correct togglestates
for iside=1:length(options.sides)
    side=options.sides(iside);
    if ~vaton(side)
        try
            objvisible([],[],PL.vatsurfs(side,:),resultfig,'vaton',[],side,0)
        end
    end
    if ~quivon(side)
        try
            objvisible([],[],PL.quiv(side),resultfig,'quivon',[],side,0)
        end
    end
end

setappdata(resultfig,'PL',PL);
ea_dispt('');


function objvisible(hobj,ev,atls,resultfig,what,la,side,onoff)
% set visibility
try
    set(atls, 'Visible', getstate(onoff));
catch
    keyboard
end
% log visibility
tvalue=getappdata(resultfig,what);


if isempty(la)
    tvalue(side)=onoff;
else
    tvalue(la,side)=onoff;
end

setappdata(resultfig,what,tvalue);


function str=getstate(val)
switch val
    case 1
        str='on';
    case 0
        str='off';
end
