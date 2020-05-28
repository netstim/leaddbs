function [PL]=ea_showfibers_volume(resultfig,options,hmchanged)
% This function shows fiber-connectivity from a volume defined by a nx3
% point-list (volume). If stimparams.showconnectivities is set, connected
% areas to the volume are also visualized. To do so, the function uses
% inhull.m which is covered by the BSD-license (see below).
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

PL.ht=uitoolbar(resultfig);
set(0,'CurrentFigure',resultfig)
colornames='rbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywk';

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
        eval([togglenames{but},'=repmat(1,expand,length(options.sides));']);
    end

    setappdata(resultfig,togglenames{but},eval(togglenames{but}));

end
clear expand

load([options.root,options.patientname,filesep,'ea_stats']);

try
    upriorvatlength=length(ea_stats.vat)+1;
    upriorftlength=length(ea_stats.ft)+1;
catch
    upriorvatlength=1;
    upriorftlength=1;
end

% assign the place where to write stim stats into struct

if isfield(options,'groupmode')
    if options.groupmode
        S.label=['gs_',options.groupid]; % hard reset label to group stim label
    end
end

[ea_stats,thisstim]=ea_assignstimcnt(ea_stats,S);

if isstruct(VAT{1}.VAT) || isstruct(VAT{2}.VAT) % e.g. simbio model used
    vat=1;
    for side=1:length(options.sides)
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

for side=1:length(options.sides)
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

                vatfv.faces=K(side).K{vat}; vatfv.vertices=VAT{side}.VAT{vat};
                vatfv=reducepatch(vatfv,0.05);
                Vcent=mean(vatfv.vertices);

                atlasName = options.atlasset;
                load([ea_space(options,'atlases'),atlasName,filesep,'atlas_index.mat']);
                for atlas=1:size(atlases.XYZ,1)
                    if stimparams(side).volume(vat)>0 % stimulation on in this VAT,
                        clear thisatl

                        try % for midline or combined atlases, only the right side atlas is used.
                            if isempty(atlases.XYZ{atlas,side}) % for midline or combined atlases, only the right side atlas is used.
                                thisatl=atlases.XYZ{atlas,1}.mm;
                                tpd=atlases.pixdim{atlas,1};
                            else
                                thisatl=atlases.XYZ{atlas,side}.mm;
                                tpd=atlases.pixdim{atlas,side};
                            end
                        catch
                            thisatl=atlases.XYZ{atlas,1}.mm;
                            tpd=atlases.pixdim{atlas,1};
                        end

                        tpv=abs(tpd(1))*abs(tpd(2))*abs(tpd(3)); % volume of one voxel in mm^3.

                        ea_stats.stimulation(thisstim).vat(side,vat).AtlasIntersection(atlas)=sum(ea_intriangulation(vatfv.vertices,vatfv.faces,thisatl))*tpv;
                        ea_stats.stimulation(thisstim).vat(side,vat).nAtlasIntersection(atlas)=ea_stats.stimulation(thisstim).vat(side,vat).AtlasIntersection(atlas)/stimparams(1,side).volume(vat);
                        ea_stats.stimulation(thisstim).vat(side,vat).volume=stimparams(1,side).volume(vat);

                        % now also add efield overlap:
                        if exist(ea_niigz([options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options)...
                                S.label,filesep,'vat_efield_',sidec]),'file')
                            Vefield=ea_load_nii(ea_niigz([options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options)...
                                S.label,filesep,'vat_efield_',sidec]));
                            atlasvoxels=Vefield.mat\[thisatl,ones(length(thisatl),1)]';
                            ea_stats.stimulation(thisstim).efield(side,vat).AtlasIntersection(atlas)=...
                                mean(spm_sample_vol(Vefield,atlasvoxels(1,:),atlasvoxels(2,:),atlasvoxels(3,:),1));
                            ea_stats.stimulation(thisstim).efield(side,vat).nAtlasIntersection(atlas)=...
                                mean(spm_sample_vol(Vefield,atlasvoxels(1,:),atlasvoxels(2,:),atlasvoxels(3,:),1))./...
                                sum(Vefield.img(:));
                            ea_stats.stimulation(thisstim).efield(side,vat).volume=sum(Vefield.img(:));
                        end
                    else % no voltage on this vat, simply set vi to zero.
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
for side=1:length(options.sides)
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
