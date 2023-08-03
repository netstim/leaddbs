function varargout=ea_genvat_horn(varargin)
% This function generates a volume of activated tissue around for each
% electrode based on the methodology described in Horn 2017 AoN.

useSI=1;
vizz=0;
if nargin==5
    % coords=varargin{1}; % Not used anymore, will reload coords using ea_load_reconstruction
    S=varargin{2};
    side=varargin{3};
    options=varargin{4};
    stimname=varargin{5};
elseif nargin==6
    % coords=varargin{1}; % Not used anymore, will reload coords using ea_load_reconstruction
    S=varargin{2};
    side=varargin{3};
    options=varargin{4};
    stimname=varargin{5};
    lgfigure=varargin{6};
elseif nargin==1
    if ischar(varargin{1}) % return name of method.
        varargout{1} = 'SimBio/FieldTrip (see Horn 2017)';
        varargout{2} = true; % Support directed lead
        return
    end
end

thresh=options.prefs.machine.vatsettings.horn_ethresh; %0.2;

if useSI
    thresh=thresh.*(10^3);
end

% S=ea_activecontacts(S);
if ~any(S.activecontacts{side}) % empty VAT, no active contacts.
    ofv.vertices=[0,0,0
        0,0,0
        0,0,0];
    ofv.faces=[1,2,3];
    varargout{1}=ofv;
    varargout{2}=0;
    return
end

%% get electrodes handles // initial parameters:
if exist('lgfigure', 'var')
    resultfig = getappdata(lgfigure,'resultfig');
else
    resultfig = '';
end

% Important to load in reco from a new since we need to decide whether to
% use native or template coordinates. Even when running in template space,
% the native coordinates are sometimes used (VTA is then calculated in native space and ported to template).
[coords_mm,trajectory,markers]=ea_load_reconstruction(options);
elstruct(1).coords_mm=coords_mm;
elstruct(1).trajectory=trajectory;
elstruct(1).name=options.patientname;
elstruct(1).markers=markers;

if ~isempty(resultfig)
    elspec = getappdata(resultfig,'elspec');
    setappdata(resultfig,'elstruct',elstruct);
else
    elspec = options.elspec;
end

options.usediffusion=0; % set to 1 to incorporate diffusion signal (for now only possible using the mesoFT tracker).
coords = coords_mm{side};

% Add stretchfactor to elstruct simply for purpose of checking if headmodel
% changed. Larger stim amplitudes need larger bounding boxes so
% stretchfactor must be incorporated here.
if max(S.amplitude{side})>4
    elstruct.stretchfactor=0.75; %(max(S.amplitude{side})/10);
else
    elstruct.stretchfactor=0.5;
end

stimDir = fullfile(options.subj.stimDir, ea_nt(options), stimname);
ea_mkdir(stimDir);
headmodelDir = fullfile(options.subj.subjDir, 'headmodel', ea_nt(options));
ea_mkdir(headmodelDir);
filePrefix = ['sub-', options.subj.subjId, '_desc-'];

hmchanged=ea_headmodel_changed(options,side,elstruct); % can only use this test once.
assignin('caller','hmchanged',hmchanged);
if hmchanged
    disp('Headmodel needs to be re-calculated. This may take a while...');

    cnt=1;
    mesh.tet=[];
    mesh.pnt=[];
    mesh.tissue=[];
    mesh.tissuelabel={'gray','white','contacts','insulation'};
    % add gm to mesh
    if options.prefs.machine.vatsettings.horn_useatlas
        switch options.prefs.vat.gm
            case 'atlas'
                atlasName = options.prefs.machine.vatsettings.horn_atlasset;
                load([ea_space(options,'atlases'),atlasName,filesep,'atlas_index.mat']);
                if ~isfield(atlases,'tissuetypes')
                    atlases.tissuetypes=ones(length(atlases.names),1);
                end
                for atlas=1:numel(atlases.roi)
                    if isempty(atlases.roi{atlas}.fv) || (atlases.tissuetypes~=1)
                        continue
                    end
                    fv(cnt)=atlases.roi{atlas}.fv;

                    ins=surfinterior(fv(cnt).vertices,fv(cnt).faces);
                    %tissuetype(cnt)=1;
                    cnt=cnt+1;
                end
            case 'tpm'
                c1=ea_load_nii([ea_space(options),'TPM.nii,1']);
                fv=isosurface(c1.img,0.5,'noshare');
                fv.vertices=c1.mat*[fv.vertices,ones(length(fv.vertices),1)]';
                fv.vertices=fv.vertices(1:3,:)';
            case 'mask'
                fv=ea_fem_getmask(options);
                %figure, patch('faces',fv.faces,'vertices',fv.vertices)
        end
    else
        fv=[];
    end

    [elfv,ntissuetype,Y,electrode]=ea_buildelfv(elspec,elstruct,side);
    Ymod=Y;
    success=0;
    for attempt=1:4 % allow four attempts with really small jitters in case scene generates intersecting faces FIX ME this needs a better solution
        try
            [mesh.tet,mesh.pnt,activeidx,wmboundary,centroids,tissuetype]=ea_mesh_electrode(fv,elfv,ntissuetype,electrode,options,S,side,electrode.numel,Ymod,elspec);
            if ~isempty(mesh.tet)
                success=1;
                break
            end
        catch ME
            % The VTA model has led to an intersection of meshes, which
            % can sometimes happen. We will introduce a small jitter to
            % the electrode and try again.
            Ymod=Y+(randn(4)/700); % Very small jitter on transformation which will be used on electrode. - should not exceed ~700. Use vizz below to see effects.
            if attempt>2
                fv=ea_fem_getmask(options,1); % try no surface smoothing for atlas fv - in rare cases this has led to conflicts for tetgen.
            end
            if vizz
                h=figure;
                telfv=elfv;
                for c=1:length(elfv)
                       telfv(c).vertices=Y*[telfv(c).vertices,ones(size(telfv(c).vertices,1),1)]';
                    telfv(c).vertices=telfv(c).vertices(1:3,:)';
                patch(telfv(c),'edgecolor','m','facecolor','none');
                end
                telfv=elfv;
                for c=1:length(elfv)
                    telfv(c).vertices=Ymod*[telfv(c).vertices,ones(size(telfv(c).vertices,1),1)]';
                    telfv(c).vertices=telfv(c).vertices(1:3,:)';
                    patch(telfv(c),'edgecolor','g','facecolor','none');
                end
                axis equal
                view(0,0)
                h.Position=[1000          85         253        1253];
            end
        end
        [~, tetgenName, tetgenExt] = fileparts(ea_getExec(mcpath('tetgen')));
        ea_kill('name', [tetgenName, tetgenExt]);
    end

    if ~success
        if exist('ME', 'var')
            ea_cprintf('CmdWinErrors', '%s\n', ME.message);
        end
        ea_error(['Despite all attempts the VTA model could not be created.\n' ...
            'Please check MATLAB Command Window for detailed error information.\n' ...
            'Ideas: try estimating the VTA model directly in template space and/or without using an atlas to define gray matter.']);
    end

    % replace wmboundary
    tess = mesh.tet(:,1:4);
    tess = sort(tess,2);

    % all faces
    faces=[tess(:,[1 2 3]);tess(:,[1 2 4]); ...
        tess(:,[1 3 4]);tess(:,[2 3 4])];

    % find replicate faces
    faces = sortrows(faces);
    k = find(all(diff(faces)==0,2));

    % delete the internal (shared) faces
    faces([k;k+1],:) = [];

    wmboundary = unique(faces(:))';
    % end replace.

    if vizz
        figure
        hold on
        plot3(mesh.pnt(wmboundary,1),mesh.pnt(wmboundary,2),mesh.pnt(wmboundary,3),'r*');
        plot3(mesh.pnt(:,1),mesh.pnt(:,2),mesh.pnt(:,3),'b.');
    end

    if ~success
        ea_error('Lead-DBS could not solve the current estimation.');
    end

    mesh.tissue=tissuetype;
    meshregions=mesh.tet(:,5);
    mesh.tet=mesh.tet(:,1:4);

    if useSI
        mesh.pnt=mesh.pnt/1000; % in meter
        mesh.unit='m';
    end
	% plot3(mesh.pnt(:,1),mesh.pnt(:,2),mesh.pnt(:,3),'c.');
    %% calculate volume conductor
    ea_dispt('Creating volume conductor...');

    if useSI
        SIfx=1;
    else
        SIfx=1000;
    end

    try
        vol=ea_ft_headmodel_simbio(mesh,'conductivity',SIfx*[options.prefs.machine.vatsettings.horn_cgm options.prefs.machine.vatsettings.horn_cwm 1/(10^(-8)) 1/(10^16)]); % multiply by thousand to use S/mm
        %vol=ea_ft_headmodel_simbio(mesh,'conductivity',SIfx*[0.0915 0.059 1/(10^(-8)) 1/(10^16)]); % multiply by thousand to use S/mm
    catch % reorder elements so not to be degenerated.
        tmesh=mesh;
        tmesh.tet=tmesh.tet(:,[1 2 4 3]);
        vol=ea_ft_headmodel_simbio(tmesh,'conductivity',SIfx*[options.prefs.machine.vatsettings.horn_cgm options.prefs.machine.vatsettings.horn_cwm 1/(10^(-8)) 1/(10^16)]); % multiply by thousand to use S/mm
    end
    if useSI % convert back before saving headmodel
        mesh.pnt=mesh.pnt*1000; % in meter
        mesh.unit='mm';
    end

    save(fullfile(headmodelDir, [filePrefix, 'headmodel', num2str(side),'.mat']), 'vol','mesh','centroids','wmboundary','elfv','meshregions','-v7.3');
    ea_save_hmprotocol(options,side,elstruct,1);
else
    % simply load vol.
    ea_dispt('Loading headmodel...');
    load(fullfile(headmodelDir, [filePrefix, 'headmodel', num2str(side),'.mat']));
    activeidx=ea_getactiveidx(S,side,centroids,mesh,elfv,elspec,meshregions);
end

switch side
    case 1
        sidec='R';
        cnts={'k0','k1','k2','k3','k4','k5','k6','k7'};
    case 2
        sidec='L';
        cnts={'k8','k9','k10','k11','k12','k13','k14','k15'};
end

if ~isfield(S, 'sources')
    S.sources=1:4;
end

for source=S.sources
    stimsource=S.([sidec,'s',num2str(source)]);
    constvol=stimsource.va==1; % constvol is 1 for constant voltage and 0 for constant current.

    for cnt=1:length(cnts)
        if constvol
            U(cnt)=(logical(stimsource.(cnts{cnt}).perc))*stimsource.amp; % do not split amplitude in constant voltage setting.
        else
            U(cnt)=(stimsource.(cnts{cnt}).perc/100)*stimsource.amp;
        end
        if stimsource.(cnts{cnt}).pol==1
            U(cnt)=U(cnt)*-1;
        end
    end

    actInd=find(U); % active contact indices
    if ~isempty(actInd)
        actContact=coords(actInd,:);

        volts=U(U~=0);

        %% calculate voltage distribution based on dipole
        ea_dispt('Calculating voltage distribution...');
        if useSI
            SIfx=1000;
        else
            SIfx=1;
        end
        %ix=knnsearch(vol.pos,dpvx/SIfx); % add dpvx/1000 for m

        if any(volts>0)
            unipolar=0;
            U=U/2;
            %volts=volts/2;
        else
            unipolar=1;
        end
        ix=[];
        voltix=[];

        cnt=1;
        for ac=actInd
            ix=[ix;activeidx(source).con(ac).ix];
            voltix=[voltix;repmat(U(ac),length(activeidx(source).con(ac).ix),1),...
                repmat(cnt,length(activeidx(source).con(ac).ix),1)];
            cnt=cnt+1;
        end

        if isempty(ix)
            rmdir(fullfile(options.subj.subjDir, 'headmodel'),'s'); % the least I can do at this point is to clean up the faulty headmodel.
            ea_error('Something went wrong. Active vertex index not found.');
        end

        if ~constvol
            voltix(:,1)=voltix(:,1)/1000; % from mA to A
            %voltix=voltix;
        end

        potential = ea_apply_dbs(vol,ix,voltix,unipolar,constvol,wmboundary); % output in V. 4 indexes insulating material.
        % save('results','mesh','vol','ix','voltix','unipolar','constvol','wmboundary','potential3v','potential3ma','gradient3v','gradient3ma');

        voltix=voltix(:,1); % get rid of index column
        if vizz
            figure
            hold on
            plot3(mesh.pnt(wmboundary,1),mesh.pnt(wmboundary,2),mesh.pnt(wmboundary,3),'r*');
            plot3(mesh.pnt(:,1),mesh.pnt(:,2),mesh.pnt(:,3),'b.');
        end

        ea_dispt('Calculating E-Field...');
        gradient{source} = ea_calc_gradient(vol,potential); % output in V/m.
        %% stuff by Till to get high EF values for active electrodes
        % this can be adjusted by assigning all tetrahedar belonging to the
        % active electrode a new value:
        % gradient{source}(elec_tet_ix,:) = new_value;
        elec_tet_ix = sub2ind(size(mesh.pnt),vertcat(ix,ix,ix),vertcat(ones(length(ix),1),ones(length(ix),1).*2,ones(length(ix),1).*3));
        elec_tet_ix = find(sum(ismember(mesh.tet,elec_tet_ix),2)==4);

        % gradient{source}(elec_tet_ix,:) = repmat(max(gradient{source}),[length(elec_tet_ix),1]); %assign maximum efield value
        tmp = sort(abs(gradient{source}),'descend');
        gradient{source}(elec_tet_ix,:) = repmat(mean(tmp(1:ceil(length(tmp(:,1))*0.001),:)),[length(elec_tet_ix),1]); % choose mean of highest 0.1% as new efield value
        clear tmp
    else % empty source..
        gradient{source}=zeros(size(vol.tet,1),3);
    end
end

gradient=gradient{1}+gradient{2}+gradient{3}+gradient{4}; % combined gradient from all sources.

vol.pos=vol.pos*SIfx; % convert back to mm.

% midpts can be in either MNI or native space depending on options.native
midpts=mean(cat(3,vol.pos(vol.tet(:,1),:),vol.pos(vol.tet(:,2),:),vol.pos(vol.tet(:,3),:),vol.pos(vol.tet(:,4),:)),3); % midpoints of each pyramid

reduc=10;

% generate flowfield visualization:
% generate a jittered indices vector to be used to reduce flowfield
% display by ~factor reduc.
ea_dispt('Calculating quiver field of gradient for display purposes...');

indices=zeros(length(1:reduc:length(midpts)),1);
cnt=1;
for idx=1:reduc:length(midpts)
    indices(cnt)=idx+round(randn(1)*(reduc/3));
    cnt=cnt+1;
end
indices=unique(indices(2:end-1));
indices(indices==0)=[];
indices(indices>length(midpts))=[];

if ~options.native % VTA calculated in MNI space directly
    [vatfv,vatvolume,radius]=ea_write_vta_nii(S,stimname,midpts,indices,elspec,actContact,voltix,constvol,thresh,mesh,gradient,side,resultfig,options);
    varargout{1}=vatfv;
    varargout{2}=vatvolume;
    varargout{3}=radius;
    ea_dispt('');
else % VTA calculated in native space and then transformed back to MNI
        % Write out native space VTA
        [vatfv,vatvolume,radius]=ea_write_vta_nii(S,stimname,midpts,indices,elspec,actContact,voltix,constvol,thresh,mesh,gradient,side,resultfig,options);

        % If visualizing in native space -> output native space results first before conversion to MNI
        if options.orignative==1
            varargout{1}=vatfv;
            varargout{2}=vatvolume;
            varargout{3}=radius;
            ea_dispt('');
        end

        % Convert midpts and actContact from native space to MNI space
        ptsvx_native = ea_mm2vox([midpts;actContact], options.subj.preopAnat.(options.subj.AnchorModality).coreg)';
        ptsmm_mni = ea_map_coords(ptsvx_native, options.subj.preopAnat.(options.subj.AnchorModality).coreg, ...
            [options.subj.subjDir,filesep,'inverseTransform'], '')';
        midpts_mni = ptsmm_mni(1:size(midpts,1),:);
        actContact_mni = ptsmm_mni(size(midpts,1)+1:end,:);
        options.native=0; % go back to template space for export
        [vatfv,vatvolume,radius]=ea_write_vta_nii(S,stimname,midpts_mni,indices,elspec,actContact_mni,voltix,constvol,thresh,mesh,gradient,side,resultfig,options);
        options.native=options.orignative; % go back to originally set space

        % If visualizing in MNI space -> output MNI space results
        if options.orignative==0
            varargout{1}=vatfv;
            varargout{2}=vatvolume;
            varargout{3}=radius;
            ea_dispt('');
        end
end


function changed=ea_headmodel_changed(options,side,elstruct)
% function that checked if anything (user settings) has changed and
% headmodel needs to be recalculated..
changed=1; % in doubt always reconstruct headmodel

if isequaln(ea_load_hmprotocol(options,side),ea_save_hmprotocol(options,side,elstruct,0)) % important to use isequaln if not nans are treated as not equal...
    changed=0;
end


function protocol=ea_save_hmprotocol(options,side,elstruct,sv)
% function to construct and/or save protocol.
protocol=struct; % default for errors
protocol.elmodel=options.elmodel;
protocol.elstruct=elstruct;
protocol.usediffusion=options.usediffusion;
protocol.atlas=options.atlasset;
protocol.version=1.1;
protocol.vatsettings=options.prefs.machine.vatsettings;

if sv % save protocol to disk
    headmodelDir = fullfile(options.subj.subjDir, 'headmodel', ea_nt(options));
    filePrefix = ['sub-', options.subj.subjId, '_desc-'];
    save(fullfile(headmodelDir, [filePrefix, 'hmprotocol',num2str(side),'.mat']), 'protocol');
end


function protocol=ea_load_hmprotocol(options,side)
% function that loads protocol
try
    headmodelDir = fullfile(options.subj.subjDir, 'headmodel', ea_nt(options));
    filePrefix = ['sub-', options.subj.subjId, '_desc-'];
    load(fullfile(headmodelDir, [filePrefix, 'hmprotocol',num2str(side),'.mat']));
catch
    protocol=struct; % default for errors or if not present
end


%% begin FieldTrip/SimBio functions:
function gradient = ea_calc_gradient(vol,potential)
normal = cross(vol.pos(vol.tet(:,4),:)-vol.pos(vol.tet(:,3),:),vol.pos(vol.tet(:,3),:)-vol.pos(vol.tet(:,2),:));
gradient = repmat(potential(vol.tet(:,1))./sum(normal.*(vol.pos(vol.tet(:,1),:)-(vol.pos(vol.tet(:,2),:)+vol.pos(vol.tet(:,3),:)+vol.pos(vol.tet(:,4),:))/3),2),1,3).*normal;
normal = cross(vol.pos(vol.tet(:,1),:)-vol.pos(vol.tet(:,4),:),vol.pos(vol.tet(:,4),:)-vol.pos(vol.tet(:,3),:));
gradient = gradient + repmat(potential(vol.tet(:,2))./sum(normal.*(vol.pos(vol.tet(:,2),:)-(vol.pos(vol.tet(:,3),:)+vol.pos(vol.tet(:,4),:)+vol.pos(vol.tet(:,1),:))/3),2),1,3).*normal;
normal = cross(vol.pos(vol.tet(:,2),:)-vol.pos(vol.tet(:,1),:),vol.pos(vol.tet(:,1),:)-vol.pos(vol.tet(:,4),:));
gradient = gradient + repmat(potential(vol.tet(:,3))./sum(normal.*(vol.pos(vol.tet(:,3),:)-(vol.pos(vol.tet(:,4),:)+vol.pos(vol.tet(:,1),:)+vol.pos(vol.tet(:,2),:))/3),2),1,3).*normal;
normal = cross(vol.pos(vol.tet(:,3),:)-vol.pos(vol.tet(:,2),:),vol.pos(vol.tet(:,2),:)-vol.pos(vol.tet(:,1),:));
gradient = gradient + repmat(potential(vol.tet(:,4))./sum(normal.*(vol.pos(vol.tet(:,4),:)-(vol.pos(vol.tet(:,1),:)+vol.pos(vol.tet(:,2),:)+vol.pos(vol.tet(:,3),:))/3),2),1,3).*normal;


function potential = ea_apply_dbs(vol,elec,val,unipolar,constvol,boundarynodes)
if constvol
    if unipolar
        dirinodes = [boundarynodes, elec'];
    else
        dirinodes = elec;
    end

    rhs = zeros(length(vol.pos),1);
    dirival = zeros(size(vol.pos,1),1);
    dirival(elec) = val(:,1);
else
    if unipolar
        dirinodes = boundarynodes;
    else
        dirinodes = 1;
    end
    dirival = zeros(size(vol.pos,1),1);

    rhs = zeros(size(vol.pos,1),1);
    uvals=unique(val(:,2));
    if unipolar && length(uvals)==1
        elec_center_id = ea_find_elec_center(elec,vol.pos);
        rhs(elec_center_id) = val(1,1);
    else
        for v=1:length(uvals)
        elec_center_id = ea_find_elec_center(elec(val(:,2)==uvals(v)),vol.pos);
        thesevals=val(val(:,2)==uvals(v),1);
        rhs(elec_center_id) = thesevals(1);
        end
        %warning('Bipolar constant current stimulation currently not implemented!');
    end
end

[stiff, rhs] = ea_dbs(vol.stiff,rhs,dirinodes,dirival);

potential = ea_sb_solve(stiff,rhs);


function center_id = ea_find_elec_center(elec, pos)

center = mean(pos(elec,:));

dist_center = sqrt(sum((pos(elec,:)-repmat(center,length(elec),1)).^2,2));
[dist, elec_id] = min(dist_center);
center_id = elec(elec_id);


function [stiff,rhs] = ea_dbs(stiff,rhs,dirinodes,dirival)

diagonal = diag(stiff);
stiff = stiff + stiff';
rhs = rhs - stiff*dirival;
stiff(dirinodes,:) = 0.0;
stiff(:,dirinodes) = 0.0;
diagonal = -diagonal;
%diagonal(:) = 0;
diagonal(dirinodes) = 1.0;
stiff = stiff + spdiags(diagonal(:),0,length(diagonal),length(diagonal));
rhs(dirinodes) = dirival(dirinodes);


function [warped] = ea_ft_warp_apply(M, input, method, tol)

% ea_ft_warp_apply performs a 3D linear or nonlinear transformation on the input
% coordinates, similar to those in AIR 3.08. You can find technical
% documentation on warping in general at http://bishopw.loni.ucla.edu/AIR3
%
% Use as
%   [warped] = ea_ft_warp_apply(M, input, method, tol)
% where
%   M        vector or matrix with warping parameters
%   input    Nx3 matrix with coordinates
%   warped   Nx3 matrix with coordinates
%   method   string describing the warping method
%   tol      (optional) value determining the numerical precision of the
%             output, to deal with numerical round off imprecisions due to
%             the warping
%
% The methods 'nonlin0', 'nonlin2' ... 'nonlin5' specify a
% polynomial transformation. The size of the transformation matrix
% depends on the order of the warp
%   zeroth order :  1 parameter  per coordinate (translation)
%   first  order :  4 parameters per coordinate (total 12, affine)
%   second order : 10 parameters per coordinate
%   third  order : 20 parameters per coordinate
%   fourth order : 35 parameters per coordinate
%   fifth  order : 56 parameters per coordinate (total 168)
% The size of M should be 3xP, where P is the number of parameters
% per coordinate. Alternatively, you can specify the method to be
% 'nonlinear', where the order will be determined from the size of
% the matrix M.
%
% If the method 'homogeneous' is selected, the input matrix M should be
% a 4x4 homogenous transformation matrix.
%
% If the method 'sn2individual' or 'individual2sn' is selected, the input
% M should be a structure based on nonlinear (warping) normalisation parameters
% created by SPM8 for alignment between an individual structural MRI and the
% template MNI brain.  These options call private functions of the same name.
% M will have subfields like this:
%     Affine: [4x4 double]
%         Tr: [4-D double]
%         VF: [1x1 struct]
%         VG: [1x1 struct]
%      flags: [1x1 struct]
%
% If any other method is selected, it is assumed that it specifies
% the name of an auxiliary function that will, when given the input
% parameter vector M, return an 4x4 homogenous transformation
% matrix. Supplied functions in the warping toolbox are translate,
% rotate, scale, rigidbody, globalrescale, traditional, affine,
% perspective.
%
% See also FT_WARP_OPTIM, FT_WARP_ERROR

% Copyright (C) 2000-2013, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ea_ft_warp_apply.m 10132 2015-01-27 16:08:29Z johzum $

if nargin<4
    tol = [];
end

if nargin<3 && all(size(M)==4)
    % no specific transformation mode has been selected
    % it looks like a homogenous transformation matrix
    method = 'homogeneous';
elseif nargin<3
    % the default method is 'nonlinear'
    method = 'nonlinear';
end

if size(input,2)==2
    % convert the input points from 2D to 3D representation
    input(:,3) = 0;
    input3d = false;
else
    input3d = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nonlinear warping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(strcmp(method, {'nonlinear', 'nonlin0', 'nonlin1', 'nonlin2', 'nonlin3', 'nonlin4', 'nonlin5'}))
    x = input(:,1);
    y = input(:,2);
    z = input(:,3);
    s = size(M);

    if s(1)~=3
        error('invalid size of nonlinear transformation matrix');
    elseif strcmp(method, 'nonlin0') && s(2)~=1
        error('invalid size of nonlinear transformation matrix');
    elseif strcmp(method, 'nonlin1') && s(2)~=4
        error('invalid size of nonlinear transformation matrix');
    elseif strcmp(method, 'nonlin2') && s(2)~=10
        error('invalid size of nonlinear transformation matrix');
    elseif strcmp(method, 'nonlin3') && s(2)~=20
        error('invalid size of nonlinear transformation matrix');
    elseif strcmp(method, 'nonlin4') && s(2)~=35
        error('invalid size of nonlinear transformation matrix');
    elseif strcmp(method, 'nonlin5') && s(2)~=56
        error('invalid size of nonlinear transformation matrix');
    end

    if s(2)==1
        % this is a translation, which in a strict sense is not the 0th order nonlinear transformation
        xx = M(1,1) + x;
        yy = M(2,1) + y;
        zz = M(3,1) + z;
    elseif s(2)==4
        xx = M(1,1) + M(1,2)*x + M(1,3)*y + M(1,4)*z;
        yy = M(2,1) + M(2,2)*x + M(2,3)*y + M(2,4)*z;
        zz = M(3,1) + M(3,2)*x + M(3,3)*y + M(3,4)*z;
    elseif s(2)==10
        xx = M(1,1) + M(1,2)*x + M(1,3)*y + M(1,4)*z + M(1,5)*x.*x + M(1,6)*x.*y + M(1,7)*x.*z + M(1,8)*y.*y + M(1,9)*y.*z + M(1,10)*z.*z;
        yy = M(2,1) + M(2,2)*x + M(2,3)*y + M(2,4)*z + M(2,5)*x.*x + M(2,6)*x.*y + M(2,7)*x.*z + M(2,8)*y.*y + M(2,9)*y.*z + M(2,10)*z.*z;
        zz = M(3,1) + M(3,2)*x + M(3,3)*y + M(3,4)*z + M(3,5)*x.*x + M(3,6)*x.*y + M(3,7)*x.*z + M(3,8)*y.*y + M(3,9)*y.*z + M(3,10)*z.*z;
    elseif s(2)==20
        xx = M(1,1) + M(1,2)*x + M(1,3)*y + M(1,4)*z + M(1,5)*x.*x + M(1,6)*x.*y + M(1,7)*x.*z + M(1,8)*y.*y + M(1,9)*y.*z + M(1,10)*z.*z + M(1,11)*x.*x.*x + M(1,12)*x.*x.*y + M(1,13)*x.*x.*z + M(1,14)*x.*y.*y + M(1,15)*x.*y.*z + M(1,16)*x.*z.*z + M(1,17)*y.*y.*y + M(1,18)*y.*y.*z + M(1,19)*y.*z.*z + M(1,20)*z.*z.*z;
        yy = M(2,1) + M(2,2)*x + M(2,3)*y + M(2,4)*z + M(2,5)*x.*x + M(2,6)*x.*y + M(2,7)*x.*z + M(2,8)*y.*y + M(2,9)*y.*z + M(2,10)*z.*z + M(2,11)*x.*x.*x + M(2,12)*x.*x.*y + M(2,13)*x.*x.*z + M(2,14)*x.*y.*y + M(2,15)*x.*y.*z + M(2,16)*x.*z.*z + M(2,17)*y.*y.*y + M(2,18)*y.*y.*z + M(2,19)*y.*z.*z + M(2,20)*z.*z.*z;
        zz = M(3,1) + M(3,2)*x + M(3,3)*y + M(3,4)*z + M(3,5)*x.*x + M(3,6)*x.*y + M(3,7)*x.*z + M(3,8)*y.*y + M(3,9)*y.*z + M(3,10)*z.*z + M(3,11)*x.*x.*x + M(3,12)*x.*x.*y + M(3,13)*x.*x.*z + M(3,14)*x.*y.*y + M(3,15)*x.*y.*z + M(3,16)*x.*z.*z + M(3,17)*y.*y.*y + M(3,18)*y.*y.*z + M(3,19)*y.*z.*z + M(3,20)*z.*z.*z;
    elseif s(2)==35
        xx = M(1,1) + M(1,2)*x + M(1,3)*y + M(1,4)*z + M(1,5)*x.*x + M(1,6)*x.*y + M(1,7)*x.*z + M(1,8)*y.*y + M(1,9)*y.*z + M(1,10)*z.*z + M(1,11)*x.*x.*x + M(1,12)*x.*x.*y + M(1,13)*x.*x.*z + M(1,14)*x.*y.*y + M(1,15)*x.*y.*z + M(1,16)*x.*z.*z + M(1,17)*y.*y.*y + M(1,18)*y.*y.*z + M(1,19)*y.*z.*z + M(1,20)*z.*z.*z + M(1,21)*x.*x.*x.*x + M(1,22)*x.*x.*x.*y + M(1,23)*x.*x.*x.*z + M(1,24)*x.*x.*y.*y + M(1,25)*x.*x.*y.*z + M(1,26)*x.*x.*z.*z + M(1,27)*x.*y.*y.*y + M(1,28)*x.*y.*y.*z + M(1,29)*x.*y.*z.*z + M(1,30)*x.*z.*z.*z + M(1,31)*y.*y.*y.*y + M(1,32)*y.*y.*y.*z + M(1,33)*y.*y.*z.*z + M(1,34)*y.*z.*z.*z + M(1,35)*z.*z.*z.*z;
        yy = M(2,1) + M(2,2)*x + M(2,3)*y + M(2,4)*z + M(2,5)*x.*x + M(2,6)*x.*y + M(2,7)*x.*z + M(2,8)*y.*y + M(2,9)*y.*z + M(2,10)*z.*z + M(2,11)*x.*x.*x + M(2,12)*x.*x.*y + M(2,13)*x.*x.*z + M(2,14)*x.*y.*y + M(2,15)*x.*y.*z + M(2,16)*x.*z.*z + M(2,17)*y.*y.*y + M(2,18)*y.*y.*z + M(2,19)*y.*z.*z + M(2,20)*z.*z.*z + M(2,21)*x.*x.*x.*x + M(2,22)*x.*x.*x.*y + M(2,23)*x.*x.*x.*z + M(2,24)*x.*x.*y.*y + M(2,25)*x.*x.*y.*z + M(2,26)*x.*x.*z.*z + M(2,27)*x.*y.*y.*y + M(2,28)*x.*y.*y.*z + M(2,29)*x.*y.*z.*z + M(2,30)*x.*z.*z.*z + M(2,31)*y.*y.*y.*y + M(2,32)*y.*y.*y.*z + M(2,33)*y.*y.*z.*z + M(2,34)*y.*z.*z.*z + M(2,35)*z.*z.*z.*z;
        zz = M(3,1) + M(3,2)*x + M(3,3)*y + M(3,4)*z + M(3,5)*x.*x + M(3,6)*x.*y + M(3,7)*x.*z + M(3,8)*y.*y + M(3,9)*y.*z + M(3,10)*z.*z + M(3,11)*x.*x.*x + M(3,12)*x.*x.*y + M(3,13)*x.*x.*z + M(3,14)*x.*y.*y + M(3,15)*x.*y.*z + M(3,16)*x.*z.*z + M(3,17)*y.*y.*y + M(3,18)*y.*y.*z + M(3,19)*y.*z.*z + M(3,20)*z.*z.*z + M(3,21)*x.*x.*x.*x + M(3,22)*x.*x.*x.*y + M(3,23)*x.*x.*x.*z + M(3,24)*x.*x.*y.*y + M(3,25)*x.*x.*y.*z + M(3,26)*x.*x.*z.*z + M(3,27)*x.*y.*y.*y + M(3,28)*x.*y.*y.*z + M(3,29)*x.*y.*z.*z + M(3,30)*x.*z.*z.*z + M(3,31)*y.*y.*y.*y + M(3,32)*y.*y.*y.*z + M(3,33)*y.*y.*z.*z + M(3,34)*y.*z.*z.*z + M(3,35)*z.*z.*z.*z;
    elseif s(2)==56
        xx = M(1,1) + M(1,2)*x + M(1,3)*y + M(1,4)*z + M(1,5)*x.*x + M(1,6)*x.*y + M(1,7)*x.*z + M(1,8)*y.*y + M(1,9)*y.*z + M(1,10)*z.*z + M(1,11)*x.*x.*x + M(1,12)*x.*x.*y + M(1,13)*x.*x.*z + M(1,14)*x.*y.*y + M(1,15)*x.*y.*z + M(1,16)*x.*z.*z + M(1,17)*y.*y.*y + M(1,18)*y.*y.*z + M(1,19)*y.*z.*z + M(1,20)*z.*z.*z + M(1,21)*x.*x.*x.*x + M(1,22)*x.*x.*x.*y + M(1,23)*x.*x.*x.*z + M(1,24)*x.*x.*y.*y + M(1,25)*x.*x.*y.*z + M(1,26)*x.*x.*z.*z + M(1,27)*x.*y.*y.*y + M(1,28)*x.*y.*y.*z + M(1,29)*x.*y.*z.*z + M(1,30)*x.*z.*z.*z + M(1,31)*y.*y.*y.*y + M(1,32)*y.*y.*y.*z + M(1,33)*y.*y.*z.*z + M(1,34)*y.*z.*z.*z + M(1,35)*z.*z.*z.*z + M(1,36)*x.*x.*x.*x.*x + M(1,37)*x.*x.*x.*x.*y + M(1,38)*x.*x.*x.*x.*z + M(1,39)*x.*x.*x.*y.*y + M(1,40)*x.*x.*x.*y.*z + M(1,41)*x.*x.*x.*z.*z + M(1,42)*x.*x.*y.*y.*y + M(1,43)*x.*x.*y.*y.*z + M(1,44)*x.*x.*y.*z.*z + M(1,45)*x.*x.*z.*z.*z + M(1,46)*x.*y.*y.*y.*y + M(1,47)*x.*y.*y.*y.*z + M(1,48)*x.*y.*y.*z.*z + M(1,49)*x.*y.*z.*z.*z + M(1,50)*x.*z.*z.*z.*z + M(1,51)*y.*y.*y.*y.*y + M(1,52)*y.*y.*y.*y.*z + M(1,53)*y.*y.*y.*z.*z + M(1,54)*y.*y.*z.*z.*z + M(1,55)*y.*z.*z.*z.*z + M(1,56)*z.*z.*z.*z.*z;
        yy = M(2,1) + M(2,2)*x + M(2,3)*y + M(2,4)*z + M(2,5)*x.*x + M(2,6)*x.*y + M(2,7)*x.*z + M(2,8)*y.*y + M(2,9)*y.*z + M(2,10)*z.*z + M(2,11)*x.*x.*x + M(2,12)*x.*x.*y + M(2,13)*x.*x.*z + M(2,14)*x.*y.*y + M(2,15)*x.*y.*z + M(2,16)*x.*z.*z + M(2,17)*y.*y.*y + M(2,18)*y.*y.*z + M(2,19)*y.*z.*z + M(2,20)*z.*z.*z + M(2,21)*x.*x.*x.*x + M(2,22)*x.*x.*x.*y + M(2,23)*x.*x.*x.*z + M(2,24)*x.*x.*y.*y + M(2,25)*x.*x.*y.*z + M(2,26)*x.*x.*z.*z + M(2,27)*x.*y.*y.*y + M(2,28)*x.*y.*y.*z + M(2,29)*x.*y.*z.*z + M(2,30)*x.*z.*z.*z + M(2,31)*y.*y.*y.*y + M(2,32)*y.*y.*y.*z + M(2,33)*y.*y.*z.*z + M(2,34)*y.*z.*z.*z + M(2,35)*z.*z.*z.*z + M(2,36)*x.*x.*x.*x.*x + M(2,37)*x.*x.*x.*x.*y + M(2,38)*x.*x.*x.*x.*z + M(2,39)*x.*x.*x.*y.*y + M(2,40)*x.*x.*x.*y.*z + M(2,41)*x.*x.*x.*z.*z + M(2,42)*x.*x.*y.*y.*y + M(2,43)*x.*x.*y.*y.*z + M(2,44)*x.*x.*y.*z.*z + M(2,45)*x.*x.*z.*z.*z + M(2,46)*x.*y.*y.*y.*y + M(2,47)*x.*y.*y.*y.*z + M(2,48)*x.*y.*y.*z.*z + M(2,49)*x.*y.*z.*z.*z + M(2,50)*x.*z.*z.*z.*z + M(2,51)*y.*y.*y.*y.*y + M(2,52)*y.*y.*y.*y.*z + M(2,53)*y.*y.*y.*z.*z + M(2,54)*y.*y.*z.*z.*z + M(2,55)*y.*z.*z.*z.*z + M(2,56)*z.*z.*z.*z.*z;
        zz = M(3,1) + M(3,2)*x + M(3,3)*y + M(3,4)*z + M(3,5)*x.*x + M(3,6)*x.*y + M(3,7)*x.*z + M(3,8)*y.*y + M(3,9)*y.*z + M(3,10)*z.*z + M(3,11)*x.*x.*x + M(3,12)*x.*x.*y + M(3,13)*x.*x.*z + M(3,14)*x.*y.*y + M(3,15)*x.*y.*z + M(3,16)*x.*z.*z + M(3,17)*y.*y.*y + M(3,18)*y.*y.*z + M(3,19)*y.*z.*z + M(3,20)*z.*z.*z + M(3,21)*x.*x.*x.*x + M(3,22)*x.*x.*x.*y + M(3,23)*x.*x.*x.*z + M(3,24)*x.*x.*y.*y + M(3,25)*x.*x.*y.*z + M(3,26)*x.*x.*z.*z + M(3,27)*x.*y.*y.*y + M(3,28)*x.*y.*y.*z + M(3,29)*x.*y.*z.*z + M(3,30)*x.*z.*z.*z + M(3,31)*y.*y.*y.*y + M(3,32)*y.*y.*y.*z + M(3,33)*y.*y.*z.*z + M(3,34)*y.*z.*z.*z + M(3,35)*z.*z.*z.*z + M(3,36)*x.*x.*x.*x.*x + M(3,37)*x.*x.*x.*x.*y + M(3,38)*x.*x.*x.*x.*z + M(3,39)*x.*x.*x.*y.*y + M(3,40)*x.*x.*x.*y.*z + M(3,41)*x.*x.*x.*z.*z + M(3,42)*x.*x.*y.*y.*y + M(3,43)*x.*x.*y.*y.*z + M(3,44)*x.*x.*y.*z.*z + M(3,45)*x.*x.*z.*z.*z + M(3,46)*x.*y.*y.*y.*y + M(3,47)*x.*y.*y.*y.*z + M(3,48)*x.*y.*y.*z.*z + M(3,49)*x.*y.*z.*z.*z + M(3,50)*x.*z.*z.*z.*z + M(3,51)*y.*y.*y.*y.*y + M(3,52)*y.*y.*y.*y.*z + M(3,53)*y.*y.*y.*z.*z + M(3,54)*y.*y.*z.*z.*z + M(3,55)*y.*z.*z.*z.*z + M(3,56)*z.*z.*z.*z.*z;
    else
        error('invalid size of nonlinear transformation matrix');
    end

    warped = [xx yy zz];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % linear warping using homogenous coordinate transformation matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(method, 'homogenous') || strcmp(method, 'homogeneous')
    if all(size(M)==3)
        % convert the 3x3 homogenous transformation matrix (corresponding with 2D)
        % into a 4x4 homogenous transformation matrix (corresponding with 3D)
        M = [
            M(1,1) M(1,2)  0  M(1,3)
            M(2,1) M(2,2)  0  M(2,3)
            0      0       0  0
            M(3,1) M(3,2)  0  M(3,3)
            ];
    end

    %warped = M * [input'; ones(1, size(input, 1))];
    %warped = warped(1:3,:)';

    % below achieves the same as lines 154-155
    warped = [input ones(size(input, 1),1)]*M(1:3,:)';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % using external function that returns a homogeneous transformation matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif exist(method, 'file') && ~isa(M, 'struct')
    % get the homogenous transformation matrix
    H = feval(method, M);
    warped = ea_ft_warp_apply(H, input, 'homogeneous');

elseif strcmp(method, 'sn2individual') && isa(M, 'struct')
    % use SPM structure with parameters for an inverse warp
    % from normalized space to individual, can be non-linear
    warped = sn2individual(M, input);

elseif strcmp(method, 'individual2sn') && isa(M, 'struct')
    % use SPM structure with parameters for a warp from
    % individual space to normalized space, can be non-linear
    %error('individual2sn is not yet implemented');
    warped = individual2sn(M, input);
else
    error('unrecognized transformation method');
end

if ~input3d
    % convert from 3D back to 2D representation
    warped = warped(:,1:2);
end

if ~isempty(tol)
    if tol>0
        warped = fix(warped./tol)*tol;
    end
end


function x = ea_sb_solve(sysmat,vecb)

% SB_SOLVE
%
% $Id: sb_solve.m 8776 2013-11-14 09:04:48Z roboos $
try
    L = ichol(sysmat);
catch
    alpha = max(sum(abs(sysmat),2)./diag(sysmat))-2;
    L = ichol(sysmat, struct('type','ict','droptol',1e-3,'diagcomp',alpha));
end

%scalen
[~,x]=evalc('pcg(sysmat,vecb,10e-10,5000,L,L'',vecb)');


function [type] = ea_ft_voltype(vol, desired)

% FT_VOLTYPE determines the type of volume conduction model of the head
%
% Use as
%   [type] = ft_voltype(vol)
% to get a string describing the type, or
%   [flag] = ft_voltype(vol, desired)
% to get a boolean value.
%
% For EEG the following volume conduction models are recognized
%   singlesphere       analytical single sphere model
%   concentricspheres  analytical concentric sphere model with up to 4 spheres
%   halfspace          infinite homogenous medium on one side, vacuum on the other
%   openmeeg           boundary element method, based on the OpenMEEG software
%   bemcp              boundary element method, based on the implementation from Christophe Phillips
%   dipoli             boundary element method, based on the implementation from Thom Oostendorp
%   asa                boundary element method, based on the (commercial) ASA software
%   simbio             finite element method, based on the SimBio software
%   fns                finite difference method, based on the FNS software
%   interpolate        interpolate the potential based on pre-computed leadfields
%
% and for MEG the following volume conduction models are recognized
%   singlesphere       analytical single sphere model
%   localspheres       local spheres model for MEG, one sphere per channel
%   singleshell        realisically shaped single shell approximation, based on the implementation from Guido Nolte
%   infinite           magnetic dipole in an infinite vacuum
%   interpolate        interpolate the potential based on pre-computed leadfields
%
% See also FT_COMPUTE_LEADFIELD, FT_READ_VOL, FT_HEADMODEL_BEMCP,
% FT_HEADMODEL_ASA, FT_HEADMODEL_DIPOLI, FT_HEADMODEL_SIMBIO,
% FT_HEADMODEL_FNS, FT_HEADMODEL_HALFSPACE, FT_HEADMODEL_INFINITE,
% FT_HEADMODEL_OPENMEEG, FT_HEADMODEL_SINGLESPHERE,
% FT_HEADMODEL_CONCENTRICSPHERES, FT_HEADMODEL_LOCALSPHERES,
% FT_HEADMODEL_SINGLESHELL, FT_HEADMODEL_INTERPOLATE

% Copyright (C) 2007-2013, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ft_voltype.m 10012 2014-12-03 09:14:18Z roboos $

% these are for remembering the type on subsequent calls with the same input arguments
persistent previous_argin previous_argout

if iscell(vol) && numel(vol)<4
    % this might represent combined EEG, ECoG and/or MEG
    type = cell(size(vol));
    if nargin<2
        desired = cell(size(vol)); % empty elements
    end
    for i=1:numel(vol)
        type{i} = ea_ft_voltype(vol{i}, desired{i});
    end
    return
end

if nargin<2
    % ensure that all input arguments are defined
    desired = [];
end

current_argin = {vol, desired};
if isequal(current_argin, previous_argin)
    % don't do the type detection again, but return the previous values from cache
    type = previous_argout{1};
    return
end

if isfield(vol, 'type') && ~(ea_ft_datatype(vol, 'grad') || ea_ft_datatype(vol, 'sens')) % grad and sens also contain .type fields
    % preferably the structure specifies its own type
    type = vol.type;

elseif isfield(vol, 'r') && numel(vol.r)==1 && ~isfield(vol, 'label')
    type = 'singlesphere';

elseif isfield(vol, 'r') && isfield(vol, 'o') && isfield(vol, 'label')
    % this is before the spheres have been assigned to the coils
    % and every sphere is still associated with a channel
    type = 'localspheres';

elseif isfield(vol, 'r') && isfield(vol, 'o') && size(vol.r,1)==size(vol.o,1) && size(vol.r,1)>4
    % this is after the spheres have been assigned to the coils
    % note that this one is easy to confuse with the concentric one
    type = 'localspheres';

elseif isfield(vol, 'r') && numel(vol.r)>=2 && ~isfield(vol, 'label')
    type = 'concentricspheres';

elseif isfield(vol, 'bnd') && isfield(vol, 'mat')
    type = 'bem'; % it could be dipoli, asa, bemcp or openmeeg

elseif isfield(vol, 'bnd') && isfield(vol, 'forwpar')
    type = 'singleshell';

elseif isfield(vol, 'bnd') && numel(vol.bnd)==1
    type = 'singleshell';

elseif isempty(vol) || (isstruct(vol) && isequal(fieldnames(vol), {'unit'}))
    % it is empty, or only contains a specification of geometrical units
    type = 'infinite';

else
    type = 'unknown';

end % if isfield(vol, 'type')

if ~isempty(desired)
    % return a boolean flag
    switch desired
        case 'bem'
            type = any(strcmp(type, {'bem', 'dipoli', 'asa', 'bemcp', 'openmeeg'}));
        otherwise
            type = any(strcmp(type, desired));
    end % switch desired
end % determine the correspondence to the desired type

% remember the current input and output arguments, so that they can be
% reused on a subsequent call in case the same input argument is given
current_argout  = {type};
previous_argin  = current_argin;
previous_argout = current_argout;

return % voltype main()

function sens = ea_undobalancing(sens)

% UNDOBALANCING removes all balancing coefficients from the gradiometer sensor array
%
% This is used in CHANNELPOSITION, FT_PREPARE_LAYOUT, FT_SENSTYPE

while isfield(sens, 'balance') && isfield(sens.balance, 'current') && ~strcmp(sens.balance.current, 'none')
    fnames = setdiff(fieldnames(sens.balance), 'current');
    indx   = find(ismember(fnames, sens.balance.current));

    if length(indx)==1
        % undo the synthetic gradient balancing
        fprintf('undoing the %s balancing for the gradiometer definition\n', sens.balance.current);

        % if componentanalysis was followed by rejectcomponent, the balancing matrix is rank deficient
        % leading to problems in the correct allocation of the coils to the channels
        if strcmp(sens.balance.current, 'invcomp') && strcmp(sens.balance.previous{1}, 'comp')
            tra1 = full(sens.balance.invcomp.tra);
            tra2 = full(sens.balance.comp.tra);
            tra3 = tra1;
            tmp  = tra1*tra2;
            tmp  = null(tmp); % nullspace after componentanalysis and rejectcomponent
            tmp  = tmp*tmp';  % this is the part which was removed at some point
            [ix,iy]     = match_str(sens.balance.comp.labelorg, sens.balance.invcomp.labelnew);
            tra3(iy,iy) = (eye(numel(ix))+tmp(ix,ix))*tra1(iy,iy);
            sens.balance.invcomp.tra = tra3;
            % FIXME check whether this is robust
        end

        if strcmp(sens.balance.current, 'planar')
            if isfield(sens, 'type') && contains(sens.type, '_planar')
                % remove the planar postfox from the sensor type
                sens.type = sens.type(1:(end-7));
            end
        end

        sens = ft_apply_montage(sens, sens.balance.(sens.balance.current), 'inverse', 'yes', 'keepunused', 'yes', 'warning', 'no');

        if ~isfield(sens, 'chanpos') || any(isnan(sens.chanpos(:))) || any(isnan(sens.chanori(:)))
            % this happens if the data has been component-analyzed
            % try to reconstruct the channel position and orientation
            [pos, ori, lab] = channelposition(sens);
            [sel1, sel2] = match_str(sens.label, lab);
            sens.chanpos(sel1,:) = pos(sel2,:);
            sens.chanori(sel1,:) = ori(sel2,:);
        end

    else
        warning('cannot undo %s balancing in the gradiometer definition\n', sens.balance.current);
        break
    end
end

function [data] = ea_ft_checkdata(data, varargin)

% FT_CHECKDATA checks the input data of the main FieldTrip functions, e.g. whether
% the type of data strucure corresponds with the required data. If neccessary
% and possible, this function will adjust the data structure to the input
% requirements (e.g. change dimord, average over trials, convert inside from
% index into logical).
%
% If the input data does NOT correspond to the requirements, this function
% is supposed to give a elaborate warning message and if applicable point
% the user to external documentation (link to website).
%
% Use as
%   [data] = ft_checkdata(data, ...)
%
% Optional input arguments should be specified as key-value pairs and can include
%   feedback           = yes, no
%   datatype           = raw, freq, timelock, comp, spike, source,  dip, volume, segmentation, parcellation
%   dimord             = any combination of time, freq, chan, refchan, rpt, subj, chancmb, rpttap, pos
%   senstype           = ctf151, ctf275, ctf151_planar, ctf275_planar, neuromag122, neuromag306, bti148, bti248, bti248_planar, magnetometer, electrode
%   inside             = logical, index
%   ismeg              = yes, no
%   hasunit            = yes, no
%   hascoordsys        = yes, no
%   hassampleinfo      = yes, no, ifmakessense (only applies to raw data)
%   hascumtapcnt       = yes, no (only applies to freq data)
%   hasdim             = yes, no
%   hasdof             = yes, no
%   cmbrepresentation  = sparse, full (applies to covariance and cross-spectral density)
%   fsample            = sampling frequency to use to go from SPIKE to RAW representation
%   segmentationstyle  = indexed, probabilistic (only applies to segmentation)
%   parcellationstyle  = indexed, probabilistic (only applies to parcellation)
%   hasbrain           = yes, no (only applies to segmentation)
%
% For some options you can specify multiple values, e.g.
%   [data] = ft_checkdata(data, 'senstype', {'ctf151', 'ctf275'}), e.g. in megrealign
%   [data] = ft_checkdata(data, 'datatype', {'timelock', 'freq'}), e.g. in sourceanalysis

% Copyright (C) 2007-2015, Robert Oostenveld
% Copyright (C) 2010-2012, Martin Vinck
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ft_checkdata.m 10287 2015-03-28 13:04:14Z roboos $

% in case of an error this function could use dbstack for more detailled
% user feedback
%
% this function should replace/encapsulate
%   fixdimord
%   fixinside
%   fixprecision
%   fixvolume
%   data2raw
%   raw2data
%   grid2transform
%   transform2grid
%   fourier2crsspctrm
%   freq2cumtapcnt
%   sensortype
%   time2offset
%   offset2time
%   fixsens -> this is kept a separate function because it should also be
%              called from other modules
%
% other potential uses for this function:
%   time -> offset in freqanalysis
%   average over trials
%   csd as matrix

% FIXME the following is difficult, if not impossible, to support without knowing the parameter
% FIXME it is presently (dec 2014) not being used anywhere in FT, so can be removed
%   hastrials          = yes, no

% get the optional input arguments
feedback             = ea_ft_getopt(varargin, 'feedback', 'no');
dtype                = ea_ft_getopt(varargin, 'datatype'); % should not conflict with the ea_ft_datatype function
dimord               = ea_ft_getopt(varargin, 'dimord');
stype                = ea_ft_getopt(varargin, 'senstype'); % senstype is a function name which should not be masked
ismeg                = ea_ft_getopt(varargin, 'ismeg');
inside               = ea_ft_getopt(varargin, 'inside'); % can be 'logical' or 'index'
hastrials            = ea_ft_getopt(varargin, 'hastrials');
hasunit              = ea_ft_getopt(varargin, 'hasunit', 'no');
hascoordsys          = ea_ft_getopt(varargin, 'hascoordsys', 'no');
hassampleinfo        = ea_ft_getopt(varargin, 'hassampleinfo', 'ifmakessense');
hasdimord            = ea_ft_getopt(varargin, 'hasdimord', 'no');
hasdim               = ea_ft_getopt(varargin, 'hasdim');
hascumtapcnt         = ea_ft_getopt(varargin, 'hascumtapcnt');
hasdof               = ea_ft_getopt(varargin, 'hasdof', 'no');
haspow               = ea_ft_getopt(varargin, 'haspow', 'no');
cmbrepresentation    = ea_ft_getopt(varargin, 'cmbrepresentation');
channelcmb           = ea_ft_getopt(varargin, 'channelcmb');
sourcedimord         = ea_ft_getopt(varargin, 'sourcedimord');
sourcerepresentation = ea_ft_getopt(varargin, 'sourcerepresentation');
fsample              = ea_ft_getopt(varargin, 'fsample');
segmentationstyle    = ea_ft_getopt(varargin, 'segmentationstyle'); % this will be passed on to the corresponding ft_datatype_xxx function
parcellationstyle    = ea_ft_getopt(varargin, 'parcellationstyle'); % this will be passed on to the corresponding ft_datatype_xxx function
hasbrain             = ea_ft_getopt(varargin, 'hasbrain');

% check whether people are using deprecated stuff
depHastrialdef = ea_ft_getopt(varargin, 'hastrialdef');
if (~isempty(depHastrialdef))
    warning_once('ft_checkdata option ''hastrialdef'' is deprecated; use ''hassampleinfo'' instead');
    hassampleinfo = depHastrialdef;
end
if (~isempty(ea_ft_getopt(varargin, 'hasoffset')))
    warning_once('ft_checkdata option ''hasoffset'' has been removed and will be ignored');
end

% determine the type of input data
% this can be raw, freq, timelock, comp, spike, source, volume, dip
israw           = ea_ft_datatype(data, 'raw');
isfreq          = ea_ft_datatype(data, 'freq');
istimelock      = ea_ft_datatype(data, 'timelock');
iscomp          = ea_ft_datatype(data, 'comp');
isspike         = ea_ft_datatype(data, 'spike');
isvolume        = ea_ft_datatype(data, 'volume');
issegmentation  = ea_ft_datatype(data, 'segmentation');
isparcellation  = ea_ft_datatype(data, 'parcellation');
issource        = ea_ft_datatype(data, 'source');
isdip           = ea_ft_datatype(data, 'dip');
ismvar          = ea_ft_datatype(data, 'mvar');
isfreqmvar      = ea_ft_datatype(data, 'freqmvar');
ischan          = ea_ft_datatype(data, 'chan');
% FIXME use the istrue function on ismeg and hasxxx options

if ~isequal(feedback, 'no')
    if iscomp
        % it can be comp and raw/timelock/freq at the same time, therefore this has to go first
        nchan = size(data.topo,1);
        ncomp = size(data.topo,2);
        fprintf('the input is component data with %d components and %d original channels\n', ncomp, nchan);
    end

    if israw
        nchan = length(data.label);
        ntrial = length(data.trial);
        fprintf('the input is raw data with %d channels and %d trials\n', nchan, ntrial);
    elseif istimelock
        nchan = length(data.label);
        ntime = length(data.time);
        fprintf('the input is timelock data with %d channels and %d timebins\n', nchan, ntime);
    elseif isfreq
        if isfield(data, 'label')
            nchan = length(data.label);
            nfreq = length(data.freq);
            if isfield(data, 'time'), ntime = num2str(length(data.time)); else ntime = 'no'; end
            fprintf('the input is freq data with %d channels, %d frequencybins and %s timebins\n', nchan, nfreq, ntime);
        elseif isfield(data, 'labelcmb')
            nchan = length(data.labelcmb);
            nfreq = length(data.freq);
            if isfield(data, 'time'), ntime = num2str(length(data.time)); else ntime = 'no'; end
            fprintf('the input is freq data with %d channel combinations, %d frequencybins and %s timebins\n', nchan, nfreq, ntime);
        else
            error('cannot infer freq dimensions');
        end
    elseif isspike
        nchan  = length(data.label);
        fprintf('the input is spike data with %d channels\n', nchan);
    elseif isvolume
        if issegmentation
            subtype = 'segmented volume';
        else
            subtype = 'volume';
        end
        fprintf('the input is %s data with dimensions [%d %d %d]\n', subtype, data.dim(1), data.dim(2), data.dim(3));
        clear subtype
    elseif issource
        nsource = size(data.pos, 1);
        if isparcellation
            subtype = 'parcellated source';
        else
            subtype = 'source';
        end
        if isfield(data, 'dim')
            fprintf('the input is %s data with %d brainordinates on a [%d %d %d] grid\n', subtype, nsource, data.dim(1), data.dim(2), data.dim(3));
        elseif isfield(data, 'tri')
            fprintf('the input is %s data with %d vertex positions and %d triangles\n', subtype, nsource, size(data.tri, 1));
        else
            fprintf('the input is %s data with %d brainordinates\n', subtype, nsource);
        end
        clear subtype
    elseif isdip
        fprintf('the input is dipole data\n');
    elseif ismvar
        fprintf('the input is mvar data\n');
    elseif isfreqmvar
        fprintf('the input is freqmvar data\n');
    elseif ischan
        nchan = length(data.label);
        if isfield(data, 'brainordinate')
            fprintf('the input is parcellated data with %d parcels\n', nchan);
        else
            fprintf('the input is chan data with %d channels\n', nchan);
        end
    end
end % give feedback

if issource && isvolume
    % it should be either one or the other: the choice here is to represent it as volume description since that is simpler to handle
    % the conversion is done by removing the grid positions
    data = rmfield(data, 'pos');
    issource = false;
end

% the ft_datatype_XXX functions ensures the consistency of the XXX datatype
% and provides a detailed description of the dataformat and its history
if iscomp % this should go before israw/istimelock/isfreq
    data = ft_datatype_comp(data, 'hassampleinfo', hassampleinfo);
elseif israw
    data = ft_datatype_raw(data, 'hassampleinfo', hassampleinfo);
elseif istimelock
    data = ft_datatype_timelock(data);
elseif isfreq
    data = ft_datatype_freq(data);
elseif isspike
    data = ft_datatype_spike(data);
elseif issegmentation % this should go before isvolume
    data = ea_ft_datatype_segmentation(data, 'segmentationstyle', segmentationstyle, 'hasbrain', hasbrain);
elseif isvolume
    data = ft_datatype_volume(data);
elseif isparcellation % this should go before issource
    data = ft_datatype_parcellation(data, 'parcellationstyle', parcellationstyle);
elseif issource
    data = ft_datatype_source(data);
elseif isdip
    data = ft_datatype_dip(data);
elseif ismvar || isfreqmvar
    data = ft_datatype_mvar(data);
end

if ~isempty(dtype)
    if ~isa(dtype, 'cell')
        dtype = {dtype};
    end

    okflag = 0;
    for i=1:length(dtype)
        % check that the data matches with one or more of the required ft_datatypes
        switch dtype{i}
            case 'raw+comp'
                okflag = okflag + (israw & iscomp);
            case 'freq+comp'
                okflag = okflag + (isfreq & iscomp);
            case 'timelock+comp'
                okflag = okflag + (istimelock & iscomp);
            case 'raw'
                okflag = okflag + (israw & ~iscomp);
            case 'freq'
                okflag = okflag + (isfreq & ~iscomp);
            case 'timelock'
                okflag = okflag + (istimelock & ~iscomp);
            case 'comp'
                okflag = okflag + (iscomp & ~(israw | istimelock | isfreq));
            case 'spike'
                okflag = okflag + isspike;
            case 'volume'
                okflag = okflag + isvolume;
            case 'source'
                okflag = okflag + issource;
            case 'dip'
                okflag = okflag + isdip;
            case 'mvar'
                okflag = okflag + ismvar;
            case 'freqmvar'
                okflag = okflag + isfreqmvar;
            case 'chan'
                okflag = okflag + ischan;
            case 'segmentation'
                okflag = okflag + issegmentation;
            case 'parcellation'
                okflag = okflag + isparcellation;
        end % switch dtype
    end % for dtype

    % try to convert the data if needed
    for iCell = 1:length(dtype)
        if okflag
            % the requested datatype is specified in descending order of
            % preference (if there is a preference at all), so don't bother
            % checking the rest of the requested data types if we already
            % succeeded in converting
            break;
        end
        if isequal(dtype(iCell), {'parcellation'}) && issegmentation
            data = eavolume2source(data); % segmentation=volume, parcellation=source
            data = ft_datatype_parcellation(data);
            issegmentation = 0;
            isvolume = 0;
            isparcellation = 1;
            issource = 1;
            okflag = 1;
        elseif isequal(dtype(iCell), {'segmentation'}) && isparcellation
            data = ea_source2volume(data); % segmentation=volume, parcellation=source
            data = ea_ft_datatype_segmentation(data);
            isparcellation = 0;
            issource = 0;
            issegmentation = 1;
            isvolume = 1;
            okflag = 1;
        elseif isequal(dtype(iCell), {'source'}) && isvolume
            data = eavolume2source(data);
            data = ft_datatype_source(data);
            isvolume = 0;
            issource = 1;
            okflag = 1;
        elseif isequal(dtype(iCell), {'volume'}) && (ischan || istimelock || isfreq)
            data = ea_parcellated2source(data);
            data = ft_datatype_volume(data);
            ischan = 0;
            isvolume = 1;
            okflag = 1;
        elseif isequal(dtype(iCell), {'source'}) && (ischan || istimelock || isfreq)
            data = ea_parcellated2source(data);
            data = ft_datatype_source(data);
            ischan = 0;
            issource = 1;
            okflag = 1;
        elseif isequal(dtype(iCell), {'volume'}) && issource
            data = ea_source2volume(data);
            data = ft_datatype_volume(data);
            isvolume = 1;
            issource = 0;
            okflag = 1;
        elseif isequal(dtype(iCell), {'raw+comp'}) && istimelock && iscomp
            data = eatimelock2raw(data);
            data = ft_datatype_raw(data, 'hassampleinfo', hassampleinfo);
            istimelock = 0;
            iscomp = 1;
            israw = 1;
            okflag = 1;
        elseif isequal(dtype(iCell), {'raw'}) && issource
            data = ea_source2raw(data);
            data = ft_datatype_raw(data, 'hassampleinfo', hassampleinfo);
            issource = 0;
            israw = 1;
            okflag = 1;
        elseif isequal(dtype(iCell), {'raw'}) && istimelock
            if iscomp
                data = removefields(data, {'topo', 'topolabel', 'unmixing'}); % these fields are not desired
                iscomp = 0;
            end
            data = eatimelock2raw(data);
            data = ft_datatype_raw(data, 'hassampleinfo', hassampleinfo);
            istimelock = 0;
            israw = 1;
            okflag = 1;
        elseif isequal(dtype(iCell), {'comp'}) && israw
            data = keepfields(data, {'label', 'topo', 'topolabel', 'unmixing', 'elec', 'grad', 'cfg'}); % these are the only relevant fields
            data = ft_datatype_comp(data);
            israw = 0;
            iscomp = 1;
            okflag = 1;
        elseif isequal(dtype(iCell), {'comp'}) && istimelock
            data = keepfields(data, {'label', 'topo', 'topolabel', 'unmixing', 'elec', 'grad', 'cfg'}); % these are the only relevant fields
            data = ft_datatype_comp(data);
            istimelock = 0;
            iscomp = 1;
            okflag = 1;
        elseif isequal(dtype(iCell), {'comp'}) && isfreq
            data = keepfields(data, {'label', 'topo', 'topolabel', 'unmixing', 'elec', 'grad', 'cfg'}); % these are the only relevant fields
            data = ft_datatype_comp(data);
            isfreq = 0;
            iscomp = 1;
            okflag = 1;
        elseif isequal(dtype(iCell), {'raw'}) && israw
            if iscomp
                data = removefields(data, {'topo', 'topolabel', 'unmixing'}); % these fields are not desired
                iscomp = 0;
            end
            data = ft_datatype_raw(data);
            okflag = 1;
        elseif isequal(dtype(iCell), {'timelock'}) && istimelock
            if iscomp
                data = removefields(data, {'topo', 'topolabel', 'unmixing'}); % these fields are not desired
                iscomp = 0;
            end
            data = ft_datatype_timelock(data);
            okflag = 1;
        elseif isequal(dtype(iCell), {'freq'}) && isfreq
            if iscomp
                data = removefields(data, {'topo', 'topolabel', 'unmixing'}); % these fields are not desired
                iscomp = 0;
            end
            data = ft_datatype_freq(data);
            okflag = 1;
        elseif isequal(dtype(iCell), {'timelock'}) && israw
            if iscomp
                data = removefields(data, {'topo', 'topolabel', 'unmixing'}); % these fields are not desired
                iscomp = 0;
            end
            data = earaw2timelock(data);
            data = ft_datatype_timelock(data);
            israw = 0;
            istimelock = 1;
            okflag = 1;
        elseif isequal(dtype(iCell), {'raw'}) && isfreq
            if iscomp
                data = removefields(data, {'topo', 'topolabel', 'unmixing'}); % these fields are not desired
                iscomp = 0;
            end
            data = ea_freq2raw(data);
            data = ft_datatype_raw(data, 'hassampleinfo', hassampleinfo);
            isfreq = 0;
            israw = 1;
            okflag = 1;

        elseif isequal(dtype(iCell), {'raw'}) && ischan
            data = ea_chan2timelock(data);
            data = eatimelock2raw(data);
            data = ft_datatype_raw(data);
            ischan = 0;
            israw = 1;
            okflag = 1;
        elseif isequal(dtype(iCell), {'timelock'}) && ischan
            data = ea_chan2timelock(data);
            data = ft_datatype_timelock(data);
            ischan = 0;
            istimelock = 1;
            okflag = 1;
        elseif isequal(dtype(iCell), {'freq'}) && ischan
            data = ea_chan2freq(data);
            data = ft_datatype_freq(data);
            ischan = 0;
            isfreq = 1;
            okflag = 1;
        elseif isequal(dtype(iCell), {'spike'}) && israw
            data = ea_raw2spike(data);
            data = ft_datatype_spike(data);
            israw = 0;
            isspike = 1;
            okflag = 1;
        elseif isequal(dtype(iCell), {'raw'}) && isspike
            data = ea_spike2raw(data,fsample);
            data = ft_datatype_raw(data, 'hassampleinfo', hassampleinfo);
            isspike = 0;
            israw   = 1;
            okflag  = 1;
        end
    end % for iCell

    if ~okflag
        % construct an error message
        if length(dtype)>1
            str = sprintf('%s, ', dtype{1:(end-2)});
            str = sprintf('%s%s or %s', str, dtype{end-1}, dtype{end});
        else
            str = dtype{1};
        end
        error('This function requires %s data as input.', str);
    end % if okflag
end

if ~isempty(dimord)
    if ~isa(dimord, 'cell')
        dimord = {dimord};
    end

    if isfield(data, 'dimord')
        okflag = any(strcmp(data.dimord, dimord));
    else
        okflag = 0;
    end

    if ~okflag
        % construct an error message
        if length(dimord)>1
            str = sprintf('%s, ', dimord{1:(end-2)});
            str = sprintf('%s%s or %s', str, dimord{end-1}, dimord{end});
        else
            str = dimord{1};
        end
        error('This function requires data with a dimord of %s.', str);
    end % if okflag
end

if ~isempty(stype)
    if ~isa(stype, 'cell')
        stype = {stype};
    end

    if isfield(data, 'grad') || isfield(data, 'elec')
        if any(strcmp(ft_senstype(data), stype))
            okflag = 1;
        elseif any(cellfun(@ft_senstype, repmat({data}, size(stype)), stype))
            % this is req uired to detect more general types, such as "meg" or "ctf" rather than "ctf275"
            okflag = 1;
        else
            okflag = 0;
        end
    end

    if ~okflag
        % construct an error message
        if length(stype)>1
            str = sprintf('%s, ', stype{1:(end-2)});
            str = sprintf('%s%s or %s', str, stype{end-1}, stype{end});
        else
            str = stype{1};
        end
        error('This function requires %s data as input, but you are giving %s data.', str, ft_senstype(data));
    end % if okflag
end

if ~isempty(ismeg)
    if isequal(ismeg, 'yes')
        okflag = isfield(data, 'grad');
    elseif isequal(ismeg, 'no')
        okflag = ~isfield(data, 'grad');
    end

    if ~okflag && isequal(ismeg, 'yes')
        error('This function requires MEG data with a ''grad'' field');
    elseif ~okflag && isequal(ismeg, 'no')
        error('This function should not be given MEG data with a ''grad'' field');
    end % if okflag
end

if ~isempty(inside)
    if strcmp(inside, 'index')
        warning('the indexed representation of inside/outside source locations is deprecated');
    end
    % TODO absorb the fixinside function into this code
    data   = fixinside(data, inside);
    okflag = isfield(data, 'inside');

    if ~okflag
        % construct an error message
        error('This function requires data with an ''inside'' field.');
    end % if okflag
end

%if isvolume
%  % ensure consistent dimensions of the volumetric data
%  % reshape each of the volumes that is found into a 3D array
%  param = parameterselection('all', data);
%  dim   = data.dim;
%  for i=1:length(param)
%    tmp  = getsubfield(data, param{i});
%    tmp  = reshape(tmp, dim);
%    data = setsubfield(data, param{i}, tmp);
%  end
%end

if ea_istrue(hasunit) && ~isfield(data, 'unit')
    % calling convert_units with only the input data adds the units without converting
    data = ft_convert_units(data);
end

if ea_istrue(hascoordsys) && ~isfield(data, 'coordsys')
    data = ft_determine_coordsys(data);
end

if issource || isvolume
    % the following section is to make a dimord-consistent representation of
    % volume and source data, taking trials, time and frequency into account
    if isequal(hasdimord, 'yes') && (~isfield(data, 'dimord') || ~strcmp(data.dimord,sourcedimord))

        % determine the size of the data
        if isfield(data, 'dimord')
            dimtok = tokenize(data.dimord, '_');
            if ~isempty(strmatch('time', dimtok)), Ntime = length(data.time); else Ntime = 1; end
            if ~isempty(strmatch('freq', dimtok)), Nfreq = length(data.freq); else Nfreq = 1; end
        else
            Nfreq = 1;
            Ntime = 1;
        end

        %convert old style source representation into new style
        if isfield(data, 'avg') && isfield(data.avg, 'mom') && (isfield(data, 'freq') || isfield(data, 'frequency')) && strcmp(sourcedimord, 'rpt_pos')
            %frequency domain source representation convert to single trial power
            Npos   = size(data.pos,1);
            Nrpt   = size(data.cumtapcnt,1);
            tmpmom = zeros(Npos, size(data.avg.mom{data.inside(1)},2));
            tmpmom(data.inside,:) = cat(1,data.avg.mom{data.inside});
            tmppow = zeros(Npos, Nrpt);
            tapcnt = [0;cumsum(data.cumtapcnt)];
            for k = 1:Nrpt
                Ntap = tapcnt(k+1)-tapcnt(k);
                tmppow(data.inside,k) = sum(abs(tmpmom(data.inside,(tapcnt(k)+1):tapcnt(k+1))).^2,2)./Ntap;
            end
            data.pow = tmppow';
            data     = rmfield(data, 'avg');
            if strcmp(inside, 'logical')
                data     = fixinside(data, 'logical');
                data.inside = repmat(data.inside(:)',[Nrpt 1]);
            end
        elseif isfield(data, 'avg') && isfield(data.avg, 'mom') && (isfield(data, 'freq') || isfield(data, 'frequency')) && strcmp(sourcedimord, 'rpttap_pos')
            %frequency domain source representation convert to single taper fourier coefficients
            Npos   = size(data.pos,1);
            Nrpt   = sum(data.cumtapcnt);
            data.fourierspctrm = complex(zeros(Nrpt, Npos), zeros(Nrpt, Npos));
            data.fourierspctrm(:, data.inside) = transpose(cat(1, data.avg.mom{data.inside}));
            data   = rmfield(data, 'avg');
        elseif isfield(data, 'avg') && isfield(data.avg, 'mom') && isfield(data, 'time') && strcmp(sourcedimord, 'pos_time')
            Npos   = size(data.pos,1);
            Nrpt   = 1;
            tmpmom = zeros(Npos, size(data.avg.mom{data.inside(1)},2));
            tmpmom(data.inside,:) = cat(1,data.avg.mom{data.inside});
            data.mom = tmpmom;
            if isfield(data.avg, 'noise')
                tmpnoise = data.avg.noise(:);
                data.noise = tmpnoise(:,ones(1,size(tmpmom,2)));
            end
            data = rmfield(data, 'avg');
            Ntime = length(data.time);
        elseif isfield(data, 'trial') && isfield(data.trial(1), 'mom') && isfield(data, 'time') && strcmp(sourcedimord, 'rpt_pos_time')
            Npos   = size(data.pos,1);
            Nrpt   = length(data.trial);
            Ntime  = length(data.time);
            tmpmom = zeros(Nrpt, Npos, Ntime);
            for k = 1:Nrpt
                tmpmom(k,data.inside,:) = cat(1,data.trial(k).mom{data.inside});
            end
            data     = rmfield(data, 'trial');
            data.mom = tmpmom;
        elseif isfield(data, 'trial') && isstruct(data.trial)
            Nrpt = length(data.trial);
        else
            Nrpt = 1;
        end

        % start with an initial specification of the dimord and dim
        if (~isfield(data, 'dim') || ~isfield(data, 'dimord'))
            if issource
                % at least it should have a Nx3 pos
                data.dim    = size(data.pos, 1);
                data.dimord = 'pos';
            elseif isvolume
                % at least it should have a 1x3 dim
                data.dim    = data.dim;
                data.dimord = 'dim1_dim2_dim3';
            end
        end

        % add the additional dimensions
        if Nfreq>1
            data.dimord = [data.dimord '_freq'];
            data.dim    = [data.dim     Nfreq];
        end
        if Ntime>1
            data.dimord = [data.dimord '_time'];
            data.dim    = [data.dim     Ntime];
        end
        if Nrpt>1 && strcmp(sourcedimord, 'rpt_pos')
            data.dimord = ['rpt_' data.dimord];
            data.dim    = [Nrpt   data.dim ];
        elseif Nrpt>1 && strcmp(sourcedimord, 'rpttap_pos')
            data.dimord = ['rpttap_' data.dimord];
            data.dim    = [Nrpt   data.dim ];
        end

        % the nested trial structure is not compatible with dimord
        if isfield(data, 'trial') && isstruct(data.trial)
            param = fieldnames(data.trial);
            for i=1:length(param)
                if isa(data.trial(1).(param{i}), 'cell')
                    concat = cell(data.dim(1), prod(data.dim(2:end)));
                else
                    concat = zeros(data.dim(1), prod(data.dim(2:end)));
                end
                for j=1:length(data.trial)
                    tmp = data.trial(j).(param{i});
                    concat(j,:) = tmp(:);
                end % for each trial
                data.trial = rmfield(data.trial, param{i});
                data.(param{i}) = reshape(concat, data.dim);
            end % for each param
            data = rmfield(data, 'trial');
        end
    end

    % ensure consistent dimensions of the source reconstructed data
    % reshape each of the source reconstructed parameters
    if issource && isfield(data, 'dim') && prod(data.dim)==size(data.pos,1)
        dim = [prod(data.dim) 1];
        %elseif issource && any(~cellfun('isempty',strfind(fieldnames(data), 'dimord')))
        %  dim = [size(data.pos,1) 1]; %sparsely represented source structure new style
    elseif isfield(data, 'dim')
        dim = [data.dim 1];
    elseif issource && ~isfield(data, 'dimord')
        dim = [size(data.pos,1) 1];
    elseif isfield(data, 'dimord')
        %HACK
        dimtok = tokenize(data.dimord, '_');
        for i=1:length(dimtok)
            if strcmp(dimtok(i), 'pos')
                dim(1,i) = size(data.pos,1);
            elseif strcmp(dimtok(i), '{pos}')
                dim(1,i) = size(data.pos,1);
            elseif ea_issubfield(data,dimtok{i})
                dim(1,i) = length(ea_getsubfield(data,dimtok{i}));
            else
                dim(1,i) = nan; % this applies to rpt, ori
            end
        end
        try
            % the following only works for rpt, not for ori
            i = find(isnan(dim));
            if ~isempty(i)
                n = fieldnames(data);
                for ii=1:length(n)
                    numels(1,ii) = numel(getfield(data,n{ii}));
                end
                nrpt = numels./prod(dim(setdiff(1:length(dim),i)));
                nrpt = nrpt(nrpt==round(nrpt));
                dim(i) = max(nrpt);
            end
        end % try
        if numel(dim)==1, dim(1,2) = 1; end;
    end

    % these fields should not be reshaped
    exclude = {'cfg' 'fwhm' 'leadfield' 'q' 'rough' 'pos'};
    if ~isempty(inside) && ~strcmp(inside, 'logical')
        % also exclude the inside/outside from being reshaped
        exclude = cat(2, exclude, {'inside' 'outside'});
    end

    param = setdiff(ea_parameterselection('all', data), exclude);
    for i=1:length(param)
        if any(param{i}=='.')
            % the parameter is nested in a substructure, which can have multiple elements (e.g. source.trial(1).pow, source.trial(2).pow, ...)
            % loop over the substructure array and reshape for every element
            tok  = tokenize(param{i}, '.');
            sub1 = tok{1};  % i.e. this would be 'trial'
            sub2 = tok{2};  % i.e. this would be 'pow'
            tmp1 = getfield(data, sub1);
            for j=1:numel(tmp1)
                tmp2 = getfield(tmp1(j), sub2);
                if prod(dim)==numel(tmp2)
                    tmp2 = reshape(tmp2, dim);
                end
                tmp1(j) = setfield(tmp1(j), sub2, tmp2);
            end
            data = setfield(data, sub1, tmp1);
        else
            tmp  = getfield(data, param{i});
            if prod(dim)==numel(tmp)
                tmp  = reshape(tmp, dim);
            end
            data = setfield(data, param{i}, tmp);
        end
    end

end

if isequal(hastrials, 'yes')
    okflag = isfield(data, 'trial');
    if ~okflag && isfield(data, 'dimord')
        % instead look in the dimord for rpt or subj
        okflag = contains(data.dimord, 'rpt') || ...
            contains(data.dimord, 'rpttap') || ...
            contains(data.dimord, 'subj');
    end
    if ~okflag
        error('This function requires data with a ''trial'' field');
    end % if okflag
end

if isequal(hasdim, 'yes') && ~isfield(data, 'dim')
    data.dim = pos2dim(data.pos);
elseif isequal(hasdim, 'no') && isfield(data, 'dim')
    data = rmfield(data, 'dim');
end % if hasdim

if isequal(hascumtapcnt, 'yes') && ~isfield(data, 'cumtapcnt')
    error('This function requires data with a ''cumtapcnt'' field');
elseif isequal(hascumtapcnt, 'no') && isfield(data, 'cumtapcnt')
    data = rmfield(data, 'cumtapcnt');
end % if hascumtapcnt

if isequal(hasdof, 'yes') && ~isfield(data, 'hasdof')
    error('This function requires data with a ''dof'' field');
elseif isequal(hasdof, 'no') && isfield(data, 'hasdof')
    data = rmfield(data, 'cumtapcnt');
end % if hasdof

if ~isempty(cmbrepresentation)
    if istimelock
        data = ea_fixcov(data, cmbrepresentation);
    elseif isfreq
        data = ea_fixcsd(data, cmbrepresentation, channelcmb);
    elseif isfreqmvar
        data = ea_fixcsd(data, cmbrepresentation, channelcmb);
    else
        error('This function requires data with a covariance, coherence or cross-spectrum');
    end
end % cmbrepresentation

if isfield(data, 'grad')
    % ensure that the gradiometer structure is up to date
    data.grad = ft_datatype_sens(data.grad);
end

if isfield(data, 'elec')
    % ensure that the electrode structure is up to date
    data.elec = ft_datatype_sens(data.elec);
end

function [select] = ea_parameterselection(param, data)

% PARAMETERSELECTION selects the parameters that are present as a volume in the data
% add that have a dimension that is compatible with the specified dimensions of the
% volume, i.e. either as a vector or as a 3D volume.
%
% Use as
%   [select] = parameterselection(param, data)
% where
%   param    cell-array, or single string, can be 'all'
%   data     structure with anatomical or functional data
%   select   returns the selected parameters as a cell-array

% Copyright (C) 2005-2008, Robert oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: parameterselection.m 9766 2014-08-06 09:51:26Z eelspa $

if ischar(param)
    param = {param};   % it should be a cell-array
elseif isempty(param)
    param = {};        % even being empty, it should be a cell-array
end

sel = find(strcmp(param, 'all'));
if ~isempty(sel)
    % the old default was a list of all known volume parameters
    % the new default is to try all fields present in the data
    allparam = fieldnames(data);
    % fields can be nested in source.avg
    if isfield(data, 'avg') && isstruct(data.avg)
        tmp = fieldnames(data.avg);
        for i=1:length(tmp)
            tmp{i} = ['avg.' tmp{i}];
        end
        allparam = cat(1, allparam, tmp);
    end
    % fields can be nested in source.trial
    if isfield(data, 'trial') && isstruct(data.trial)
        tmp = fieldnames(data.trial);
        for i=1:length(tmp)
            tmp{i} = ['trial.' tmp{i}];
        end
        allparam = cat(1, allparam, tmp);
    end
    param(sel) = [];                          % remove the 'all'
    param      = [param(:)' allparam(:)'];    % add the list of all possible parameters, these will be tested later
else
    % check all specified parameters and give support for some parameters like 'pow' and 'coh'
    % which most often will indicate 'avg.pow' and 'avg.coh'
    for i=1:length(param)
        if ~ea_issubfield(data, param{i}) && ea_issubfield(data, ['avg.' param{i}])
            % replace the parameter xxx by avg.xxx
            param{i} = ['avg.' param{i}];
        end
    end
end

% remove empty fields
param(cellfun('isempty', param)) = [];

% ensure that there are no double entries
param = unique(param);

select = {};
for i=1:length(param)
    if ea_issubfield(data, param{i})
        % the field is present, check whether the dimension is correct
        dim = size(ea_getsubfield(data, param{i}));
        if isfield(data, 'dim') && isequal(dim(:), data.dim(:))
            select{end+1} = param{i};
        elseif isfield(data, 'dim') && prod(dim)==prod(data.dim)
            select{end+1} = param{i};
        elseif isfield(data, 'dim') && numel(dim)==3 && isequal(dim(1:3)', data.dim(:))
            select{end+1} = param{i};
        elseif isfield(data, 'pos') && (prod(dim)==size(data.pos, 1) || dim(1)==size(data.pos,1))
            select{end+1} = param{i};
        elseif isfield(data, 'dimord') && (isfield(data, 'pos') || isfield(data, 'transform'))
            dimtok = tokenize(data.dimord, '_');
            nels   = 1;
            for k=1:numel(dimtok)
                if strcmp(dimtok{k}, 'rpt') || strcmp(dimtok{k}, 'rpttap')
                    nels = nels*dim(k);
                elseif strcmp(dimtok{k}, 'pos') && isfield(data, 'pos')
                    nels = nels*size(data.pos,1);
                elseif strcmp(dimtok{k}, '{pos}') && isfield(data, 'pos')
                    nels = nels*size(data.pos,1);
                elseif isfield(data, dimtok{k})
                    nels = nels*numel(getfield(data, dimtok{k}));
                end
            end
            if nels==prod(dim)
                select{end+1} = param{i};
            end
        end
    end
end

function [s] = ea_getsubfield(s, f)

% GETSUBFIELD returns a field from a structure just like the standard
% GETFIELD function, except that you can also specify nested fields
% using a '.' in the fieldname. The nesting can be arbitrary deep.
%
% Use as
%   f = getsubfield(s, 'fieldname')
% or as
%   f = getsubfield(s, 'fieldname.subfieldname')
%
% See also GETFIELD, ISSUBFIELD, SETSUBFIELD

% Copyright (C) 2005-2013, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: getsubfield.m 10198 2015-02-11 09:36:13Z roboos $

if iscell(f)
    f = f{1};
end

if ~ischar(f)
    error('incorrect input argument for fieldname');
end

t = textscan(f,'%s','delimiter','.');
t = t{1};
for k = 1:numel(t)
    s = s.(t{k});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% represent the covariance matrix in a particular manner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = ea_fixcov(data, desired)
if any(isfield(data, {'cov', 'corr'}))
    if ~isfield(data, 'labelcmb')
        current = 'full';
    else
        current = 'sparse';
    end
else
    error('Could not determine the current representation of the covariance matrix');
end
if isequal(current, desired)
    % nothing to do
elseif strcmp(current, 'full') && strcmp(desired, 'sparse')
    % FIXME should be implemented
    error('not yet implemented');
elseif strcmp(current, 'sparse') && strcmp(desired, 'full')
    % FIXME should be implemented
    error('not yet implemented');
end

function y = ea_istrue(x)

% ISTRUE ensures that a true/false input argument like "yes", "true"
% or "on" is converted into a boolean

% Copyright (C) 2009-2012, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: istrue.m 7123 2012-12-06 21:21:38Z roboos $

true_list  = {'yes' 'true' 'on' 'y' };
false_list = {'no' 'false' 'off' 'n' 'none'};

if ischar(x)
    % convert string to boolean value
    if any(strcmpi(x, true_list))
        y = true;
    elseif any(strcmpi(x, false_list))
        y = false;
    else
        error('cannot determine whether "%s" should be interpreted as true or false', x);
    end
else
    % convert numerical value to boolean
    y = logical(x);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% represent the cross-spectral density matrix in a particular manner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = ea_fixcsd(data, desired, channelcmb)

% FIXCSD converts univariate frequency domain data (fourierspctrm) into a bivariate
% representation (crsspctrm), or changes the representation of bivariate frequency
% domain data (sparse/full/sparsewithpow, sparsewithpow only works for crsspctrm or
% fourierspctrm)

% Copyright (C) 2010, Jan-Mathijs Schoffelen, Robert Oostenveld

if isfield(data, 'crsspctrm') && isfield(data, 'powspctrm')
    current = 'sparsewithpow';
elseif isfield(data, 'powspctrm')
    current = 'sparsewithpow';
elseif isfield(data, 'fourierspctrm') && ~isfield(data, 'labelcmb')
    current = 'fourier';
elseif ~isfield(data, 'labelcmb')
    current = 'full';
elseif isfield(data, 'labelcmb')
    current = 'sparse';
else
    error('Could not determine the current representation of the %s matrix', param);
end

% first go from univariate fourier to the required bivariate representation
if isequal(current, desired)
    % nothing to do

elseif strcmp(current, 'fourier') && strcmp(desired, 'sparsewithpow')
    dimtok = tokenize(data.dimord, '_');
    if ~isempty(strmatch('rpttap',   dimtok))
        nrpt = size(data.cumtapcnt,1);
        flag = 0;
    else
        nrpt = 1;
    end
    if ~isempty(strmatch('freq',  dimtok)), nfrq=length(data.freq);      else nfrq = 1; end
    if ~isempty(strmatch('time',  dimtok)), ntim=length(data.time);      else ntim = 1; end

    fastflag = all(data.cumtapcnt(:)==data.cumtapcnt(1));
    flag     = nrpt==1; % needed to truncate the singleton dimension upfront

    %create auto-spectra
    nchan     = length(data.label);
    if fastflag
        % all trials have the same amount of tapers
        powspctrm = zeros(nrpt,nchan,nfrq,ntim);
        ntap      = data.cumtapcnt(1);
        for p = 1:ntap
            powspctrm = powspctrm + abs(data.fourierspctrm(p:ntap:end,:,:,:,:)).^2;
        end
        powspctrm = powspctrm./ntap;
    else
        % different amount of tapers
        powspctrm = zeros(nrpt,nchan,nfrq,ntim)+i.*zeros(nrpt,nchan,nfrq,ntim);
        sumtapcnt = [0;cumsum(data.cumtapcnt(:))];
        for p = 1:nrpt
            indx   = (sumtapcnt(p)+1):sumtapcnt(p+1);
            tmpdat = data.fourierspctrm(indx,:,:,:);
            powspctrm(p,:,:,:) = (sum(tmpdat.*conj(tmpdat),1))./data.cumtapcnt(p);
        end
    end

    %create cross-spectra
    if ~isempty(channelcmb)
        ncmb      = size(channelcmb,1);
        cmbindx   = zeros(ncmb,2);
        labelcmb  = cell(ncmb,2);
        for k = 1:ncmb
            ch1 = find(strcmp(data.label, channelcmb(k,1)));
            ch2 = find(strcmp(data.label, channelcmb(k,2)));
            if ~isempty(ch1) && ~isempty(ch2)
                cmbindx(k,:)  = [ch1 ch2];
                labelcmb(k,:) = data.label([ch1 ch2])';
            end
        end

        crsspctrm = zeros(nrpt,ncmb,nfrq,ntim)+i.*zeros(nrpt,ncmb,nfrq,ntim);
        if fastflag
            for p = 1:ntap
                tmpdat1   = data.fourierspctrm(p:ntap:end,cmbindx(:,1),:,:,:);
                tmpdat2   = data.fourierspctrm(p:ntap:end,cmbindx(:,2),:,:,:);
                crsspctrm = crsspctrm + tmpdat1.*conj(tmpdat2);
            end
            crsspctrm = crsspctrm./ntap;
        else
            for p = 1:nrpt
                indx    = (sumtapcnt(p)+1):sumtapcnt(p+1);
                tmpdat1 = data.fourierspctrm(indx,cmbindx(:,1),:,:);
                tmpdat2 = data.fourierspctrm(indx,cmbindx(:,2),:,:);
                crsspctrm(p,:,:,:) = (sum(tmpdat1.*conj(tmpdat2),1))./data.cumtapcnt(p);
            end
        end
        data.crsspctrm = crsspctrm;
        data.labelcmb  = labelcmb;
    end
    data.powspctrm = powspctrm;
    data           = rmfield(data, 'fourierspctrm');
    if ntim>1
        data.dimord = 'chan_freq_time';
    else
        data.dimord = 'chan_freq';
    end

    if nrpt>1
        data.dimord = ['rpt_',data.dimord];
    end

    if flag, siz = size(data.crsspctrm); data.crsspctrm = reshape(data.crsspctrm, [siz(2:end) 1]); end
elseif strcmp(current, 'fourier') && strcmp(desired, 'sparse')

    if isempty(channelcmb), error('no channel combinations are specified'); end
    dimtok = tokenize(data.dimord, '_');
    if ~isempty(strmatch('rpttap',   dimtok))
        nrpt = size(data.cumtapcnt,1);
        flag = 0;
    else
        nrpt = 1;
    end
    if ~isempty(strmatch('freq',  dimtok)), nfrq=length(data.freq); else nfrq = 1; end
    if ~isempty(strmatch('time',  dimtok)), ntim=length(data.time); else ntim = 1; end

    flag      = nrpt==1; % flag needed to squeeze first dimension if singleton
    ncmb      = size(channelcmb,1);
    cmbindx   = zeros(ncmb,2);
    labelcmb  = cell(ncmb,2);
    for k = 1:ncmb
        ch1 = find(strcmp(data.label, channelcmb(k,1)));
        ch2 = find(strcmp(data.label, channelcmb(k,2)));
        if ~isempty(ch1) && ~isempty(ch2)
            cmbindx(k,:)  = [ch1 ch2];
            labelcmb(k,:) = data.label([ch1 ch2])';
        end
    end

    sumtapcnt = [0;cumsum(data.cumtapcnt(:))];
    fastflag  = all(data.cumtapcnt(:)==data.cumtapcnt(1));

    if fastflag && nrpt>1
        ntap = data.cumtapcnt(1);

        % compute running sum across tapers
        siz = [size(data.fourierspctrm) 1];

        for p = 1:ntap
            indx      = p:ntap:nrpt*ntap;

            if p==1.

                tmpc = zeros(numel(indx), size(cmbindx,1), siz(3), siz(4)) + ...
                    1i.*zeros(numel(indx), size(cmbindx,1), siz(3), siz(4));
            end

            for k = 1:size(cmbindx,1)
                tmpc(:,k,:,:) = data.fourierspctrm(indx,cmbindx(k,1),:,:).*  ...
                    conj(data.fourierspctrm(indx,cmbindx(k,2),:,:));
            end

            if p==1
                crsspctrm = tmpc;
            else
                crsspctrm = tmpc + crsspctrm;
            end
        end
        crsspctrm = crsspctrm./ntap;
    else
        crsspctrm = zeros(nrpt, ncmb, nfrq, ntim);
        for p = 1:nrpt
            indx    = (sumtapcnt(p)+1):sumtapcnt(p+1);
            tmpdat1 = data.fourierspctrm(indx,cmbindx(:,1),:,:);
            tmpdat2 = data.fourierspctrm(indx,cmbindx(:,2),:,:);
            crsspctrm(p,:,:,:) = (sum(tmpdat1.*conj(tmpdat2),1))./data.cumtapcnt(p);
        end
    end
    data.crsspctrm = crsspctrm;
    data.labelcmb  = labelcmb;
    data           = rmfield(data, 'fourierspctrm');
    data           = rmfield(data, 'label');
    if ntim>1
        data.dimord = 'chancmb_freq_time';
    else
        data.dimord = 'chancmb_freq';
    end

    if nrpt>1
        data.dimord = ['rpt_',data.dimord];
    end

    if flag, siz = size(data.crsspctrm); data.crsspctrm = reshape(data.crsspctrm, [siz(2:end) 1]); end
elseif strcmp(current, 'fourier') && strcmp(desired, 'full')

    % this is how it is currently and the desired functionality of prepare_freq_matrices
    dimtok = tokenize(data.dimord, '_');
    if ~isempty(strmatch('rpttap',   dimtok))
        nrpt = size(data.cumtapcnt, 1);
        flag = 0;
    else
        nrpt = 1;
        flag = 1;
    end
    if ~isempty(strmatch('rpttap',dimtok)), nrpt=size(data.cumtapcnt, 1); else nrpt = 1; end
    if ~isempty(strmatch('freq',  dimtok)), nfrq=length(data.freq);       else nfrq = 1; end
    if ~isempty(strmatch('time',  dimtok)), ntim=length(data.time);       else ntim = 1; end
    if any(data.cumtapcnt(1,:) ~= data.cumtapcnt(1,1)), error('this only works when all frequencies have the same number of tapers'); end
    nchan     = length(data.label);
    crsspctrm = zeros(nrpt,nchan,nchan,nfrq,ntim);
    sumtapcnt = [0;cumsum(data.cumtapcnt(:,1))];
    for k = 1:ntim
        for m = 1:nfrq
            for p = 1:nrpt
                %FIXME speed this up in the case that all trials have equal number of tapers
                indx   = (sumtapcnt(p)+1):sumtapcnt(p+1);
                tmpdat = transpose(data.fourierspctrm(indx,:,m,k));
                crsspctrm(p,:,:,m,k) = (tmpdat*tmpdat')./data.cumtapcnt(p);
                clear tmpdat;
            end
        end
    end
    data.crsspctrm = crsspctrm;
    data           = rmfield(data, 'fourierspctrm');

    if ntim>1
        data.dimord = 'chan_chan_freq_time';
    else
        data.dimord = 'chan_chan_freq';
    end

    if nrpt>1
        data.dimord = ['rpt_',data.dimord];
    end

    % remove first singleton dimension
    if flag || nrpt==1, siz = size(data.crsspctrm); data.crsspctrm = reshape(data.crsspctrm, siz(2:end)); end

elseif strcmp(current, 'fourier') && strcmp(desired, 'fullfast')

    dimtok = tokenize(data.dimord, '_');
    nrpt = size(data.fourierspctrm, 1);
    nchn = numel(data.label);
    nfrq = numel(data.freq);
    if ~isempty(strmatch('time',  dimtok)), ntim=numel(data.time); else ntim = 1; end

    data.fourierspctrm = reshape(data.fourierspctrm, [nrpt nchn nfrq*ntim]);
    data.fourierspctrm(~isfinite(data.fourierspctrm)) = 0;
    crsspctrm = complex(zeros(nchn,nchn,nfrq*ntim));
    for k = 1:nfrq*ntim
        tmp = transpose(data.fourierspctrm(:,:,k));
        n   = sum(tmp~=0,2);
        crsspctrm(:,:,k) = tmp*tmp'./n(1);
    end
    data           = rmfield(data, 'fourierspctrm');
    data.crsspctrm = reshape(crsspctrm, [nchn nchn nfrq ntim]);
    if isfield(data, 'time')
        data.dimord = 'chan_chan_freq_time';
    else
        data.dimord = 'chan_chan_freq';
    end

    if isfield(data, 'trialinfo'),  data = rmfield(data, 'trialinfo'); end;
    if isfield(data, 'sampleinfo'), data = rmfield(data, 'sampleinfo'); end;
    if isfield(data, 'cumsumcnt'),  data = rmfield(data, 'cumsumcnt');  end;
    if isfield(data, 'cumtapcnt'),  data = rmfield(data, 'cumtapcnt');  end;

end % convert to the requested bivariate representation

% from one bivariate representation to another
if isequal(current, desired)
    % nothing to do

elseif (strcmp(current, 'full')       && strcmp(desired, 'fourier')) || ...
        (strcmp(current, 'sparse')        && strcmp(desired, 'fourier')) || ...
        (strcmp(current, 'sparsewithpow') && strcmp(desired, 'fourier'))
    % this is not possible
    error('converting the cross-spectrum into a Fourier representation is not possible');

elseif strcmp(current, 'full') && strcmp(desired, 'sparsewithpow')
    error('not yet implemented');

elseif strcmp(current, 'sparse') && strcmp(desired, 'sparsewithpow')
    % convert back to crsspctrm/powspctrm representation: useful for plotting functions etc
    indx     = labelcmb2indx(data.labelcmb);
    autoindx = indx(indx(:,1)==indx(:,2), 1);
    cmbindx  = setdiff([1:size(indx,1)]', autoindx);

    if startsWith(data.dimord, 'rpt')
        data.powspctrm = data.crsspctrm(:, autoindx, :, :);
        data.crsspctrm = data.crsspctrm(:, cmbindx,  :, :);
    else
        data.powspctrm = data.crsspctrm(autoindx, :, :);
        data.crsspctrm = data.crsspctrm(cmbindx,  :, :);
    end
    data.label    = data.labelcmb(autoindx,1);
    data.labelcmb = data.labelcmb(cmbindx, :);

    if isempty(cmbindx)
        data = rmfield(data, 'crsspctrm');
        data = rmfield(data, 'labelcmb');
    end

elseif strcmp(current, 'full') && strcmp(desired, 'sparse')
    dimtok = tokenize(data.dimord, '_');
    if ~isempty(strmatch('rpt',   dimtok)), nrpt=size(data.cumtapcnt,1); else nrpt = 1; end
    if ~isempty(strmatch('freq',  dimtok)), nfrq=numel(data.freq);      else nfrq = 1; end
    if ~isempty(strmatch('time',  dimtok)), ntim=numel(data.time);      else ntim = 1; end
    nchan    = length(data.label);
    ncmb     = nchan*nchan;
    labelcmb = cell(ncmb, 2);
    cmbindx  = zeros(nchan, nchan);
    k = 1;
    for j=1:nchan
        for m=1:nchan
            labelcmb{k, 1} = data.label{m};
            labelcmb{k, 2} = data.label{j};
            cmbindx(m,j)   = k;
            k = k+1;
        end
    end

    % reshape all possible fields
    fn = fieldnames(data);
    for ii=1:numel(fn)
        if numel(data.(fn{ii})) == nrpt*ncmb*nfrq*ntim;
            if nrpt>1
                data.(fn{ii}) = reshape(data.(fn{ii}), nrpt, ncmb, nfrq, ntim);
            else
                data.(fn{ii}) = reshape(data.(fn{ii}), ncmb, nfrq, ntim);
            end
        end
    end
    % remove obsolete fields
    data           = rmfield(data, 'label');
    try data      = rmfield(data, 'dof'); end
    % replace updated fields
    data.labelcmb  = labelcmb;
    if ntim>1
        data.dimord = 'chancmb_freq_time';
    else
        data.dimord = 'chancmb_freq';
    end

    if nrpt>1
        data.dimord = ['rpt_',data.dimord];
    end

elseif strcmp(current, 'sparsewithpow') && strcmp(desired, 'sparse')
    % this representation for sparse data contains autospectra as e.g. {'A' 'A'} in labelcmb
    if isfield(data, 'crsspctrm')
        dimtok         = tokenize(data.dimord, '_');
        catdim         = match_str(dimtok, {'chan' 'chancmb'});
        data.crsspctrm = cat(catdim, data.powspctrm, data.crsspctrm);
        data.labelcmb  = [data.label(:) data.label(:); data.labelcmb];
        data           = rmfield(data, 'powspctrm');
    else
        data.crsspctrm = data.powspctrm;
        data.labelcmb  = [data.label(:) data.label(:)];
        data           = rmfield(data, 'powspctrm');
    end
    data = rmfield(data, 'label');

elseif strcmp(current, 'sparse') && strcmp(desired, 'full')
    dimtok = tokenize(data.dimord, '_');
    if ~isempty(strmatch('rpt',   dimtok)), nrpt=size(data.cumtapcnt,1); else nrpt = 1; end
    if ~isempty(strmatch('freq',  dimtok)), nfrq=numel(data.freq);      else nfrq = 1; end
    if ~isempty(strmatch('time',  dimtok)), ntim=numel(data.time);      else ntim = 1; end

    if ~isfield(data, 'label')
        % ensure that the bivariate spectral factorization results can be
        % processed. FIXME this is experimental and will not work if the user
        % did something weird before
        for k = 1:numel(data.labelcmb)
            tmp = tokenize(data.labelcmb{k}, '[');
            data.labelcmb{k} = tmp{1};
        end
        data.label = unique(data.labelcmb(:));
    end

    nchan     = length(data.label);
    ncmb      = size(data.labelcmb,1);
    cmbindx   = zeros(nchan,nchan);

    for k = 1:size(data.labelcmb,1)
        ch1 = find(strcmp(data.label, data.labelcmb(k,1)));
        ch2 = find(strcmp(data.label, data.labelcmb(k,2)));
        if ~isempty(ch1) && ~isempty(ch2)
            cmbindx(ch1,ch2) = k;
        end
    end

    complete = all(cmbindx(:)~=0);

    fn = fieldnames(data);
    for ii=1:numel(fn)
        if numel(data.(fn{ii})) == nrpt*ncmb*nfrq*ntim
            if nrpt==1
                data.(fn{ii}) = reshape(data.(fn{ii}), [nrpt ncmb nfrq ntim]);
            end

            tmpall = nan(nrpt,nchan,nchan,nfrq,ntim);

            for j = 1:nrpt
                for k = 1:ntim
                    for m = 1:nfrq
                        tmpdat = nan(nchan,nchan);
                        indx   = find(cmbindx);
                        if ~complete
                            % this realizes the missing combinations to be represented as the
                            % conjugate of the corresponding combination across the diagonal
                            tmpdat(indx) = reshape(data.(fn{ii})(j,cmbindx(indx),m,k),[numel(indx) 1]);
                            tmpdat       = ctranspose(tmpdat);
                        end
                        tmpdat(indx)    = reshape(data.(fn{ii})(j,cmbindx(indx),m,k),[numel(indx) 1]);
                        tmpall(j,:,:,m,k) = tmpdat;
                    end % for m
                end % for k
            end % for j

            % replace the data in the old representation with the new representation
            if nrpt>1
                data.(fn{ii}) = tmpall;
            else
                data.(fn{ii}) = reshape(tmpall, [nchan nchan nfrq ntim]);
            end
        end % if numel
    end % for ii

    % remove obsolete fields
    try data      = rmfield(data, 'powspctrm');  end
    try data      = rmfield(data, 'labelcmb');   end
    try data      = rmfield(data, 'dof');        end

    if ntim>1
        data.dimord = 'chan_chan_freq_time';
    else
        data.dimord = 'chan_chan_freq';
    end

    if nrpt>1
        data.dimord = ['rpt_',data.dimord];
    end

elseif strcmp(current, 'sparse') && strcmp(desired, 'fullfast')
    dimtok = tokenize(data.dimord, '_');
    if ~isempty(strmatch('rpt',   dimtok)), nrpt=size(data.cumtapcnt,1); else nrpt = 1; end
    if ~isempty(strmatch('freq',  dimtok)), nfrq=numel(data.freq);      else nfrq = 1; end
    if ~isempty(strmatch('time',  dimtok)), ntim=numel(data.time);      else ntim = 1; end

    if ~isfield(data, 'label')
        data.label = unique(data.labelcmb(:));
    end

    nchan     = length(data.label);
    ncmb      = size(data.labelcmb,1);
    cmbindx   = zeros(nchan,nchan);

    for k = 1:size(data.labelcmb,1)
        ch1 = find(strcmp(data.label, data.labelcmb(k,1)));
        ch2 = find(strcmp(data.label, data.labelcmb(k,2)));
        if ~isempty(ch1) && ~isempty(ch2)
            cmbindx(ch1,ch2) = k;
        end
    end

    complete = all(cmbindx(:)~=0);

    fn = fieldnames(data);
    for ii=1:numel(fn)
        if numel(data.(fn{ii})) == nrpt*ncmb*nfrq*ntim
            if nrpt==1
                data.(fn{ii}) = reshape(data.(fn{ii}), [nrpt ncmb nfrq ntim]);
            end

            tmpall = nan(nchan,nchan,nfrq,ntim);

            for k = 1:ntim
                for m = 1:nfrq
                    tmpdat = nan(nchan,nchan);
                    indx   = find(cmbindx);
                    if ~complete
                        % this realizes the missing combinations to be represented as the
                        % conjugate of the corresponding combination across the diagonal
                        tmpdat(indx) = reshape(nanmean(data.(fn{ii})(:,cmbindx(indx),m,k)),[numel(indx) 1]);
                        tmpdat       = ctranspose(tmpdat);
                    end
                    tmpdat(indx)    = reshape(nanmean(data.(fn{ii})(:,cmbindx(indx),m,k)),[numel(indx) 1]);
                    tmpall(:,:,m,k) = tmpdat;
                end % for m
            end % for k

            % replace the data in the old representation with the new representation
            if nrpt>1
                data.(fn{ii}) = tmpall;
            else
                data.(fn{ii}) = reshape(tmpall, [nchan nchan nfrq ntim]);
            end
        end % if numel
    end % for ii

    % remove obsolete fields
    try data      = rmfield(data, 'powspctrm');  end
    try data      = rmfield(data, 'labelcmb');   end
    try data      = rmfield(data, 'dof');        end

    if ntim>1
        data.dimord = 'chan_chan_freq_time';
    else
        data.dimord = 'chan_chan_freq';
    end

elseif strcmp(current, 'sparsewithpow') && any(strcmp(desired, {'full', 'fullfast'}))
    % this is how is currently done in prepare_freq_matrices
    data = ea_ft_checkdata(data, 'cmbrepresentation', 'sparse');
    data = ea_ft_checkdata(data, 'cmbrepresentation', 'full');

end % convert from one to another bivariate representation



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function source = ea_parcellated2source(data)
if ~isfield(data, 'brainordinate')
    error('projecting parcellated data onto the full brain model geometry requires the specification of brainordinates');
end
% the main structure contains the functional data on the parcels
% the brainordinate sub-structure contains the original geometrical model
source = data.brainordinate;
data   = rmfield(data, 'brainordinate');
if isfield(data, 'cfg')
    source.cfg = data.cfg;
end

fn = fieldnames(data);
fn = setdiff(fn, {'label', 'time', 'freq', 'hdr', 'cfg', 'grad', 'elec', 'dimord', 'unit'}); % remove irrelevant fields
fn(~cellfun(@isempty, regexp(fn, 'dimord$'))) = []; % remove irrelevant (dimord) fields
sel = false(size(fn));
for i=1:numel(fn)
    try
        sel(i) = ismember(getdimord(data, fn{i}), {'chan', 'chan_time', 'chan_freq', 'chan_freq_time', 'chan_chan'});
    end
end
parameter = fn(sel);

fn = fieldnames(source);
sel = false(size(fn));
for i=1:numel(fn)
    tmp = source.(fn{i});
    sel(i) = iscell(tmp) && isequal(tmp(:), data.label(:));
end
parcelparam = fn(sel);
if numel(parcelparam)~=1
    error('cannot determine which parcellation to use');
else
    parcelparam = parcelparam{1}(1:(end-5)); % minus the 'label'
end

for i=1:numel(parameter)
    source.(parameter{i}) = unparcellate(data, source, parameter{i}, parcelparam);
end

% copy over fields (these are necessary for visualising the data in ft_sourceplot)
source = copyfields(data, source, {'time', 'freq'});



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = ea_source2volume(data)

if isfield(data, 'dimord')
    % it is a modern source description

    %this part depends on the assumption that the list of positions is describing a full 3D volume in
    %an ordered way which allows for the extraction of a transformation matrix
    %i.e. slice by slice
    try
        if isfield(data, 'dim')
            data.dim = pos2dim(data.pos, data.dim);
        else
            data.dim = pos2dim(data);
        end
    catch
    end
end

if isfield(data, 'dim') && length(data.dim)>=3
    data.transform = pos2transform(data.pos, data.dim);
end

% remove the unwanted fields
data = removefields(data, {'pos', 'xgrid', 'ygrid', 'zgrid', 'tri', 'tet', 'hex'});

% make inside a volume
data = fixinside(data, 'logical');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = ea_freq2raw(freq)

if isfield(freq, 'powspctrm')
    param = 'powspctrm';
elseif isfield(freq, 'fourierspctrm')
    param = 'fourierspctrm';
else
    error('not supported for this data representation');
end

if strcmp(freq.dimord, 'rpt_chan_freq_time') || strcmp(freq.dimord, 'rpttap_chan_freq_time')
    dat = freq.(param);
elseif strcmp(freq.dimord, 'chan_freq_time')
    dat = freq.(param);
    dat = reshape(dat, [1 size(dat)]); % add a singleton dimension
else
    error('not supported for dimord %s', freq.dimord);
end

nrpt  = size(dat,1);
nchan = size(dat,2);
nfreq = size(dat,3);
ntime = size(dat,4);
data = [];
% create the channel labels like "MLP11@12Hz""
k = 0;
for i=1:nfreq
    for j=1:nchan
        k = k+1;
        data.label{k} = sprintf('%s@%dHz', freq.label{j}, freq.freq(i));
    end
end
% reshape and copy the data as if it were timecourses only
for i=1:nrpt
    data.time{i}  = freq.time;
    data.trial{i} = reshape(dat(i,:,:,:), nchan*nfreq, ntime);
    if any(isnan(data.trial{i}(1,:)))
        tmp = data.trial{i}(1,:);
        begsmp = find(isfinite(tmp),1, 'first');
        endsmp = find(isfinite(tmp),1, 'last' );
        data.trial{i} = data.trial{i}(:, begsmp:endsmp);
        data.time{i}  = data.time{i}(begsmp:endsmp);
    end
end

if isfield(freq, 'trialinfo')
    data.trialinfo = freq.trialinfo;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = ea_chan2freq(data)
data.dimord = [data.dimord '_freq'];
data.freq   = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = ea_chan2timelock(data)
data.dimord = [data.dimord '_time'];
data.time   = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spike] = ea_raw2spike(data)
fprintf('converting raw data into spike data\n');
nTrials 	 = length(data.trial);
[spikelabel] = ea_detectspikechan(data);
spikesel     = match_str(data.label, spikelabel);
nUnits       = length(spikesel);
if nUnits==0
    error('cannot convert raw data to spike format since the raw data structure does not contain spike channels');
end

trialTimes  = zeros(nTrials,2);
for iUnit = 1:nUnits
    unitIndx = spikesel(iUnit);
    spikeTimes  = []; % we dont know how large it will be, so use concatenation inside loop
    trialInds   = [];
    for iTrial = 1:nTrials

        % read in the spike times
        [spikeTimesTrial]    = ea_getspiketimes(data, iTrial, unitIndx);
        nSpikes              = length(spikeTimesTrial);
        spikeTimes           = [spikeTimes; spikeTimesTrial(:)];
        trialInds            = [trialInds; ones(nSpikes,1)*iTrial];

        % get the begs and ends of trials
        hasNum = find(~isnan(data.time{iTrial}));
        if iUnit==1, trialTimes(iTrial,:) = data.time{iTrial}([hasNum(1) hasNum(end)]); end
    end

    spike.label{iUnit}     = data.label{unitIndx};
    spike.waveform{iUnit}  = [];
    spike.time{iUnit}      = spikeTimes(:)';
    spike.trial{iUnit}     = trialInds(:)';

    if iUnit==1, spike.trialtime             = trialTimes; end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for detection of channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spikelabel, eeglabel] = ea_detectspikechan(data)

maxRate = 2000; % default on what we still consider a neuronal signal: this firing rate should never be exceeded

% autodetect the spike channels
ntrial = length(data.trial);
nchans  = length(data.label);
spikechan = zeros(nchans,1);
for i=1:ntrial
    for j=1:nchans
        hasAllInts    = all(isnan(data.trial{i}(j,:)) | data.trial{i}(j,:) == round(data.trial{i}(j,:)));
        hasAllPosInts = all(isnan(data.trial{i}(j,:)) | data.trial{i}(j,:)>=0);
        T = nansum(diff(data.time{i}),2); % total time
        fr            = nansum(data.trial{i}(j,:),2) ./ T;
        spikechan(j)  = spikechan(j) + double(hasAllInts & hasAllPosInts & fr<=maxRate);
    end
end
spikechan = (spikechan==ntrial);

spikelabel = data.label(spikechan);
eeglabel   = data.label(~spikechan);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spikeTimes] = ea_getspiketimes(data, trial, unit)
spikeIndx       = logical(data.trial{trial}(unit,:));
spikeCount      = data.trial{trial}(unit,spikeIndx);
spikeTimes      = data.time{trial}(spikeIndx);
if isempty(spikeTimes), return; end
multiSpikes     = find(spikeCount>1);
% get the additional samples and spike times, we need only loop through the bins
[addSamples, addTimes]   = deal([]);
for iBin = multiSpikes(:)' % looping over row vector
    addTimes     = [addTimes ones(1,spikeCount(iBin))*spikeTimes(iBin)];
    addSamples   = [addSamples ones(1,spikeCount(iBin))*spikeIndx(iBin)];
end
% before adding these times, first remove the old ones
spikeTimes(multiSpikes) = [];
spikeTimes              = sort([spikeTimes(:); addTimes(:)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = ea_spike2raw(spike, fsample)

if nargin<2 || isempty(fsample)
    timeDiff = abs(diff(sort([spike.time{:}])));
    fsample  = 1/min(timeDiff(timeDiff>0));
    warning('Desired sampling rate for spike data not specified, automatically resampled to %f', fsample);
end

% get some sizes
nUnits  = length(spike.label);
nTrials = size(spike.trialtime,1);

% preallocate
data.trial(1:nTrials) = {[]};
data.time(1:nTrials)  = {[]};
for iTrial = 1:nTrials

    % make bins: note that the spike.time is already within spike.trialtime
    x = [spike.trialtime(iTrial,1):(1/fsample):spike.trialtime(iTrial,2)];
    timeBins   = [x x(end)+1/fsample] - (0.5/fsample);
    time       = (spike.trialtime(iTrial,1):(1/fsample):spike.trialtime(iTrial,2));

    % convert to continuous
    trialData = zeros(nUnits,length(time));
    for iUnit = 1:nUnits

        % get the timestamps and only select those timestamps that are in the trial
        ts       = spike.time{iUnit};
        hasTrial = spike.trial{iUnit}==iTrial;
        ts       = ts(hasTrial);

        N = histc(ts,timeBins);
        if isempty(N)
            N = zeros(1,length(timeBins)-1);
        else
            N(end) = [];
        end

        % store it in a matrix
        trialData(iUnit,:) = N;
    end

    data.trial{iTrial} = trialData;
    data.time{iTrial}  = time;

end % for all trials

% create the associated labels and other aspects of data such as the header
data.label = spike.label;
data.fsample = fsample;
if isfield(spike,'hdr'), data.hdr = spike.hdr; end
if isfield(spike,'cfg'), data.cfg = spike.cfg; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = ea_source2raw(source)

fn = fieldnames(source);
fn = setdiff(fn, {'pos', 'dim', 'transform', 'time', 'freq', 'cfg'});
for i=1:length(fn)
    dimord{i} = getdimord(source, fn{i});
end
sel = strcmp(dimord, 'pos_time');
assert(sum(sel)>0, 'the source structure does not contain a suitable field to represent as raw channel-level data');
assert(sum(sel)<2, 'the source structure contains multiple fields that can be represented as raw channel-level data');
fn     = fn{sel};
dimord = dimord{sel};

switch dimord
    case 'pos_time'
        % add fake raw channel data to the original data structure
        data.trial{1} = source.(fn);
        data.time{1}  = source.time;
        % add fake channel labels
        data.label = {};
        for i=1:size(source.pos,1)
            data.label{i} = sprintf('source%d', i);
        end
        data.label = data.label(:);
        data.cfg = source.cfg;
    otherwise
        % FIXME other formats could be implemented as well
end


% function shifting the nodes

function segmentation = ea_ft_datatype_segmentation(segmentation, varargin)

% ea_ft_datatype_segmentation describes the FieldTrip MATLAB structure for segmented
% voxel-based data and atlases. A segmentation can either be indexed or probabilistic
% (see below).
%
% A segmentation is a volumetric description which is usually derived from an anatomical
% MRI, which describes for each voxel the tissue type. It for example distinguishes
% between white matter, grey matter, csf, skull and skin. It is mainly used for masking
% in visualization, construction of volume conduction models and for construction of
% cortical sheets. An volume-based atlas is basically a very detailed segmentation with
% an anatomical label for each voxel.
%
% For example, the AFNI TTatlas+tlrc segmented brain atlas (which can be created
% with FT_READ_ATLAS) looks like this
%
%              dim: [161 191 141]        the size of the 3D volume in voxels
%        transform: [4x4 double]         affine transformation matrix for mapping the voxel coordinates to head coordinate system
%         coordsys: 'tal'                the transformation matrix maps the voxels into this (head) coordinate system
%             unit: 'mm'                 the units in which the coordinate system is expressed
%           brick0: [161x191x141 uint8]  integer values from 1 to N, the value 0 means unknown
%           brick1: [161x191x141 uint8]  integer values from 1 to M, the value 0 means unknown
%      brick0label: {Nx1 cell}
%      brick1label: {Mx1 cell}
%
% An example of a whole-brain anatomical MRI that was segmented using FT_VOLUMESEGMENT
% looks like this
%
%         dim: [256 256 256]         the size of the 3D volume in voxels
%   transform: [4x4 double]          affine transformation matrix for mapping the voxel coordinates to head coordinate system
%    coordsys: 'ctf'                 the transformation matrix maps the voxels into this (head) coordinate system
%        unit: 'mm'                  the units in which the coordinate system is expressed
%        gray: [256x256x256 double]  probabilistic map of the gray matter
%       white: [256x256x256 double]  probabilistic map of the white matter
%         csf: [256x256x256 double]  probabilistic map of the cerebrospinal fluid
%
% An example segmentation with binary values that can be used for construction of a
% BEM volume conduction model of the head looks like this
%
%           dim: [256 256 256]         the dimensionality of the 3D volume
%     transform: [4x4 double]          affine transformation matrix for mapping the voxel coordinates to head coordinate system
%      coordsys: 'ctf'                 the transformation matrix maps the voxels into this (head) coordinate system
%          unit: 'mm'                  the units in which the coordinate system is expressed
%         brain: [256x256x256 logical] binary map representing the voxels which belong to the brain
%         scalp: [256x256x256 logical] binary map representing the voxels which belong to the scalp
%         skull: [256x256x256 logical] binary map representing the voxels which belong to the skull
%
% The examples above demonstrate that a segmentation can be indexed, i.e. consisting of
% subsequent integer numbers (1, 2, ...) or probabilistic, consisting of real numbers
% ranging from 0 to 1 that represent probabilities between 0% and 100%. An extreme case
% is one where the probability is either 0 or 1, in which case the probability can be
% represented as a binary or logical array.
%
% The only difference to the volume data representation is that the segmentation
% structure contains the additional fields xxx and xxxlabel. See FT_DATATYPE_VOLUME for
% further details.
%
% Required fields:
%   - dim, transform
%
% Optional fields:
%   - coordsys, unit
%
% Deprecated fields:
%   - none
%
% Obsoleted fields:
%   - none
%
% Revision history:
% (2012/latest) The explicit distunction between the indexed and probabilistic
% representation was made. For the indexed representation the additional
% xxxlabel cell-array was introduced.
%
% (2005) The initial version was defined.
%
% See also ea_ft_datatype, FT_DATATYPE_VOLUME, FT_DATATYPE_PARCELLATION

% Copyright (C) 2012, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ea_ft_datatype_segmentation.m 10212 2015-02-11 16:47:05Z roboos $

% get the optional input arguments, which should be specified as key-value pairs
version           = ea_ft_getopt(varargin, 'version', 'latest');
segmentationstyle = ea_ft_getopt(varargin, 'segmentationstyle');  % can be indexed or probabilistic
hasbrain          = ea_ft_getopt(varargin, 'hasbrain', 'no');     % no means that it is not required, if present it won't be removed

% convert from string into boolean
hasbrain = 0;

if strcmp(version, 'latest')
    segversion = '2012';
    volversion = 'latest';
    clear version
else
    segversion = version;
    volversion = version;
    clear version
end

if isempty(segmentation)
    return;
end

switch segversion
    case '2012'
        % determine whether the style of the input fields is probabilistic or indexed
        fn = fieldnames(segmentation);
        fn = setdiff(fn, 'inside'); % exclude the inside field from any conversions
        [indexed, probabilistic] = ea_determine_segmentationstyle(segmentation, fn, segmentation.dim);

        % ignore the fields that do not contain a segmentation
        sel = indexed | probabilistic;
        fn            = fn(sel);
        indexed       = indexed(sel);
        probabilistic = probabilistic(sel);

        % convert from an exclusive to cumulative representation
        % this is only only for demonstration purposes
        % for i=1:length(sel)
        %   segmentation.(fn{sel(i)}) = volumefillholes(segmentation.(fn{sel(i)}));
        % end

        [dum, i] = intersect(fn, {'scalp', 'skull', 'brain'});
        if numel(i)==3
            % put them in the preferred order
            fn(i) = {'scalp', 'skull', 'brain'};
        end
        [dum, i] = intersect(fn, {'skin', 'skull', 'brain'}); % this is not likely
        if numel(i)==3
            % put them in the preferred order
            fn(i) = {'skin', 'skull', 'brain'};
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ensure that the segmentation is internally consistent
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if any(probabilistic)
            segmentation = ea_fixsegmentation(segmentation, fn(probabilistic), 'probabilistic');
        end
        if any(indexed)
            segmentation = ea_fixsegmentation(segmentation, fn(indexed), 'indexed');
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % convert the segmentation to the desired style
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if isempty(segmentationstyle)
            % keep it as it is
        elseif strcmp(segmentationstyle, 'indexed') && any(probabilistic)
            segmentation  = convert_segmentationstyle(segmentation, fn(probabilistic), segmentation.dim, 'indexed');
            indexed(probabilistic)       = true;  % these are now indexed
            probabilistic(probabilistic) = false; % these are now indexed
        elseif strcmp(segmentationstyle, 'probabilistic') && any(indexed)
            segmentation  = convert_segmentationstyle(segmentation, fn(indexed), segmentation.dim, 'probabilistic');
            probabilistic(indexed) = true;  % these are now probabilistic
            indexed(indexed)       = false; % these are now probabilistic
        end % converting between probabilistic and indexed

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % add the brain if requested
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if hasbrain
            if all(indexed)
                fn = fieldnames(segmentation);
                sel = false(size(fn));
                for i=1:numel(fn)
                    sel(i) = any(strcmp(fn, [fn{i} 'label']));
                end
                fn = fn(sel);

                if numel(fn)>1
                    error('cannot construct a brain mask on the fly; this requires a single indexed representation');
                else
                    seg      = segmentation.(fn{1});
                    seglabel = segmentation.([fn{1} 'label']);
                    if ~any(strcmp(seglabel, 'brain'))
                        threshold = 0.5;
                        smooth    = 5;
                        % ensure that the segmentation contains the brain mask, if not then construct it from gray+white+csf
                        if length(intersect(seglabel, {'gray' 'white' 'csf'}))~=3
                            error('cannot construct a brain mask on the fly; this requires gray, white and csf');
                        end
                        gray  = seg==find(strcmp(seglabel, 'gray'));
                        white = seg==find(strcmp(seglabel, 'white'));
                        csf   = seg==find(strcmp(seglabel, 'csf'));
                        brain = gray + white + csf;
                        clear gray white csf seg
                        brain = volumesmooth(brain,    smooth,    'brain');
                        brain = volumethreshold(brain, threshold, 'brain');
                        % store it in the output
                        segmentation.brain = brain;
                    end % try to construct the brain
                end

            elseif all(probabilistic)
                if ~isfield(segmentation, 'brain')
                    if ~all(isfield(segmentation, {'gray' 'white' 'csf'}))
                        error('cannot construct a brain mask on the fly; this requires gray, white and csf');
                    end
                    threshold = 0.5;
                    smooth    = 5;
                    % ensure that the segmentation contains the brain mask, if not then construct it from gray+white+csf tissue probability maps
                    gray  = segmentation.gray;
                    white = segmentation.white;
                    csf   = segmentation.csf;
                    brain = gray + white + csf;
                    clear gray white csf
                    brain = volumesmooth(brain,    smooth,    'brain');
                    brain = volumethreshold(brain, threshold, 'brain');
                    % store it in the output
                    segmentation.brain = brain;
                end
            else
                error('cannot construct a brain mask on the fly; this requires a uniquely indexed or a uniquely probabilitic representation');
            end
        end % if hasbrain

    case '2005'
        % the only difference is that the indexed representation for xxx did not have the xxxlabel field prior to the 2012 version
        fn = fieldnames(segmentation);
        sel = ~cellfun(@isempty, regexp(fn, 'label$'));
        segmentation = rmfield(segmentation, fn(sel));
        % furthermore it corresponds to the oldest version of the volume representation
        volversion = '2003';

    otherwise
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        error('unsupported version "%s" for segmentation datatype', segversion);
end

% the segmentation is a speciat type of volume structure, so ensure that it also fulfills the requirements for that
segmentation = ea_ft_datatype_volume(segmentation, 'version', volversion);


function volume = ea_ft_datatype_volume(volume, varargin)

% FT_DATATYPE_VOLUME describes the FieldTrip MATLAB structure for volumetric data.
%
% The volume data structure represents data on a regular volumetric
% 3-D grid, like an anatomical MRI, a functional MRI, etc. It can
% also represent a source reconstructed estimate of the activity
% measured with MEG. In this case the source reconstruction is estimated
% or interpolated on the regular 3-D dipole grid (like a box).
%
% An example volume structure is
%       anatomy: [181x217x181 double]  the numeric data, in this case anatomical information
%           dim: [181 217 181]         the dimensionality of the 3D volume
%     transform: [4x4 double]          affine transformation matrix for mapping the voxel coordinates to the head coordinate system
%          unit: 'mm'                  geometrical units of the coordinate system
%      coordsys: 'ctf'                 description of the coordinate system
%
% Required fields:
%   - transform, dim
%
% Optional fields:
%   - anatomy, prob, stat, grey, white, csf, or any other field with dimensions that are consistent with dim
%   - size, coordsys
%
% Deprecated fields:
%   - dimord
%
% Obsoleted fields:
%   - none
%
% Revision history:
%
% (2014) The subfields in the avg and trial fields are now present in the
% main structure, e.g. source.avg.pow is now source.pow. Furthermore, the
% inside is always represented as logical array.
%
% (2012b) Ensure that the anatomy-field (if present) does not contain
% infinite values.
%
% (2012) A placeholder 2012 version was created that ensured the axes
% of the coordinate system to be right-handed. This actually never
% has made it to the default version. An executive decision regarding
% this has not been made as far as I (JM) am aware, and probably it's
% a more principled approach to keep the handedness free, so don't mess
% with it here. However, keep this snippet of code for reference.
%
% (2011) The dimord field was deprecated and we agreed that volume
% data should be 3-dimensional and not N-dimensional with arbitary
% dimensions. In case time-frequency recolved data has to be represented
% on a 3-d grid, the source representation should be used.
%
% (2010) The dimord field was added by some functions, but not by all
%
% (2003) The initial version was defined
%
% See also ea_ft_datatype, FT_DATATYPE_COMP, FT_DATATYPE_DIP, FT_DATATYPE_FREQ,
% FT_DATATYPE_MVAR, FT_DATATYPE_RAW, FT_DATATYPE_SOURCE, FT_DATATYPE_SPIKE,
% FT_DATATYPE_TIMELOCK, FT_DATATYPE_VOLUME

% Copyright (C) 2011-2015, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ft_datatype_volume.m 10174 2015-02-06 14:41:21Z roboos $

% get the optional input arguments, which should be specified as key-value pairs
version = ea_ft_getopt(varargin, 'version', 'latest');

if strcmp(version, 'latest')
    version = '2014';
end

if isempty(volume)
    return;
end

% it should  never have contained these, but they might be present due to an unclear
% distinction between the volume and the source representation
if isfield(volume, 'xgrid'),     volume = rmfield(volume, 'xgrid');     end
if isfield(volume, 'ygrid'),     volume = rmfield(volume, 'ygrid');     end
if isfield(volume, 'zgrid'),     volume = rmfield(volume, 'zgrid');     end
if isfield(volume, 'frequency'), volume = rmfield(volume, 'frequency'); end
if isfield(volume, 'latency'),   volume = rmfield(volume, 'latency');   end

if isfield(volume, 'pos')
    if ~isfield(volume, 'dim')
        volume.dim = pos2dim(volume.pos);
    end
    assert(prod(volume.dim)==size(volume.pos,1), 'dimensions are inconsistent with number of grid positions');
    if  ~isfield(volume, 'transform')
        volume.transform = pos2transform(volume.pos, volume.dim);
    end
    volume = rmfield(volume, 'pos');
end

switch version
    case '2014'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isfield(volume, 'dimord')
            volume = rmfield(volume, 'dimord');
        end

        if isfield(volume, 'anatomy')
            volume.anatomy(~isfinite(volume.anatomy)) = 0;
        end

        if isfield(volume, 'avg') && isstruct(volume.avg)
            % move the average fields to the main structure
            fn = fieldnames(volume.avg);
            for i=1:length(fn)
                dat = volume.avg.(fn{i});
                try
                    volume.(fn{i}) = reshape(dat, volume.dim);
                catch
                    warning('could not reshape %s to expected dimensions');
                    volume.(fn{i}) = dat;
                end
                clear dat
            end
            volume = rmfield(volume, 'avg');
        end

        % ensure that it is always logical
        volume = ea_fixinside(volume, 'logical');

    case '2012b'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isfield(volume, 'dimord')
            volume = rmfield(volume, 'dimord');
        end

        if isfield(volume, 'anatomy')
            volume.anatomy(~isfinite(volume.anatomy)) = 0;
        end

    case '2012'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % THIS ONE DOES NOT SEEM TO HAVE EVER BEEN USED
        % HOWEVER, KEEP IT FOR DOCUMENTATION PURPOSES

        if isfield(volume, 'dimord')
            volume = rmfield(volume, 'dimord');
        end

        % ensure the axes system in the transformation matrix to be
        % right-handed
        volume = volumeflip(volume, 'right');

    case '2011'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isfield(volume, 'dimord')
            volume = rmfield(volume, 'dimord');
        end

    case '2010'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % this might have been N-dimensional and contained a dimord, but in general cannot
        % be reconstructed on the fly

    case '2003'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isfield(volume, 'dimord')
            volume = rmfield(volume, 'dimord');
        end

    otherwise
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        error('unsupported version "%s" for volume datatype', version);
end

function [source] = ea_fixinside(source, opt)

% FIXINSIDE ensures that the region of interest (which is indicated by the
% field "inside") is consistently defined for source structures and volume
% structures. Furthermore, it solves backward compatibility problems.
%
% Use as
%   [source] = fixinside(source, 'logical');
% or
%   [source] = fixinside(source, 'index');

% Copyright (C) 2006, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: fixinside.m 9663 2014-06-22 07:06:19Z roboos $


if nargin<2
    opt = 'logical';
end

if ~isfield(source, 'inside')
    if isfield(source, 'pos')
        % assume that all positions are inside the region of interest
        source.inside  = [1:size(source.pos,1)]';
        source.outside = [];
    elseif isfield(source, 'dim')
        source.inside  = [1:prod(source.dim)]';
        source.outside = [];
    end
end

if ~isfield(source, 'inside')
    % nothing to do
    return;
end

% determine the format
if isa(source.inside, 'logical')
    logicalfmt = 1;
elseif all(source.inside(:)==0 | source.inside(:)==1)
    source.inside = logical(source.inside);
    logicalfmt = 1;
else
    logicalfmt = 0;
end

if ~logicalfmt && strcmp(opt, 'logical')
    % convert to a logical array
    if ~isfield(source, 'outside')
        source.outside = [];
    end
    inside(source.inside)  = (1==1);  % true
    inside(source.outside) = (1==0);  % false
    source.inside = inside(:);
    if isfield(source, 'outside')
        source = rmfield(source, 'outside');
    end
elseif logicalfmt && strcmp(opt, 'index')
    % convert to a vectors with indices
    tmp = source.inside;
    source.inside  = find( tmp(:));
    source.outside = find(~tmp(:));
else
    % nothing to do
end

function segmentation = ea_fixsegmentation(segmentation, fn, style)

% FIXSEGMENTATION is a helper function that ensures the segmentation to be internally
% consistent. It is used by ea_ft_datatype_segmentation and ft_datatype_parcellation.
%
% % See also CONVERT_SEGMENTATIONSTYLE, DETERMINE_SEGMENTATIONSTYLE

switch style
    case 'indexed'

        for i=1:length(fn)
            indexval = unique(segmentation.(fn{i})(:));  % find the unique tissue types
            indexval = indexval(indexval~=0);            % these are the only ones that matter

            if any(indexval<0)
                error('an indexed representation cannot contain negative numbers');
            end

            if ~isfield(segmentation, [fn{i} 'label'])
                % ensure that the tissues have labels
                indexlabel = {};
                for j=1:length(indexval)
                    indexlabel{indexval(j)} = sprintf('tissue %d', indexval(j));
                end
                segmentation.([fn{i} 'label']) = indexlabel;
            else
                % ensure that the tissue labels are consistent with the index values
                indexlabel = segmentation.([fn{i} 'label']);
                if numel(indexval)>numel(indexlabel)
                    error('each index value should have a corresponding entry in %s', [fn{i} 'label']);
                elseif any(cellfun(@isempty, indexlabel(indexval)))
                    error('each index value should have a corresponding entry in %s', [fn{i} 'label']);
                end
                % the checks above allow for the situation where
                %   indexval   = [1 2 4]
                %   indexlabel = {'a', 'b', 'c', 'd'} or {'a', 'b', [], 'd'}
                % which happens if the segmentation unexpectedly does not contain a certain tissue type
            end

            % ensure that the indices are subsequent integers, i.e. [1 2 3] rather than [1 2 4]
            for j=1:length(indexval)
                tmp = segmentation.(fn{i});
                tmp(tmp==indexval(j)) = j;
                segmentation.(fn{i}) = tmp;
            end
            segmentation.([fn{i} 'label']) = segmentation.([fn{i} 'label'])(indexval);
        end
        clear tmp indexval indexlabel

    case 'probabilistic'

        % convert from a cumulative to an exclusive representation
        within = false(length(fn));
        if length(fn)>4
            % test for each tissue whether it is overlapping with or contained in each other tissue
            warning('more than 4 tissue types, this may take a while');
        end
        for i=1:length(fn)
            segi = segmentation.(fn{i})>0;
            for j=1:length(fn)
                if i==j
                    % don't test for self-overlap
                    continue
                end
                if ~any(segi(:))
                    % don't bother to test completely empty segmentations
                    continue
                end
                segj = segmentation.(fn{j})>0;
                within(i,j) = all(segj(segi(:))); % segi is fully contained in segj
                if i~=j && within(i,j)
                    fprintf('the %s is fully contained in the %s, removing it from the %s\n', fn{i}, fn{j}, fn{j});
                    segmentation.(fn{j})(segi) = 0;
                end
            end
        end
        clear segi segj within

    otherwise
        error('unsupported style "%s"', style);
end


function val = ea_ft_getopt(opt, key, default, emptymeaningful)

% ea_ft_getopt gets the value of a specified option from a configuration structure
% or from a cell-array with key-value pairs.
%
% Use as
%   val = ea_ft_getopt(s, key, default, emptymeaningful)
% where the input values are
%   s               = structure or cell-array
%   key             = string
%   default         = any valid MATLAB data type
%   emptymeaningful = boolean value (optional, default = 0)
%
% If the key is present as field in the structure, or as key-value
% pair in the cell-array, the corresponding value will be returned.
%
% If the key is not present, ea_ft_getopt will return an empty array.
%
% If the key is present but has an empty value, then the emptymeaningful
% flag specifies whether the empty value or the default value should
% be returned. If emptymeaningful==true, then an empty array will be
% returned. If emptymeaningful==false, then the specified default will
% be returned.
%
% See also FT_SETOPT, FT_CHECKOPT

% Copyright (C) 2011-2012, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ea_ft_getopt.m 7123 2012-12-06 21:21:38Z roboos $

if nargin<3
    default = [];
end

if nargin < 4
    emptymeaningful = 0;
end

if isa(opt, 'struct') || isa(opt, 'config')
    % get the key-value from the structure
    fn = fieldnames(opt);
    if ~any(strcmp(key, fn))
        val = default;
    else
        val = opt.(key);
    end

elseif isa(opt, 'cell')
    % get the key-value from the cell-array
    if mod(length(opt),2)
        error('optional input arguments should come in key-value pairs, i.e. there should be an even number');
    end

    % the 1st, 3rd, etc. contain the keys, the 2nd, 4th, etc. contain the values
    keys = opt(1:2:end);
    vals = opt(2:2:end);

    % the following may be faster than cellfun(@ischar, keys)
    valid = false(size(keys));
    for i=1:numel(keys)
        valid(i) = ischar(keys{i});
    end

    if ~all(valid)
        error('optional input arguments should come in key-value pairs, the optional input argument %d is invalid (should be a string)', i);
    end

    hit = find(strcmpi(key, keys));
    if isempty(hit)
        % the requested key was not found
        val = default;
    elseif length(hit)==1
        % the requested key was found
        val = vals{hit};
    else
        error('multiple input arguments with the same name');
    end

elseif isempty(opt)
    % no options are specified, return default
    val = default;
end % isstruct or iscell or isempty

if isempty(val) && ~isempty(default) && ~emptymeaningful
    % use the default value instead of the empty input that was specified:
    % this applies for example if you do functionname('key', []), where
    % the empty is meant to indicate that the user does not know or care
    % what the value is
    val = default;
end





function [output] = ea_ft_transform_geometry(transform, input)

% FT_TRANSFORM_GEOMETRY applies a homogeneous coordinate transformation to
% a structure with geometric information. These objects include:
%  - volume conductor geometry, consisting of a mesh, a set of meshes, a
%      single sphere, or multiple spheres.
%  - gradiometer of electrode structure containing sensor positions and
%      coil orientations (for MEG).
%  - headshape description containing positions in 3D space.
%  - sourcemodel description containing positions and optional orientations
%      in 3D space.
%
% The units in which the transformation matrix is expressed are assumed to
% be the same units as the units in which the geometric object is
% expressed. Depending on the input object, the homogeneous transformation
% matrix should be limited to a rigid-body translation plus rotation
% (MEG-gradiometer array), or to a rigid-body translation plus rotation
% plus a global rescaling (volume conductor geometry).
%
% Use as
%   output = ft_transform_geometry(transform, input)

% Copyright (C) 2011, Jan-Mathijs Schoffelen
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ft_transform_geometry.m$

% flg rescaling check
allowscaling = ~ea_ft_senstype(input, 'meg');

% determine the rotation matrix
rotation = eye(4);
rotation(1:3,1:3) = transform(1:3,1:3);

if any(abs(transform(4,:)-[0 0 0 1])>100*eps)
    error('invalid transformation matrix');
end

if ~allowscaling
    % allow for some numerical imprecision
    if abs(det(rotation)-1)>1e-6%100*eps
        %if abs(det(rotation)-1)>100*eps  % allow for some numerical imprecision
        error('only a rigid body transformation without rescaling is allowed');
    end
end

if allowscaling
    % FIXME build in a check for uniform rescaling probably do svd or so
    % FIXME insert check for nonuniform scaling, should give an error
end

tfields   = {'pos' 'pnt' 'o' 'chanpos' 'chanposorg' 'coilpos' 'elecpos', 'nas', 'lpa', 'rpa', 'zpoint'}; % apply rotation plus translation
rfields   = {'ori' 'nrm' 'coilori'}; % only apply rotation
mfields   = {'transform'};           % plain matrix multiplication
recfields = {'fid' 'bnd' 'orig'};    % recurse into these fields
% the field 'r' is not included here, because it applies to a volume
% conductor model, and scaling is not allowed, so r will not change.

fnames    = fieldnames(input);
for k = 1:numel(fnames)
    if ~isempty(input.(fnames{k}))
        if any(strcmp(fnames{k}, tfields))
            input.(fnames{k}) = ea_apply(transform, input.(fnames{k}));
        elseif any(strcmp(fnames{k}, rfields))
            input.(fnames{k}) = ea_apply(rotation, input.(fnames{k}));
        elseif any(strcmp(fnames{k}, mfields))
            input.(fnames{k}) = transform*input.(fnames{k});
        elseif any(strcmp(fnames{k}, recfields))
            for j = 1:numel(input.(fnames{k}))
                input.(fnames{k})(j) = ea_ft_transform_geometry(transform, input.(fnames{k})(j));
            end
        else
            % do nothing
        end
    end
end
output = input;
return;


function [new] = ea_apply(transform, old)
old(:,4) = 1;
new = old * transform';
new = new(:,1:3);


function [type] = ea_ft_senstype(input, desired)

% FT_SENSTYPE determines the type of acquisition device by looking at the
% channel names and comparing them with predefined lists.
%
% Use as
%   [type] = ft_senstype(sens)
%   [flag] = ft_senstype(sens, desired)
%
% The output type can be any of the following
%   'ctf64'
%   'ctf151'
%   'ctf151_planar'
%   'ctf275'
%   'ctf275_planar'
%   'bti148'
%   'bti148_planar'
%   'bti248'
%   'bti248_planar'
%   'bti248grad'
%   'bti248grad_planar'
%   'itab28'
%   'itab153'
%   'itab153_planar'
%   'yokogawa9'
%   'yokogawa64'
%   'yokogawa64_planar'
%   'yokogawa160'
%   'yokogawa160_planar'
%   'yokogawa440'
%   'yokogawa440'_planar
%   'ext1020' (this includes eeg1020, eeg1010 and eeg1005)
%   'neuromag122'
%   'neuromag306'
%   'babysquid74'
%   'egi32'
%   'egi64'
%   'egi128'
%   'egi256'
%   'biosemi64'
%   'biosemi128'
%   'biosemi256'
%   'neuralynx'
%   'plexon'
%   'artinis'
%   'eeg' (this was called 'electrode' in older versions)
%   'meg' (this was called 'magnetometer' in older versions)
%   'nirs'
%
% The optional input argument for the desired type can be any of the above,
% or any of the following
%   'eeg'
%   'meg'
%   'meg_planar'
%   'meg_axial'
%   'ctf'
%   'bti'
%   'neuromag'
%   'yokogawa'
% If you specify the desired type, this function will return a boolean
% true/false depending on the input data.
%
% Besides specifiying a sensor definition (i.e. a grad or elec structure,
% see FT_DATATYPE_SENS), it is also possible to give a data structure
% containing a grad or elec field, or giving a list of channel names (as
% cell-arrray). So assuming that you have a FieldTrip data structure, any
% of the following calls would also be fine.
%   ft_senstype(hdr)
%   ft_senstype(data)
%   ft_senstype(data.label)
%   ft_senstype(data.grad)
%   ft_senstype(data.grad.label)
%
% See also FT_SENSLABEL, FT_CHANTYPE, FT_READ_SENS, FT_COMPUTE_LEADFIELD, FT_DATATYPE_SENS

% Copyright (C) 2007-2014, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ft_senstype.m 10340 2015-04-17 14:10:04Z jorhor $

% these are for remembering the type on subsequent calls with the same input arguments
persistent previous_argin previous_argout

% this is to avoid a recursion loop
persistent recursion
if isempty(recursion)
    recursion = false;
end

if iscell(input) && numel(input)<4 && ~all(cellfun(@ischar, input))
    % this might represent combined EEG, ECoG and/or MEG
    type = cell(size(input));
    if nargin<2
        desired = cell(size(input)); % empty elements
    end
    for i=1:numel(input)
        type{i} = ea_ft_senstype(input{i}, desired{i});
    end
    return
end

if nargin<2
    % ensure that all input arguments are defined
    desired = [];
end

current_argin = {input, desired};
if isequal(current_argin, previous_argin)
    % don't do the type detection again, but return the previous output from cache
    type = previous_argout{1};
    return
end

isdata   = isa(input, 'struct')  && (isfield(input, 'hdr') || isfield(input, 'time') || isfield(input, 'freq') || isfield(input, 'grad') || isfield(input, 'elec') || isfield(input, 'opto'));
isheader = isa(input, 'struct')  && isfield(input, 'label') && isfield(input, 'Fs');
isgrad   = isa(input, 'struct')  && isfield(input, 'label') && isfield(input, 'pnt')  &&  isfield(input, 'ori'); % old style
iselec   = isa(input, 'struct')  && isfield(input, 'label') && isfield(input, 'pnt')  && ~isfield(input, 'ori'); % old style
isgrad   = (isa(input, 'struct') && isfield(input, 'label') && isfield(input, 'coilpos')) || isgrad;             % new style
iselec   = (isa(input, 'struct') && isfield(input, 'label') && isfield(input, 'elecpos')) || iselec;             % new style
isnirs   = isa(input, 'struct')  && isfield(input, 'label') && isfield(input, 'transceiver');
islabel  = isa(input, 'cell')    && ~isempty(input) && isa(input{1}, 'char');
haslabel = isa(input, 'struct')  && isfield(input, 'label');

if ~(isdata || isheader || isgrad || iselec || isnirs || islabel || haslabel) && isfield(input, 'hdr')
    input    = input.hdr;
    isheader = true;
end

if isdata
    % the input may be a data structure which then contains a grad/elec structure, a header or only the labels
    % preferably look at the data and not the header for the grad, because it might be re-balanced and/or planar
    if isfield(input, 'grad')
        sens   = input.grad;
        isgrad = true;
    elseif isfield(input, 'elec')
        sens   = input.elec;
        iselec = true;
    elseif ea_issubfield(input, 'hdr.grad')
        sens   = input.hdr.grad;
        isgrad = true;
    elseif ea_issubfield(input, 'hdr.elec')
        sens   = input.hdr.elec;
        iselec = true;
    elseif ea_issubfield(input, 'hdr.opto')
        sens   = input.hdr.opto;
        isnirs = true;
    elseif ea_issubfield(input, 'hdr.label')
        sens.label = input.hdr.label;
        islabel    = true;
    elseif isfield(input, 'label')
        sens.label = input.label;
        islabel    = true;
    end

elseif isheader
    if isfield(input, 'grad')
        sens   = input.grad;
        isgrad = true;
    elseif isfield(input, 'elec')
        sens   = input.elec;
        iselec = true;
    elseif isfield(input, 'opto')
        sens   = input.opto;
        isnirs = true;
    elseif isfield(input, 'label')
        sens.label = input.label;
        islabel    = true;
    end

elseif isgrad
    sens = input;

elseif iselec
    sens = input;

elseif isnirs
    sens = input;

elseif islabel
    sens.label = input;

elseif haslabel
    % it does not resemble anything that we had expected at this location, but it does have channel labels
    % the channel labels can be used to determine the type of sensor array
    sens.label = input.label;
    islabel    = true;

else
    sens = [];
end


if isfield(sens, 'type')
    % preferably the structure specifies its own type
    type = sens.type;

    % do not make a distinction between the neuromag data with or without space in the channel names
    if strcmp(type, 'neuromag306alt')
        type = 'neuromag306';
    elseif strcmp(type, 'neuromag122alt')
        type = 'neuromag122';
    end

elseif isfield(input, 'nChans') && input.nChans==1 && isfield(input, 'label') && ~isempty(regexp(input.label{1}, '^csc', 'once'))
    % this is a single channel header that was read from a Neuralynx file, might be fcdc_matbin or neuralynx_nsc
    type = 'neuralynx';

elseif ea_issubfield(input, 'orig.FileHeader') &&  ea_issubfield(input, 'orig.VarHeader')
    % this is a complete header that was read from a Plexon *.nex file using read_plexon_nex
    type = 'plexon';

elseif ea_issubfield(input, 'orig.stname')
    % this is a complete header that was read from an ITAB dataset
    type = 'itab';

elseif ea_issubfield(input, 'orig.sys_name')
    % this is a complete header that was read from a Yokogawa dataset
    if strcmp(input.orig.sys_name, '9ch Biomagnetometer System') || input.orig.channel_count<20
        % this is the small animal system that is installed at the UCL Ear Institute
        % see http://www.ucl.ac.uk/news/news-articles/0907/09070101
        type = 'yokogawa9';
    elseif input.orig.channel_count<160
        type = 'yokogawa64';
    elseif input.orig.channel_count<300
        type = 'yokogawa160';
    else
        % FIXME this might fail if there are many bad channels
        type = 'yokogawa440';
    end

elseif ea_issubfield(input, 'orig.FILE.Ext') && strcmp(input.orig.FILE.Ext, 'edf')
    % this is a complete header that was read from an EDF or EDF+ dataset
    type = 'eeg';

else
    % start with unknown, then try to determine the proper type by looking at the labels
    type = 'unknown';

    if isgrad && isfield(sens, 'type')
        type = sens.type;

    elseif isgrad
        % this looks like MEG
        % revert the component balancing that was previously applied
        if isfield(sens, 'balance') && strcmp(sens.balance.current, 'comp')
            sens = ea_undobalancing(sens);
        end

        % determine the type of magnetometer/gradiometer system based on the channel names alone
        % this uses a recursive call to the "islabel" section further down
        type = ea_ft_senstype(sens.label);
        %sens_temp.type = type;
        if strcmp(type, 'unknown') %|| ~ft_senstype(sens_temp,'meg')
            % although we don't know the type, we do know that it is MEG
            type = 'meg';
        end

    elseif iselec
        % this looks like EEG

        % determine the type of eeg/acquisition system based on the channel names alone
        % this uses a recursive call to the "islabel" section further down
        type = ea_ft_senstype(sens.label);
        %sens_temp.type = type;
        if strcmp(type, 'unknown') %|| ~ft_senstype(sens_temp,'eeg')
            % although we don't know the type, we do know that it is EEG
            type = 'eeg';
        end

    elseif isnirs
        % this looks like NIRS

        % determine the type of eeg/acquisition system based on the channel names alone
        % this uses a recursive call to the "islabel" section further down
        type = ft_senstype(sens.label);
        %sens_temp.type = type;
        if strcmp(type, 'unknown') %|| ~ft_senstype(sens_temp,'eeg')
            % although we don't know the type, we do know that it is EEG
            type = 'nirs';
        end

    elseif islabel
        % look only at the channel labels
        if     (mean(ismember(ea_ft_senslabel('ant128'),         sens.label)) > 0.8)
            type = 'ant128';
        elseif (mean(ismember(ea_ft_senslabel('ctf275'),         sens.label)) > 0.8)
            type = 'ctf275';
        elseif (mean(ismember(ea_ft_senslabel('ctfheadloc'),     sens.label)) > 0.8)  % look at the head localization channels
            type = 'ctf275';
        elseif (mean(ismember(ea_ft_senslabel('ctf151'),         sens.label)) > 0.8)
            type = 'ctf151';
        elseif (mean(ismember(ea_ft_senslabel('ctf64'),          sens.label)) > 0.8)
            type = 'ctf64';
        elseif (mean(ismember(ea_ft_senslabel('ctf275_planar'),  sens.label)) > 0.8)
            type = 'ctf275_planar';
        elseif (mean(ismember(ea_ft_senslabel('ctf151_planar'),  sens.label)) > 0.8)
            type = 'ctf151_planar';
        elseif (mean(ismember(ea_ft_senslabel('bti248'),         sens.label)) > 0.8) % note that it might also be a bti248grad system
            type = 'bti248';
        elseif (mean(ismember(ea_ft_senslabel('bti148'),         sens.label)) > 0.8)
            type = 'bti148';
        elseif (mean(ismember(ea_ft_senslabel('bti248_planar'),  sens.label)) > 0.8) % note that it might also be a bti248grad_planar system
            type = 'bti248_planar';
        elseif (mean(ismember(ea_ft_senslabel('bti148_planar'),  sens.label)) > 0.8)
            type = 'bti148_planar';
        elseif (mean(ismember(ea_ft_senslabel('itab28'),         sens.label)) > 0.8)
            type = 'itab28';
        elseif (mean(ismember(ea_ft_senslabel('itab153'),        sens.label)) > 0.8)
            type = 'itab153';
        elseif (mean(ismember(ea_ft_senslabel('itab153_planar'), sens.label)) > 0.8)
            type = 'itab153_planar';

            % the order is important for the different yokogawa systems, because they all share the same channel names
        elseif (mean(ismember(ea_ft_senslabel('yokogawa440'),        sens.label)) > 0.7)
            type = 'yokogawa440';
        elseif (mean(ismember(ea_ft_senslabel('yokogawa440_planar'), sens.label)) > 0.7)
            type = 'yokogawa440_planar';
        elseif (mean(ismember(ea_ft_senslabel('yokogawa160'),        sens.label)) > 0.4)
            type = 'yokogawa160';
        elseif (mean(ismember(ea_ft_senslabel('yokogawa160_planar'), sens.label)) > 0.4)
            type = 'yokogawa160_planar';
        elseif (mean(ismember(ea_ft_senslabel('yokogawa64'),         sens.label)) > 0.4)
            type = 'yokogawa64';
        elseif (mean(ismember(ea_ft_senslabel('yokogawa64_planar'),  sens.label)) > 0.4)
            type = 'yokogawa64_planar';
        elseif all(ismember(ea_ft_senslabel('yokogawa9'),            sens.label))
            type = 'yokogawa9';

        elseif any(mean(ismember(ea_ft_senslabel('neuromag306'),     sens.label)) > 0.4) % there are two possibilities for the channel labels: with and without a space
            type = 'neuromag306';
        elseif any(mean(ismember(ea_ft_senslabel('neuromag122'),     sens.label)) > 0.4) % there are two possibilities for the channel labels: with and without a space
            type = 'neuromag122';

        elseif (mean(ismember(ea_ft_senslabel('biosemi256'),         sens.label)) > 0.8)
            type = 'biosemi256';
        elseif (mean(ismember(ea_ft_senslabel('biosemi128'),         sens.label)) > 0.8)
            type = 'biosemi128';
        elseif (mean(ismember(ea_ft_senslabel('biosemi64'),          sens.label)) > 0.8)
            type = 'biosemi64';
        elseif (mean(ismember(ea_ft_senslabel('egi256'),             sens.label)) > 0.8)
            type = 'egi256';
        elseif (mean(ismember(ea_ft_senslabel('egi128'),             sens.label)) > 0.8)
            type = 'egi128';
        elseif (mean(ismember(ea_ft_senslabel('egi64'),              sens.label)) > 0.8)
            type = 'egi64';
        elseif (mean(ismember(ea_ft_senslabel('egi32'),              sens.label)) > 0.8)
            type = 'egi32';

            % the following check on the fraction of channels in the user's data rather than on the fraction of channels in the predefined set
        elseif (mean(ismember(sens.label, ea_ft_senslabel('eeg1020'))) > 0.8)
            type = 'eeg1020';
        elseif (mean(ismember(sens.label, ea_ft_senslabel('eeg1010'))) > 0.8)
            type = 'eeg1010';
        elseif (mean(ismember(sens.label, ea_ft_senslabel('eeg1005'))) > 0.8)
            type = 'eeg1005';

        elseif (sum(ismember(sens.label, ea_ft_senslabel('eeg1005'))) > 10) % Otherwise it's not even worth recognizing
            type = 'ext1020'; % this will also cover small subsets of eeg1020, eeg1010 and eeg1005
        elseif any(ismember(ea_ft_senslabel('btiref'), sens.label))
            type = 'bti'; % it might be 148 or 248 channels
        elseif any(ismember(ea_ft_senslabel('ctfref'), sens.label))
            type = 'ctf'; % it might be 151 or 275 channels

            %     elseif (mean(ismember(sens.label,    ft_senslabel('nirs'))) > 0.8)
            %       type = 'nirs';
        end
    end % look at label, ori and/or pnt
end % if isfield(sens, 'type')

if strcmp(type, 'unknown') && ~recursion
    % try whether only lowercase channel labels makes a difference
    if islabel
        recursion = true;
        type = ea_ft_senstype(lower(input));
        recursion = false;
    elseif isfield(input, 'label')
        input.label = lower(input.label);
        recursion = true;
        type = ea_ft_senstype(input);
        recursion = false;
    end
end

if strcmp(type, 'unknown') && ~recursion
    % try whether only uppercase channel labels makes a difference
    if islabel
        recursion = true;
        type = ea_ft_senstype(upper(input));
        recursion = false;
    elseif isfield(input, 'label')
        input.label = upper(input.label);
        recursion = true;
        type = ea_ft_senstype(input);
        recursion = false;
    end
end

if ~isempty(desired)
    % return a boolean flag
    switch desired
        case 'ext1020'
            type = any(strcmp(type, {'eeg1005' 'eeg1010' 'eeg1020' 'ext1020'}));
        case {'eeg' 'electrode'}
            type = any(strcmp(type, {'eeg' 'electrode' 'ant128' 'biosemi64' 'biosemi128' 'biosemi256' 'egi32' 'egi64' 'egi128' 'egi256' 'eeg1005' 'eeg1010' 'eeg1020' 'ext1020'}));
        case 'biosemi'
            type = any(strcmp(type, {'biosemi64' 'biosemi128' 'biosemi256'}));
        case 'egi'
            type = any(strcmp(type, {'egi32' 'egi64' 'egi128' 'egi256'}));
        case 'meg'
            type = any(strcmp(type, {'meg' 'magnetometer' 'ctf' 'bti' 'ctf64' 'ctf151' 'ctf275' 'ctf151_planar' 'ctf275_planar' 'neuromag122' 'neuromag306' 'bti148' 'bti148_planar' 'bti248' 'bti248_planar' 'bti248grad' 'bti248grad_planar' 'yokogawa9' 'yokogawa160' 'yokogawa160_planar' 'yokogawa64' 'yokogawa64_planar' 'yokogawa440' 'yokogawa440_planar' 'itab' 'itab28' 'itab153' 'itab153_planar'}));
        case 'ctf'
            type = any(strcmp(type, {'ctf' 'ctf64' 'ctf151' 'ctf275' 'ctf151_planar' 'ctf275_planar'}));
        case 'bti'
            type = any(strcmp(type, {'bti' 'bti148' 'bti148_planar' 'bti248' 'bti248_planar' 'bti248grad' 'bti248grad_planar'}));
        case 'neuromag'
            type = any(strcmp(type, {'neuromag122' 'neuromag306'}));
        case 'yokogawa'
            type = any(strcmp(type, {'yokogawa160' 'yokogawa160_planar' 'yokogawa64' 'yokogawa64_planar' 'yokogawa440' 'yokogawa440_planar'}));
        case 'itab'
            type = any(strcmp(type, {'itab' 'itab28' 'itab153' 'itab153_planar'}));
        case 'meg_axial'
            % note that neuromag306 is mixed planar and axial
            type = any(strcmp(type, {'neuromag306' 'ctf64' 'ctf151' 'ctf275' 'bti148' 'bti248' 'bti248grad' 'yokogawa9' 'yokogawa64' 'yokogawa160' 'yokogawa440'}));
        case 'meg_planar'
            % note that neuromag306 is mixed planar and axial
            type = any(strcmp(type, {'neuromag122' 'neuromag306' 'ctf151_planar' 'ctf275_planar' 'bti148_planar' 'bti248_planar' 'bti248grad_planar' 'yokogawa160_planar' 'yokogawa64_planar' 'yokogawa440_planar'}));
        otherwise
            type = any(strcmp(type, desired));
    end % switch desired
end % detemine the correspondence to the desired type

% remember the current input and output arguments, so that they can be
% reused on a subsequent call in case the same input argument is given
current_argout = {type};
previous_argin  = current_argin;
previous_argout = current_argout;

return % ft_senstype main()


function [r] = ea_issubfield(s, f)

% ea_issubfield tests for the presence of a field in a structure just like the standard
% Matlab ISFIELD function, except that you can also specify nested fields
% using a '.' in the fieldname. The nesting can be arbitrary deep.
%
% Use as
%   f = ea_issubfield(s, 'fieldname')
% or as
%   f = ea_issubfield(s, 'fieldname.subfieldname')
%
% This function returns true if the field is present and false if the field
% is not present.
%
% See also ISFIELD, GETSUBFIELD, SETSUBFIELD

% Copyright (C) 2005-2013, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ea_issubfield.m 10237 2015-02-16 19:53:27Z roboos $

%try
%  getsubfield(s, f);    % if this works, then the subfield must be present
%  r = true;
%catch
%  r = false;                % apparently the subfield is not present
%end

t = textscan(f,'%s','delimiter','.');
t = t{1};
r = true;
for k = 1:numel(t)
    if isfield(s, t{k})
        s = s.(t{k});
    else
        r = false;
        return;
    end
end

function label = ea_ft_senslabel(type, varargin)

% FT_SENSLABEL returns a list of predefined sensor labels given the
% EEG or MEG system type which can be used to detect the type of data.
%
% Use as
%  label = ft_senslabel(type)
%
% The input sensor array type can be any of the following
%  'ant128'
%  'biosemi64'
%  'biosemi128'
%  'biosemi256'
%  'bti148'
%  'bti148_planar'
%  'bti248'
%  'bti248_planar'
%  'btiref'
%  'ctf151'
%  'ctf151_planar'
%  'ctf275'
%  'ctf275_planar'
%  'ctfheadloc'
%  'ctfref'
%  'eeg1005'
%  'eeg1010'
%  'eeg1020'
%  'ext1020'
%  'egi32'
%  'egi64'
%  'egi128'
%  'egi256'
%  'neuromag122'
%  'neuromag306'
%  'itab28'
%  'itab153'
%  'itab153_planar'
%  'yokogawa9'
%  'yokogawa64'
%  'yokogawa64_planar'
%  'yokogawa160'
%  'yokogawa160_planar'
%  'yokogawa440'
%  'yokogawa440_planar'
%
% It is also possible to specify
%  'eeg'
%  'electrode'
% although for these an empty set of labels (i.e. {}) will be returned.
%
% See also FT_SENSTYPE, FT_CHANNELSELECTION

% Copyright (C) 2007-2013, Robert Oostenveld
% Copyright (C) 2008, Vladimir Litvak
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%  FieldTrip is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  FieldTrip is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ft_senslabel.m 10340 2015-04-17 14:10:04Z jorhor $

% these are for speeding up subsequent calls with the same input arguments
persistent eeg electrode ant128 btiref bti148 bti148_planar bti148_planar_combined bti248 bti248_planar bti248_planar_combined ctfref ctfheadloc ctf64 ctf151 ctf151_planar ctf151_planar_combined ctf275 ctf275_planar ctf275_planar_combined neuromag122 neuromag122_combined neuromag306 neuromag306_combined eeg1020 eeg1010 eeg1005 ext1020 biosemi64 biosemi128 biosemi256 egi32 egi64 egi128 egi256 itab28 itab153 itab153_planar itab153_planar_combined yokogawa9 yokogawa64 yokogawa64_planar yokogawa64_planar_combined yokogawa160 yokogawa160_planar yokogawa160_planar_combined yokogawa440 yokogawa440_planar yokogawa440_planar_combined
% these are for backward compatibility
persistent neuromag122alt neuromag122alt_combined neuromag306alt neuromag306alt_combined

if nargin<1
    % ensure that all input arguments are defined
    type = 'none';
end

% get the optional input arguments
output  = ea_ft_getopt(varargin, 'output', 'normal'); % 'normal' or 'planarcombined'

if ~exist('type', 'var')
    error('the requested sensor type "%s" is not supported', type);

elseif isempty(eval(type))
    % assign the list of channels only once, keep it as persistent variable

    switch type
        case 'ant128'
            label = {
                'Z1'
                'Z2'
                'Z3'
                'Z4'
                'Z5'
                'Z6'
                'Z7'
                'Z8'
                'Z9'
                'Z10'
                'Z11'
                'Z12'
                'Z13'
                'Z14'
                'L1'
                'L2'
                'L3'
                'L4'
                'L5'
                'L6'
                'L7'
                'L8'
                'L9'
                'L10'
                'L11'
                'L12'
                'L13'
                'L14'
                'LL1'
                'LL2'
                'LL3'
                'LL4'
                'LL5'
                'LL6'
                'LL7'
                'LL8'
                'LL9'
                'LL10'
                'LL11'
                'LL12'
                'LL13'
                'LA1'
                'LA2'
                'LA3'
                'LA4'
                'LA5'
                'LB1'
                'LB2'
                'LB3'
                'LB4'
                'LB5'
                'LB6'
                'LC1'
                'LC2'
                'LC3'
                'LC4'
                'LC5'
                'LC6'
                'LC7'
                'LD1'
                'LD2'
                'LD3'
                'LD4'
                'LD5'
                'LD6'
                'LD7'
                'LE1'
                'LE2'
                'LE3'
                'LE4'
                'Lm'
                'R1'
                'R2'
                'R3'
                'R4'
                'R5'
                'R6'
                'R7'
                'R8'
                'R9'
                'R10'
                'R11'
                'R12'
                'R13'
                'R14'
                'RR1'
                'RR2'
                'RR3'
                'RR4'
                'RR5'
                'RR6'
                'RR7'
                'RR8'
                'RR9'
                'RR10'
                'RR11'
                'RR12'
                'RR13'
                'RA1'
                'RA2'
                'RA3'
                'RA4'
                'RA5'
                'RB1'
                'RB2'
                'RB3'
                'RB4'
                'RB5'
                'RB6'
                'RC1'
                'RC2'
                'RC3'
                'RC4'
                'RC5'
                'RC6'
                'RC7'
                'RD1'
                'RD2'
                'RD3'
                'RD4'
                'RD5'
                'RD6'
                'RD7'
                'RE1'
                'RE2'
                'RE3'
                'RE4'
                'Rm'
                };

        case 'btiref'
            label = {
                'MRxA'
                'MRyA'
                'MRzA'
                'MLxA'
                'MLyA'
                'MLzA'
                'MCxA'
                'MCyA'
                'MCzA'
                'MRxaA'
                'MRyaA'
                'MRzaA'
                'MLxaA'
                'MLyaA'
                'MLzaA'
                'MCxaA'
                'MCyaA'
                'MCzaA'
                'GxxA'
                'GyxA'
                'GzxA'
                'GyyA'
                'GzyA'
                };

        case 'bti148'
            label = cell(148,1);
            for i=1:148
                label{i,1} = sprintf('A%d', i);
            end

        case 'bti148_planar'
            label = cell(148,3);
            for i=1:148
                label{i,1} = sprintf('A%d_dH', i);
                label{i,2} = sprintf('A%d_dV', i);
                label{i,3} = sprintf('A%d', i);
            end
            bti148_planar_combined = label(:,3);
            label = label(:,1:2);

        case 'bti248'
            label = cell(248,1);
            for i=1:248
                label{i,1} = sprintf('A%d', i);
            end

        case 'bti248_planar'
            label = cell(248,3);
            for i=1:248
                label{i,1} = sprintf('A%d_dH', i);
                label{i,2} = sprintf('A%d_dV', i);
                label{i,3} = sprintf('A%d', i);
            end
            bti248_planar_combined = label(:,3);
            label = label(:,1:2);

        case 'ctfref'
            label = {
                'BG1'
                'BG2'
                'BG3'
                'BP1'
                'BP2'
                'BP3'
                'BR1'
                'BR2'
                'BR3'
                'G11'
                'G12'
                'G13'
                'G22'
                'G23'
                'P11'
                'P12'
                'P13'
                'P22'
                'P23'
                'Q11'
                'Q12'
                'Q13'
                'Q22'
                'Q23'
                'R11'
                'R12'
                'R13'
                'R22'
                'R23'
                };

        case 'ctfheadloc'
            label = {
                'HLC0011'
                'HLC0012'
                'HLC0013'
                'HLC0021'
                'HLC0022'
                'HLC0023'
                'HLC0031'
                'HLC0032'
                'HLC0033'
                'HLC0018'
                'HLC0028'
                'HLC0038'
                'HLC0014'
                'HLC0015'
                'HLC0016'
                'HLC0017'
                'HLC0024'
                'HLC0025'
                'HLC0026'
                'HLC0027'
                'HLC0034'
                'HLC0035'
                'HLC0036'
                'HLC0037'
                };

        case 'ctf64'
            label = {
                'SL11'
                'SL12'
                'SL13'
                'SL14'
                'SL15'
                'SL16'
                'SL17'
                'SL18'
                'SL19'
                'SL21'
                'SL22'
                'SL23'
                'SL24'
                'SL25'
                'SL26'
                'SL27'
                'SL28'
                'SL29'
                'SL31'
                'SL32'
                'SL33'
                'SL34'
                'SL35'
                'SL41'
                'SL42'
                'SL43'
                'SL44'
                'SL45'
                'SL46'
                'SL47'
                'SL51'
                'SL52'
                'SR11'
                'SR12'
                'SR13'
                'SR14'
                'SR15'
                'SR16'
                'SR17'
                'SR18'
                'SR19'
                'SR21'
                'SR22'
                'SR23'
                'SR24'
                'SR25'
                'SR26'
                'SR27'
                'SR28'
                'SR29'
                'SR31'
                'SR32'
                'SR33'
                'SR34'
                'SR35'
                'SR41'
                'SR42'
                'SR43'
                'SR44'
                'SR45'
                'SR46'
                'SR47'
                'SR51'
                'SR52'
                };

        case 'ctf151'
            label = {
                'MLC11'
                'MLC12'
                'MLC13'
                'MLC14'
                'MLC15'
                'MLC21'
                'MLC22'
                'MLC23'
                'MLC24'
                'MLC31'
                'MLC32'
                'MLC33'
                'MLC41'
                'MLC42'
                'MLC43'
                'MLF11'
                'MLF12'
                'MLF21'
                'MLF22'
                'MLF23'
                'MLF31'
                'MLF32'
                'MLF33'
                'MLF34'
                'MLF41'
                'MLF42'
                'MLF43'
                'MLF44'
                'MLF45'
                'MLF51'
                'MLF52'
                'MLO11'
                'MLO12'
                'MLO21'
                'MLO22'
                'MLO31'
                'MLO32'
                'MLO33'
                'MLO41'
                'MLO42'
                'MLO43'
                'MLP11'
                'MLP12'
                'MLP13'
                'MLP21'
                'MLP22'
                'MLP31'
                'MLP32'
                'MLP33'
                'MLP34'
                'MLT11'
                'MLT12'
                'MLT13'
                'MLT14'
                'MLT15'
                'MLT16'
                'MLT21'
                'MLT22'
                'MLT23'
                'MLT24'
                'MLT25'
                'MLT26'
                'MLT31'
                'MLT32'
                'MLT33'
                'MLT34'
                'MLT35'
                'MLT41'
                'MLT42'
                'MLT43'
                'MLT44'
                'MRC11'
                'MRC12'
                'MRC13'
                'MRC14'
                'MRC15'
                'MRC21'
                'MRC22'
                'MRC23'
                'MRC24'
                'MRC31'
                'MRC32'
                'MRC33'
                'MRC41'
                'MRC42'
                'MRC43'
                'MRF11'
                'MRF12'
                'MRF21'
                'MRF22'
                'MRF23'
                'MRF31'
                'MRF32'
                'MRF33'
                'MRF34'
                'MRF41'
                'MRF42'
                'MRF43'
                'MRF44'
                'MRF45'
                'MRF51'
                'MRF52'
                'MRO11'
                'MRO12'
                'MRO21'
                'MRO22'
                'MRO31'
                'MRO32'
                'MRO33'
                'MRO41'
                'MRO42'
                'MRO43'
                'MRP11'
                'MRP12'
                'MRP13'
                'MRP21'
                'MRP22'
                'MRP31'
                'MRP32'
                'MRP33'
                'MRP34'
                'MRT11'
                'MRT12'
                'MRT13'
                'MRT14'
                'MRT15'
                'MRT16'
                'MRT21'
                'MRT22'
                'MRT23'
                'MRT24'
                'MRT25'
                'MRT26'
                'MRT31'
                'MRT32'
                'MRT33'
                'MRT34'
                'MRT35'
                'MRT41'
                'MRT42'
                'MRT43'
                'MRT44'
                'MZC01'
                'MZC02'
                'MZF01'
                'MZF02'
                'MZF03'
                'MZO01'
                'MZO02'
                'MZP01'
                'MZP02'
                };

        case 'ctf151_planar'
            label = {
                'MLC11_dH'  'MLC11_dV'  'MLC11'
                'MLC12_dH'  'MLC12_dV'  'MLC12'
                'MLC13_dH'  'MLC13_dV'  'MLC13'
                'MLC14_dH'  'MLC14_dV'  'MLC14'
                'MLC15_dH'  'MLC15_dV'  'MLC15'
                'MLC21_dH'  'MLC21_dV'  'MLC21'
                'MLC22_dH'  'MLC22_dV'  'MLC22'
                'MLC23_dH'  'MLC23_dV'  'MLC23'
                'MLC24_dH'  'MLC24_dV'  'MLC24'
                'MLC31_dH'  'MLC31_dV'  'MLC31'
                'MLC32_dH'  'MLC32_dV'  'MLC32'
                'MLC33_dH'  'MLC33_dV'  'MLC33'
                'MLC41_dH'  'MLC41_dV'  'MLC41'
                'MLC42_dH'  'MLC42_dV'  'MLC42'
                'MLC43_dH'  'MLC43_dV'  'MLC43'
                'MLF11_dH'  'MLF11_dV'  'MLF11'
                'MLF12_dH'  'MLF12_dV'  'MLF12'
                'MLF21_dH'  'MLF21_dV'  'MLF21'
                'MLF22_dH'  'MLF22_dV'  'MLF22'
                'MLF23_dH'  'MLF23_dV'  'MLF23'
                'MLF31_dH'  'MLF31_dV'  'MLF31'
                'MLF32_dH'  'MLF32_dV'  'MLF32'
                'MLF33_dH'  'MLF33_dV'  'MLF33'
                'MLF34_dH'  'MLF34_dV'  'MLF34'
                'MLF41_dH'  'MLF41_dV'  'MLF41'
                'MLF42_dH'  'MLF42_dV'  'MLF42'
                'MLF43_dH'  'MLF43_dV'  'MLF43'
                'MLF44_dH'  'MLF44_dV'  'MLF44'
                'MLF45_dH'  'MLF45_dV'  'MLF45'
                'MLF51_dH'  'MLF51_dV'  'MLF51'
                'MLF52_dH'  'MLF52_dV'  'MLF52'
                'MLO11_dH'  'MLO11_dV'  'MLO11'
                'MLO12_dH'  'MLO12_dV'  'MLO12'
                'MLO21_dH'  'MLO21_dV'  'MLO21'
                'MLO22_dH'  'MLO22_dV'  'MLO22'
                'MLO31_dH'  'MLO31_dV'  'MLO31'
                'MLO32_dH'  'MLO32_dV'  'MLO32'
                'MLO33_dH'  'MLO33_dV'  'MLO33'
                'MLO41_dH'  'MLO41_dV'  'MLO41'
                'MLO42_dH'  'MLO42_dV'  'MLO42'
                'MLO43_dH'  'MLO43_dV'  'MLO43'
                'MLP11_dH'  'MLP11_dV'  'MLP11'
                'MLP12_dH'  'MLP12_dV'  'MLP12'
                'MLP13_dH'  'MLP13_dV'  'MLP13'
                'MLP21_dH'  'MLP21_dV'  'MLP21'
                'MLP22_dH'  'MLP22_dV'  'MLP22'
                'MLP31_dH'  'MLP31_dV'  'MLP31'
                'MLP32_dH'  'MLP32_dV'  'MLP32'
                'MLP33_dH'  'MLP33_dV'  'MLP33'
                'MLP34_dH'  'MLP34_dV'  'MLP34'
                'MLT11_dH'  'MLT11_dV'  'MLT11'
                'MLT12_dH'  'MLT12_dV'  'MLT12'
                'MLT13_dH'  'MLT13_dV'  'MLT13'
                'MLT14_dH'  'MLT14_dV'  'MLT14'
                'MLT15_dH'  'MLT15_dV'  'MLT15'
                'MLT16_dH'  'MLT16_dV'  'MLT16'
                'MLT21_dH'  'MLT21_dV'  'MLT21'
                'MLT22_dH'  'MLT22_dV'  'MLT22'
                'MLT23_dH'  'MLT23_dV'  'MLT23'
                'MLT24_dH'  'MLT24_dV'  'MLT24'
                'MLT25_dH'  'MLT25_dV'  'MLT25'
                'MLT26_dH'  'MLT26_dV'  'MLT26'
                'MLT31_dH'  'MLT31_dV'  'MLT31'
                'MLT32_dH'  'MLT32_dV'  'MLT32'
                'MLT33_dH'  'MLT33_dV'  'MLT33'
                'MLT34_dH'  'MLT34_dV'  'MLT34'
                'MLT35_dH'  'MLT35_dV'  'MLT35'
                'MLT41_dH'  'MLT41_dV'  'MLT41'
                'MLT42_dH'  'MLT42_dV'  'MLT42'
                'MLT43_dH'  'MLT43_dV'  'MLT43'
                'MLT44_dH'  'MLT44_dV'  'MLT44'
                'MRC11_dH'  'MRC11_dV'  'MRC11'
                'MRC12_dH'  'MRC12_dV'  'MRC12'
                'MRC13_dH'  'MRC13_dV'  'MRC13'
                'MRC14_dH'  'MRC14_dV'  'MRC14'
                'MRC15_dH'  'MRC15_dV'  'MRC15'
                'MRC21_dH'  'MRC21_dV'  'MRC21'
                'MRC22_dH'  'MRC22_dV'  'MRC22'
                'MRC23_dH'  'MRC23_dV'  'MRC23'
                'MRC24_dH'  'MRC24_dV'  'MRC24'
                'MRC31_dH'  'MRC31_dV'  'MRC31'
                'MRC32_dH'  'MRC32_dV'  'MRC32'
                'MRC33_dH'  'MRC33_dV'  'MRC33'
                'MRC41_dH'  'MRC41_dV'  'MRC41'
                'MRC42_dH'  'MRC42_dV'  'MRC42'
                'MRC43_dH'  'MRC43_dV'  'MRC43'
                'MRF11_dH'  'MRF11_dV'  'MRF11'
                'MRF12_dH'  'MRF12_dV'  'MRF12'
                'MRF21_dH'  'MRF21_dV'  'MRF21'
                'MRF22_dH'  'MRF22_dV'  'MRF22'
                'MRF23_dH'  'MRF23_dV'  'MRF23'
                'MRF31_dH'  'MRF31_dV'  'MRF31'
                'MRF32_dH'  'MRF32_dV'  'MRF32'
                'MRF33_dH'  'MRF33_dV'  'MRF33'
                'MRF34_dH'  'MRF34_dV'  'MRF34'
                'MRF41_dH'  'MRF41_dV'  'MRF41'
                'MRF42_dH'  'MRF42_dV'  'MRF42'
                'MRF43_dH'  'MRF43_dV'  'MRF43'
                'MRF44_dH'  'MRF44_dV'  'MRF44'
                'MRF45_dH'  'MRF45_dV'  'MRF45'
                'MRF51_dH'  'MRF51_dV'  'MRF51'
                'MRF52_dH'  'MRF52_dV'  'MRF52'
                'MRO11_dH'  'MRO11_dV'  'MRO11'
                'MRO12_dH'  'MRO12_dV'  'MRO12'
                'MRO21_dH'  'MRO21_dV'  'MRO21'
                'MRO22_dH'  'MRO22_dV'  'MRO22'
                'MRO31_dH'  'MRO31_dV'  'MRO31'
                'MRO32_dH'  'MRO32_dV'  'MRO32'
                'MRO33_dH'  'MRO33_dV'  'MRO33'
                'MRO41_dH'  'MRO41_dV'  'MRO41'
                'MRO42_dH'  'MRO42_dV'  'MRO42'
                'MRO43_dH'  'MRO43_dV'  'MRO43'
                'MRP11_dH'  'MRP11_dV'  'MRP11'
                'MRP12_dH'  'MRP12_dV'  'MRP12'
                'MRP13_dH'  'MRP13_dV'  'MRP13'
                'MRP21_dH'  'MRP21_dV'  'MRP21'
                'MRP22_dH'  'MRP22_dV'  'MRP22'
                'MRP31_dH'  'MRP31_dV'  'MRP31'
                'MRP32_dH'  'MRP32_dV'  'MRP32'
                'MRP33_dH'  'MRP33_dV'  'MRP33'
                'MRP34_dH'  'MRP34_dV'  'MRP34'
                'MRT11_dH'  'MRT11_dV'  'MRT11'
                'MRT12_dH'  'MRT12_dV'  'MRT12'
                'MRT13_dH'  'MRT13_dV'  'MRT13'
                'MRT14_dH'  'MRT14_dV'  'MRT14'
                'MRT15_dH'  'MRT15_dV'  'MRT15'
                'MRT16_dH'  'MRT16_dV'  'MRT16'
                'MRT21_dH'  'MRT21_dV'  'MRT21'
                'MRT22_dH'  'MRT22_dV'  'MRT22'
                'MRT23_dH'  'MRT23_dV'  'MRT23'
                'MRT24_dH'  'MRT24_dV'  'MRT24'
                'MRT25_dH'  'MRT25_dV'  'MRT25'
                'MRT26_dH'  'MRT26_dV'  'MRT26'
                'MRT31_dH'  'MRT31_dV'  'MRT31'
                'MRT32_dH'  'MRT32_dV'  'MRT32'
                'MRT33_dH'  'MRT33_dV'  'MRT33'
                'MRT34_dH'  'MRT34_dV'  'MRT34'
                'MRT35_dH'  'MRT35_dV'  'MRT35'
                'MRT41_dH'  'MRT41_dV'  'MRT41'
                'MRT42_dH'  'MRT42_dV'  'MRT42'
                'MRT43_dH'  'MRT43_dV'  'MRT43'
                'MRT44_dH'  'MRT44_dV'  'MRT44'
                'MZC01_dH'  'MZC01_dV'  'MZC01'
                'MZC02_dH'  'MZC02_dV'  'MZC02'
                'MZF01_dH'  'MZF01_dV'  'MZF01'
                'MZF02_dH'  'MZF02_dV'  'MZF02'
                'MZF03_dH'  'MZF03_dV'  'MZF03'
                'MZO01_dH'  'MZO01_dV'  'MZO01'
                'MZO02_dH'  'MZO02_dV'  'MZO02'
                'MZP01_dH'  'MZP01_dV'  'MZP01'
                'MZP02_dH'  'MZP02_dV'  'MZP02'
                };
            ctf151_planar_combined = label(:,3);
            label = label(:,1:2);

        case 'ctf275'
            label = {
                'MLC11'
                'MLC12'
                'MLC13'
                'MLC14'
                'MLC15'
                'MLC16'
                'MLC17'
                'MLC21'
                'MLC22'
                'MLC23'
                'MLC24'
                'MLC25'
                'MLC31'
                'MLC32'
                'MLC41'
                'MLC42'
                'MLC51'
                'MLC52'
                'MLC53'
                'MLC54'
                'MLC55'
                'MLC61'
                'MLC62'
                'MLC63'
                'MLF11'
                'MLF12'
                'MLF13'
                'MLF14'
                'MLF21'
                'MLF22'
                'MLF23'
                'MLF24'
                'MLF25'
                'MLF31'
                'MLF32'
                'MLF33'
                'MLF34'
                'MLF35'
                'MLF41'
                'MLF42'
                'MLF43'
                'MLF44'
                'MLF45'
                'MLF46'
                'MLF51'
                'MLF52'
                'MLF53'
                'MLF54'
                'MLF55'
                'MLF56'
                'MLF61'
                'MLF62'
                'MLF63'
                'MLF64'
                'MLF65'
                'MLF66'
                'MLF67'
                'MLO11'
                'MLO12'
                'MLO13'
                'MLO14'
                'MLO21'
                'MLO22'
                'MLO23'
                'MLO24'
                'MLO31'
                'MLO32'
                'MLO33'
                'MLO34'
                'MLO41'
                'MLO42'
                'MLO43'
                'MLO44'
                'MLO51'
                'MLO52'
                'MLO53'
                'MLP11'
                'MLP12'
                'MLP21'
                'MLP22'
                'MLP23'
                'MLP31'
                'MLP32'
                'MLP33'
                'MLP34'
                'MLP35'
                'MLP41'
                'MLP42'
                'MLP43'
                'MLP44'
                'MLP45'
                'MLP51'
                'MLP52'
                'MLP53'
                'MLP54'
                'MLP55'
                'MLP56'
                'MLP57'
                'MLT11'
                'MLT12'
                'MLT13'
                'MLT14'
                'MLT15'
                'MLT16'
                'MLT21'
                'MLT22'
                'MLT23'
                'MLT24'
                'MLT25'
                'MLT26'
                'MLT27'
                'MLT31'
                'MLT32'
                'MLT33'
                'MLT34'
                'MLT35'
                'MLT36'
                'MLT37'
                'MLT41'
                'MLT42'
                'MLT43'
                'MLT44'
                'MLT45'
                'MLT46'
                'MLT47'
                'MLT51'
                'MLT52'
                'MLT53'
                'MLT54'
                'MLT55'
                'MLT56'
                'MLT57'
                'MRC11'
                'MRC12'
                'MRC13'
                'MRC14'
                'MRC15'
                'MRC16'
                'MRC17'
                'MRC21'
                'MRC22'
                'MRC23'
                'MRC24'
                'MRC25'
                'MRC31'
                'MRC32'
                'MRC41'
                'MRC42'
                'MRC51'
                'MRC52'
                'MRC53'
                'MRC54'
                'MRC55'
                'MRC61'
                'MRC62'
                'MRC63'
                'MRF11'
                'MRF12'
                'MRF13'
                'MRF14'
                'MRF21'
                'MRF22'
                'MRF23'
                'MRF24'
                'MRF25'
                'MRF31'
                'MRF32'
                'MRF33'
                'MRF34'
                'MRF35'
                'MRF41'
                'MRF42'
                'MRF43'
                'MRF44'
                'MRF45'
                'MRF46'
                'MRF51'
                'MRF52'
                'MRF53'
                'MRF54'
                'MRF55'
                'MRF56'
                'MRF61'
                'MRF62'
                'MRF63'
                'MRF64'
                'MRF65'
                'MRF66'
                'MRF67'
                'MRO11'
                'MRO12'
                'MRO13'
                'MRO14'
                'MRO21'
                'MRO22'
                'MRO23'
                'MRO24'
                'MRO31'
                'MRO32'
                'MRO33'
                'MRO34'
                'MRO41'
                'MRO42'
                'MRO43'
                'MRO44'
                'MRO51'
                'MRO52'
                'MRO53'
                'MRP11'
                'MRP12'
                'MRP21'
                'MRP22'
                'MRP23'
                'MRP31'
                'MRP32'
                'MRP33'
                'MRP34'
                'MRP35'
                'MRP41'
                'MRP42'
                'MRP43'
                'MRP44'
                'MRP45'
                'MRP51'
                'MRP52'
                'MRP53'
                'MRP54'
                'MRP55'
                'MRP56'
                'MRP57'
                'MRT11'
                'MRT12'
                'MRT13'
                'MRT14'
                'MRT15'
                'MRT16'
                'MRT21'
                'MRT22'
                'MRT23'
                'MRT24'
                'MRT25'
                'MRT26'
                'MRT27'
                'MRT31'
                'MRT32'
                'MRT33'
                'MRT34'
                'MRT35'
                'MRT36'
                'MRT37'
                'MRT41'
                'MRT42'
                'MRT43'
                'MRT44'
                'MRT45'
                'MRT46'
                'MRT47'
                'MRT51'
                'MRT52'
                'MRT53'
                'MRT54'
                'MRT55'
                'MRT56'
                'MRT57'
                'MZC01'
                'MZC02'
                'MZC03'
                'MZC04'
                'MZF01'
                'MZF02'
                'MZF03'
                'MZO01'
                'MZO02'
                'MZO03'
                'MZP01'
                };

        case 'ctf275_planar'
            label = {
                'MLC11_dH'  'MLC11_dV'  'MLC11'
                'MLC12_dH'  'MLC12_dV'  'MLC12'
                'MLC13_dH'  'MLC13_dV'  'MLC13'
                'MLC14_dH'  'MLC14_dV'  'MLC14'
                'MLC15_dH'  'MLC15_dV'  'MLC15'
                'MLC16_dH'  'MLC16_dV'  'MLC16'
                'MLC17_dH'  'MLC17_dV'  'MLC17'
                'MLC21_dH'  'MLC21_dV'  'MLC21'
                'MLC22_dH'  'MLC22_dV'  'MLC22'
                'MLC23_dH'  'MLC23_dV'  'MLC23'
                'MLC24_dH'  'MLC24_dV'  'MLC24'
                'MLC25_dH'  'MLC25_dV'  'MLC25'
                'MLC31_dH'  'MLC31_dV'  'MLC31'
                'MLC32_dH'  'MLC32_dV'  'MLC32'
                'MLC41_dH'  'MLC41_dV'  'MLC41'
                'MLC42_dH'  'MLC42_dV'  'MLC42'
                'MLC51_dH'  'MLC51_dV'  'MLC51'
                'MLC52_dH'  'MLC52_dV'  'MLC52'
                'MLC53_dH'  'MLC53_dV'  'MLC53'
                'MLC54_dH'  'MLC54_dV'  'MLC54'
                'MLC55_dH'  'MLC55_dV'  'MLC55'
                'MLC61_dH'  'MLC61_dV'  'MLC61'
                'MLC62_dH'  'MLC62_dV'  'MLC62'
                'MLC63_dH'  'MLC63_dV'  'MLC63'
                'MLF11_dH'  'MLF11_dV'  'MLF11'
                'MLF12_dH'  'MLF12_dV'  'MLF12'
                'MLF13_dH'  'MLF13_dV'  'MLF13'
                'MLF14_dH'  'MLF14_dV'  'MLF14'
                'MLF21_dH'  'MLF21_dV'  'MLF21'
                'MLF22_dH'  'MLF22_dV'  'MLF22'
                'MLF23_dH'  'MLF23_dV'  'MLF23'
                'MLF24_dH'  'MLF24_dV'  'MLF24'
                'MLF25_dH'  'MLF25_dV'  'MLF25'
                'MLF31_dH'  'MLF31_dV'  'MLF31'
                'MLF32_dH'  'MLF32_dV'  'MLF32'
                'MLF33_dH'  'MLF33_dV'  'MLF33'
                'MLF34_dH'  'MLF34_dV'  'MLF34'
                'MLF35_dH'  'MLF35_dV'  'MLF35'
                'MLF41_dH'  'MLF41_dV'  'MLF41'
                'MLF42_dH'  'MLF42_dV'  'MLF42'
                'MLF43_dH'  'MLF43_dV'  'MLF43'
                'MLF44_dH'  'MLF44_dV'  'MLF44'
                'MLF45_dH'  'MLF45_dV'  'MLF45'
                'MLF46_dH'  'MLF46_dV'  'MLF46'
                'MLF51_dH'  'MLF51_dV'  'MLF51'
                'MLF52_dH'  'MLF52_dV'  'MLF52'
                'MLF53_dH'  'MLF53_dV'  'MLF53'
                'MLF54_dH'  'MLF54_dV'  'MLF54'
                'MLF55_dH'  'MLF55_dV'  'MLF55'
                'MLF56_dH'  'MLF56_dV'  'MLF56'
                'MLF61_dH'  'MLF61_dV'  'MLF61'
                'MLF62_dH'  'MLF62_dV'  'MLF62'
                'MLF63_dH'  'MLF63_dV'  'MLF63'
                'MLF64_dH'  'MLF64_dV'  'MLF64'
                'MLF65_dH'  'MLF65_dV'  'MLF65'
                'MLF66_dH'  'MLF66_dV'  'MLF66'
                'MLF67_dH'  'MLF67_dV'  'MLF67'
                'MLO11_dH'  'MLO11_dV'  'MLO11'
                'MLO12_dH'  'MLO12_dV'  'MLO12'
                'MLO13_dH'  'MLO13_dV'  'MLO13'
                'MLO14_dH'  'MLO14_dV'  'MLO14'
                'MLO21_dH'  'MLO21_dV'  'MLO21'
                'MLO22_dH'  'MLO22_dV'  'MLO22'
                'MLO23_dH'  'MLO23_dV'  'MLO23'
                'MLO24_dH'  'MLO24_dV'  'MLO24'
                'MLO31_dH'  'MLO31_dV'  'MLO31'
                'MLO32_dH'  'MLO32_dV'  'MLO32'
                'MLO33_dH'  'MLO33_dV'  'MLO33'
                'MLO34_dH'  'MLO34_dV'  'MLO34'
                'MLO41_dH'  'MLO41_dV'  'MLO41'
                'MLO42_dH'  'MLO42_dV'  'MLO42'
                'MLO43_dH'  'MLO43_dV'  'MLO43'
                'MLO44_dH'  'MLO44_dV'  'MLO44'
                'MLO51_dH'  'MLO51_dV'  'MLO51'
                'MLO52_dH'  'MLO52_dV'  'MLO52'
                'MLO53_dH'  'MLO53_dV'  'MLO53'
                'MLP11_dH'  'MLP11_dV'  'MLP11'
                'MLP12_dH'  'MLP12_dV'  'MLP12'
                'MLP21_dH'  'MLP21_dV'  'MLP21'
                'MLP22_dH'  'MLP22_dV'  'MLP22'
                'MLP23_dH'  'MLP23_dV'  'MLP23'
                'MLP31_dH'  'MLP31_dV'  'MLP31'
                'MLP32_dH'  'MLP32_dV'  'MLP32'
                'MLP33_dH'  'MLP33_dV'  'MLP33'
                'MLP34_dH'  'MLP34_dV'  'MLP34'
                'MLP35_dH'  'MLP35_dV'  'MLP35'
                'MLP41_dH'  'MLP41_dV'  'MLP41'
                'MLP42_dH'  'MLP42_dV'  'MLP42'
                'MLP43_dH'  'MLP43_dV'  'MLP43'
                'MLP44_dH'  'MLP44_dV'  'MLP44'
                'MLP45_dH'  'MLP45_dV'  'MLP45'
                'MLP51_dH'  'MLP51_dV'  'MLP51'
                'MLP52_dH'  'MLP52_dV'  'MLP52'
                'MLP53_dH'  'MLP53_dV'  'MLP53'
                'MLP54_dH'  'MLP54_dV'  'MLP54'
                'MLP55_dH'  'MLP55_dV'  'MLP55'
                'MLP56_dH'  'MLP56_dV'  'MLP56'
                'MLP57_dH'  'MLP57_dV'  'MLP57'
                'MLT11_dH'  'MLT11_dV'  'MLT11'
                'MLT12_dH'  'MLT12_dV'  'MLT12'
                'MLT13_dH'  'MLT13_dV'  'MLT13'
                'MLT14_dH'  'MLT14_dV'  'MLT14'
                'MLT15_dH'  'MLT15_dV'  'MLT15'
                'MLT16_dH'  'MLT16_dV'  'MLT16'
                'MLT21_dH'  'MLT21_dV'  'MLT21'
                'MLT22_dH'  'MLT22_dV'  'MLT22'
                'MLT23_dH'  'MLT23_dV'  'MLT23'
                'MLT24_dH'  'MLT24_dV'  'MLT24'
                'MLT25_dH'  'MLT25_dV'  'MLT25'
                'MLT26_dH'  'MLT26_dV'  'MLT26'
                'MLT27_dH'  'MLT27_dV'  'MLT27'
                'MLT31_dH'  'MLT31_dV'  'MLT31'
                'MLT32_dH'  'MLT32_dV'  'MLT32'
                'MLT33_dH'  'MLT33_dV'  'MLT33'
                'MLT34_dH'  'MLT34_dV'  'MLT34'
                'MLT35_dH'  'MLT35_dV'  'MLT35'
                'MLT36_dH'  'MLT36_dV'  'MLT36'
                'MLT37_dH'  'MLT37_dV'  'MLT37'
                'MLT41_dH'  'MLT41_dV'  'MLT41'
                'MLT42_dH'  'MLT42_dV'  'MLT42'
                'MLT43_dH'  'MLT43_dV'  'MLT43'
                'MLT44_dH'  'MLT44_dV'  'MLT44'
                'MLT45_dH'  'MLT45_dV'  'MLT45'
                'MLT46_dH'  'MLT46_dV'  'MLT46'
                'MLT47_dH'  'MLT47_dV'  'MLT47'
                'MLT51_dH'  'MLT51_dV'  'MLT51'
                'MLT52_dH'  'MLT52_dV'  'MLT52'
                'MLT53_dH'  'MLT53_dV'  'MLT53'
                'MLT54_dH'  'MLT54_dV'  'MLT54'
                'MLT55_dH'  'MLT55_dV'  'MLT55'
                'MLT56_dH'  'MLT56_dV'  'MLT56'
                'MLT57_dH'  'MLT57_dV'  'MLT57'
                'MRC11_dH'  'MRC11_dV'  'MRC11'
                'MRC12_dH'  'MRC12_dV'  'MRC12'
                'MRC13_dH'  'MRC13_dV'  'MRC13'
                'MRC14_dH'  'MRC14_dV'  'MRC14'
                'MRC15_dH'  'MRC15_dV'  'MRC15'
                'MRC16_dH'  'MRC16_dV'  'MRC16'
                'MRC17_dH'  'MRC17_dV'  'MRC17'
                'MRC21_dH'  'MRC21_dV'  'MRC21'
                'MRC22_dH'  'MRC22_dV'  'MRC22'
                'MRC23_dH'  'MRC23_dV'  'MRC23'
                'MRC24_dH'  'MRC24_dV'  'MRC24'
                'MRC25_dH'  'MRC25_dV'  'MRC25'
                'MRC31_dH'  'MRC31_dV'  'MRC31'
                'MRC32_dH'  'MRC32_dV'  'MRC32'
                'MRC41_dH'  'MRC41_dV'  'MRC41'
                'MRC42_dH'  'MRC42_dV'  'MRC42'
                'MRC51_dH'  'MRC51_dV'  'MRC51'
                'MRC52_dH'  'MRC52_dV'  'MRC52'
                'MRC53_dH'  'MRC53_dV'  'MRC53'
                'MRC54_dH'  'MRC54_dV'  'MRC54'
                'MRC55_dH'  'MRC55_dV'  'MRC55'
                'MRC61_dH'  'MRC61_dV'  'MRC61'
                'MRC62_dH'  'MRC62_dV'  'MRC62'
                'MRC63_dH'  'MRC63_dV'  'MRC63'
                'MRF11_dH'  'MRF11_dV'  'MRF11'
                'MRF12_dH'  'MRF12_dV'  'MRF12'
                'MRF13_dH'  'MRF13_dV'  'MRF13'
                'MRF14_dH'  'MRF14_dV'  'MRF14'
                'MRF21_dH'  'MRF21_dV'  'MRF21'
                'MRF22_dH'  'MRF22_dV'  'MRF22'
                'MRF23_dH'  'MRF23_dV'  'MRF23'
                'MRF24_dH'  'MRF24_dV'  'MRF24'
                'MRF25_dH'  'MRF25_dV'  'MRF25'
                'MRF31_dH'  'MRF31_dV'  'MRF31'
                'MRF32_dH'  'MRF32_dV'  'MRF32'
                'MRF33_dH'  'MRF33_dV'  'MRF33'
                'MRF34_dH'  'MRF34_dV'  'MRF34'
                'MRF35_dH'  'MRF35_dV'  'MRF35'
                'MRF41_dH'  'MRF41_dV'  'MRF41'
                'MRF42_dH'  'MRF42_dV'  'MRF42'
                'MRF43_dH'  'MRF43_dV'  'MRF43'
                'MRF44_dH'  'MRF44_dV'  'MRF44'
                'MRF45_dH'  'MRF45_dV'  'MRF45'
                'MRF46_dH'  'MRF46_dV'  'MRF46'
                'MRF51_dH'  'MRF51_dV'  'MRF51'
                'MRF52_dH'  'MRF52_dV'  'MRF52'
                'MRF53_dH'  'MRF53_dV'  'MRF53'
                'MRF54_dH'  'MRF54_dV'  'MRF54'
                'MRF55_dH'  'MRF55_dV'  'MRF55'
                'MRF56_dH'  'MRF56_dV'  'MRF56'
                'MRF61_dH'  'MRF61_dV'  'MRF61'
                'MRF62_dH'  'MRF62_dV'  'MRF62'
                'MRF63_dH'  'MRF63_dV'  'MRF63'
                'MRF64_dH'  'MRF64_dV'  'MRF64'
                'MRF65_dH'  'MRF65_dV'  'MRF65'
                'MRF66_dH'  'MRF66_dV'  'MRF66'
                'MRF67_dH'  'MRF67_dV'  'MRF67'
                'MRO11_dH'  'MRO11_dV'  'MRO11'
                'MRO12_dH'  'MRO12_dV'  'MRO12'
                'MRO13_dH'  'MRO13_dV'  'MRO13'
                'MRO14_dH'  'MRO14_dV'  'MRO14'
                'MRO21_dH'  'MRO21_dV'  'MRO21'
                'MRO22_dH'  'MRO22_dV'  'MRO22'
                'MRO23_dH'  'MRO23_dV'  'MRO23'
                'MRO24_dH'  'MRO24_dV'  'MRO24'
                'MRO31_dH'  'MRO31_dV'  'MRO31'
                'MRO32_dH'  'MRO32_dV'  'MRO32'
                'MRO33_dH'  'MRO33_dV'  'MRO33'
                'MRO34_dH'  'MRO34_dV'  'MRO34'
                'MRO41_dH'  'MRO41_dV'  'MRO41'
                'MRO42_dH'  'MRO42_dV'  'MRO42'
                'MRO43_dH'  'MRO43_dV'  'MRO43'
                'MRO44_dH'  'MRO44_dV'  'MRO44'
                'MRO51_dH'  'MRO51_dV'  'MRO51'
                'MRO52_dH'  'MRO52_dV'  'MRO52'
                'MRO53_dH'  'MRO53_dV'  'MRO53'
                'MRP11_dH'  'MRP11_dV'  'MRP11'
                'MRP12_dH'  'MRP12_dV'  'MRP12'
                'MRP21_dH'  'MRP21_dV'  'MRP21'
                'MRP22_dH'  'MRP22_dV'  'MRP22'
                'MRP23_dH'  'MRP23_dV'  'MRP23'
                'MRP31_dH'  'MRP31_dV'  'MRP31'
                'MRP32_dH'  'MRP32_dV'  'MRP32'
                'MRP33_dH'  'MRP33_dV'  'MRP33'
                'MRP34_dH'  'MRP34_dV'  'MRP34'
                'MRP35_dH'  'MRP35_dV'  'MRP35'
                'MRP41_dH'  'MRP41_dV'  'MRP41'
                'MRP42_dH'  'MRP42_dV'  'MRP42'
                'MRP43_dH'  'MRP43_dV'  'MRP43'
                'MRP44_dH'  'MRP44_dV'  'MRP44'
                'MRP45_dH'  'MRP45_dV'  'MRP45'
                'MRP51_dH'  'MRP51_dV'  'MRP51'
                'MRP52_dH'  'MRP52_dV'  'MRP52'
                'MRP53_dH'  'MRP53_dV'  'MRP53'
                'MRP54_dH'  'MRP54_dV'  'MRP54'
                'MRP55_dH'  'MRP55_dV'  'MRP55'
                'MRP56_dH'  'MRP56_dV'  'MRP56'
                'MRP57_dH'  'MRP57_dV'  'MRP57'
                'MRT11_dH'  'MRT11_dV'  'MRT11'
                'MRT12_dH'  'MRT12_dV'  'MRT12'
                'MRT13_dH'  'MRT13_dV'  'MRT13'
                'MRT14_dH'  'MRT14_dV'  'MRT14'
                'MRT15_dH'  'MRT15_dV'  'MRT15'
                'MRT16_dH'  'MRT16_dV'  'MRT16'
                'MRT21_dH'  'MRT21_dV'  'MRT21'
                'MRT22_dH'  'MRT22_dV'  'MRT22'
                'MRT23_dH'  'MRT23_dV'  'MRT23'
                'MRT24_dH'  'MRT24_dV'  'MRT24'
                'MRT25_dH'  'MRT25_dV'  'MRT25'
                'MRT26_dH'  'MRT26_dV'  'MRT26'
                'MRT27_dH'  'MRT27_dV'  'MRT27'
                'MRT31_dH'  'MRT31_dV'  'MRT31'
                'MRT32_dH'  'MRT32_dV'  'MRT32'
                'MRT33_dH'  'MRT33_dV'  'MRT33'
                'MRT34_dH'  'MRT34_dV'  'MRT34'
                'MRT35_dH'  'MRT35_dV'  'MRT35'
                'MRT36_dH'  'MRT36_dV'  'MRT36'
                'MRT37_dH'  'MRT37_dV'  'MRT37'
                'MRT41_dH'  'MRT41_dV'  'MRT41'
                'MRT42_dH'  'MRT42_dV'  'MRT42'
                'MRT43_dH'  'MRT43_dV'  'MRT43'
                'MRT44_dH'  'MRT44_dV'  'MRT44'
                'MRT45_dH'  'MRT45_dV'  'MRT45'
                'MRT46_dH'  'MRT46_dV'  'MRT46'
                'MRT47_dH'  'MRT47_dV'  'MRT47'
                'MRT51_dH'  'MRT51_dV'  'MRT51'
                'MRT52_dH'  'MRT52_dV'  'MRT52'
                'MRT53_dH'  'MRT53_dV'  'MRT53'
                'MRT54_dH'  'MRT54_dV'  'MRT54'
                'MRT55_dH'  'MRT55_dV'  'MRT55'
                'MRT56_dH'  'MRT56_dV'  'MRT56'
                'MRT57_dH'  'MRT57_dV'  'MRT57'
                'MZC01_dH'  'MZC01_dV'  'MZC01'
                'MZC02_dH'  'MZC02_dV'  'MZC02'
                'MZC03_dH'  'MZC03_dV'  'MZC03'
                'MZC04_dH'  'MZC04_dV'  'MZC04'
                'MZF01_dH'  'MZF01_dV'  'MZF01'
                'MZF02_dH'  'MZF02_dV'  'MZF02'
                'MZF03_dH'  'MZF03_dV'  'MZF03'
                'MZO01_dH'  'MZO01_dV'  'MZO01'
                'MZO02_dH'  'MZO02_dV'  'MZO02'
                'MZO03_dH'  'MZO03_dV'  'MZO03'
                'MZP01_dH'  'MZP01_dV'  'MZP01'
                };
            ctf275_planar_combined = label(:,3);
            label = label(:,1:2);

        case {'neuromag122' 'neuromag122alt'}
            % this is the combination of the two versions (with and without space)
            label = {
                'MEG 001'  'MEG 002'  'MEG 001+002'
                'MEG 003'  'MEG 004'  'MEG 003+004'
                'MEG 005'  'MEG 006'  'MEG 005+006'
                'MEG 007'  'MEG 008'  'MEG 007+008'
                'MEG 009'  'MEG 010'  'MEG 009+010'
                'MEG 011'  'MEG 012'  'MEG 011+012'
                'MEG 013'  'MEG 014'  'MEG 013+014'
                'MEG 015'  'MEG 016'  'MEG 015+016'
                'MEG 017'  'MEG 018'  'MEG 017+018'
                'MEG 019'  'MEG 020'  'MEG 019+020'
                'MEG 021'  'MEG 022'  'MEG 021+022'
                'MEG 023'  'MEG 024'  'MEG 023+024'
                'MEG 025'  'MEG 026'  'MEG 025+026'
                'MEG 027'  'MEG 028'  'MEG 027+028'
                'MEG 029'  'MEG 030'  'MEG 029+030'
                'MEG 031'  'MEG 032'  'MEG 031+032'
                'MEG 033'  'MEG 034'  'MEG 033+034'
                'MEG 035'  'MEG 036'  'MEG 035+036'
                'MEG 037'  'MEG 038'  'MEG 037+038'
                'MEG 039'  'MEG 040'  'MEG 039+040'
                'MEG 041'  'MEG 042'  'MEG 041+042'
                'MEG 043'  'MEG 044'  'MEG 043+044'
                'MEG 045'  'MEG 046'  'MEG 045+046'
                'MEG 047'  'MEG 048'  'MEG 047+048'
                'MEG 049'  'MEG 050'  'MEG 049+050'
                'MEG 051'  'MEG 052'  'MEG 051+052'
                'MEG 053'  'MEG 054'  'MEG 053+054'
                'MEG 055'  'MEG 056'  'MEG 055+056'
                'MEG 057'  'MEG 058'  'MEG 057+058'
                'MEG 059'  'MEG 060'  'MEG 059+060'
                'MEG 061'  'MEG 062'  'MEG 061+062'
                'MEG 063'  'MEG 064'  'MEG 063+064'
                'MEG 065'  'MEG 066'  'MEG 065+066'
                'MEG 067'  'MEG 068'  'MEG 067+068'
                'MEG 069'  'MEG 070'  'MEG 069+070'
                'MEG 071'  'MEG 072'  'MEG 071+072'
                'MEG 073'  'MEG 074'  'MEG 073+074'
                'MEG 075'  'MEG 076'  'MEG 075+076'
                'MEG 077'  'MEG 078'  'MEG 077+078'
                'MEG 079'  'MEG 080'  'MEG 079+080'
                'MEG 081'  'MEG 082'  'MEG 081+082'
                'MEG 083'  'MEG 084'  'MEG 083+084'
                'MEG 085'  'MEG 086'  'MEG 085+086'
                'MEG 087'  'MEG 088'  'MEG 087+088'
                'MEG 089'  'MEG 090'  'MEG 089+090'
                'MEG 091'  'MEG 092'  'MEG 091+092'
                'MEG 093'  'MEG 094'  'MEG 093+094'
                'MEG 095'  'MEG 096'  'MEG 095+096'
                'MEG 097'  'MEG 098'  'MEG 097+098'
                'MEG 099'  'MEG 100'  'MEG 099+100'
                'MEG 101'  'MEG 102'  'MEG 101+102'
                'MEG 103'  'MEG 104'  'MEG 103+104'
                'MEG 105'  'MEG 106'  'MEG 105+106'
                'MEG 107'  'MEG 108'  'MEG 107+108'
                'MEG 109'  'MEG 110'  'MEG 109+110'
                'MEG 111'  'MEG 112'  'MEG 111+112'
                'MEG 113'  'MEG 114'  'MEG 113+114'
                'MEG 115'  'MEG 116'  'MEG 115+116'
                'MEG 117'  'MEG 118'  'MEG 117+118'
                'MEG 119'  'MEG 120'  'MEG 119+120'
                'MEG 121'  'MEG 122'  'MEG 121+122'
                % this is an alternative set of labels without a space in them
                'MEG001'  'MEG002'  'MEG001+002'
                'MEG003'  'MEG004'  'MEG003+004'
                'MEG005'  'MEG006'  'MEG005+006'
                'MEG007'  'MEG008'  'MEG007+008'
                'MEG009'  'MEG010'  'MEG009+010'
                'MEG011'  'MEG012'  'MEG011+012'
                'MEG013'  'MEG014'  'MEG013+014'
                'MEG015'  'MEG016'  'MEG015+016'
                'MEG017'  'MEG018'  'MEG017+018'
                'MEG019'  'MEG020'  'MEG019+020'
                'MEG021'  'MEG022'  'MEG021+022'
                'MEG023'  'MEG024'  'MEG023+024'
                'MEG025'  'MEG026'  'MEG025+026'
                'MEG027'  'MEG028'  'MEG027+028'
                'MEG029'  'MEG030'  'MEG029+030'
                'MEG031'  'MEG032'  'MEG031+032'
                'MEG033'  'MEG034'  'MEG033+034'
                'MEG035'  'MEG036'  'MEG035+036'
                'MEG037'  'MEG038'  'MEG037+038'
                'MEG039'  'MEG040'  'MEG039+040'
                'MEG041'  'MEG042'  'MEG041+042'
                'MEG043'  'MEG044'  'MEG043+044'
                'MEG045'  'MEG046'  'MEG045+046'
                'MEG047'  'MEG048'  'MEG047+048'
                'MEG049'  'MEG050'  'MEG049+050'
                'MEG051'  'MEG052'  'MEG051+052'
                'MEG053'  'MEG054'  'MEG053+054'
                'MEG055'  'MEG056'  'MEG055+056'
                'MEG057'  'MEG058'  'MEG057+058'
                'MEG059'  'MEG060'  'MEG059+060'
                'MEG061'  'MEG062'  'MEG061+062'
                'MEG063'  'MEG064'  'MEG063+064'
                'MEG065'  'MEG066'  'MEG065+066'
                'MEG067'  'MEG068'  'MEG067+068'
                'MEG069'  'MEG070'  'MEG069+070'
                'MEG071'  'MEG072'  'MEG071+072'
                'MEG073'  'MEG074'  'MEG073+074'
                'MEG075'  'MEG076'  'MEG075+076'
                'MEG077'  'MEG078'  'MEG077+078'
                'MEG079'  'MEG080'  'MEG079+080'
                'MEG081'  'MEG082'  'MEG081+082'
                'MEG083'  'MEG084'  'MEG083+084'
                'MEG085'  'MEG086'  'MEG085+086'
                'MEG087'  'MEG088'  'MEG087+088'
                'MEG089'  'MEG090'  'MEG089+090'
                'MEG091'  'MEG092'  'MEG091+092'
                'MEG093'  'MEG094'  'MEG093+094'
                'MEG095'  'MEG096'  'MEG095+096'
                'MEG097'  'MEG098'  'MEG097+098'
                'MEG099'  'MEG100'  'MEG099+100'
                'MEG101'  'MEG102'  'MEG101+102'
                'MEG103'  'MEG104'  'MEG103+104'
                'MEG105'  'MEG106'  'MEG105+106'
                'MEG107'  'MEG108'  'MEG107+108'
                'MEG109'  'MEG110'  'MEG109+110'
                'MEG111'  'MEG112'  'MEG111+112'
                'MEG113'  'MEG114'  'MEG113+114'
                'MEG115'  'MEG116'  'MEG115+116'
                'MEG117'  'MEG118'  'MEG117+118'
                'MEG119'  'MEG120'  'MEG119+120'
                'MEG121'  'MEG122'  'MEG121+122'
                };
            neuromag122_combined = label(:,3);
            neuromag122alt_combined = label(:,3);
            label = label(:,1:2);

        case {'neuromag306' 'neuromag306alt'}
            % this is the combination of the two versions (with and without space)
            label = {
                'MEG 0112'  'MEG 0113'  'MEG 0111'  'MEG 0112+0113'
                'MEG 0122'  'MEG 0123'  'MEG 0121'  'MEG 0122+0123'
                'MEG 0132'  'MEG 0133'  'MEG 0131'  'MEG 0132+0133'
                'MEG 0142'  'MEG 0143'  'MEG 0141'  'MEG 0142+0143'
                'MEG 0212'  'MEG 0213'  'MEG 0211'  'MEG 0212+0213'
                'MEG 0222'  'MEG 0223'  'MEG 0221'  'MEG 0222+0223'
                'MEG 0232'  'MEG 0233'  'MEG 0231'  'MEG 0232+0233'
                'MEG 0242'  'MEG 0243'  'MEG 0241'  'MEG 0242+0243'
                'MEG 0312'  'MEG 0313'  'MEG 0311'  'MEG 0312+0313'
                'MEG 0322'  'MEG 0323'  'MEG 0321'  'MEG 0322+0323'
                'MEG 0332'  'MEG 0333'  'MEG 0331'  'MEG 0332+0333'
                'MEG 0342'  'MEG 0343'  'MEG 0341'  'MEG 0342+0343'
                'MEG 0412'  'MEG 0413'  'MEG 0411'  'MEG 0412+0413'
                'MEG 0422'  'MEG 0423'  'MEG 0421'  'MEG 0422+0423'
                'MEG 0432'  'MEG 0433'  'MEG 0431'  'MEG 0432+0433'
                'MEG 0442'  'MEG 0443'  'MEG 0441'  'MEG 0442+0443'
                'MEG 0512'  'MEG 0513'  'MEG 0511'  'MEG 0512+0513'
                'MEG 0522'  'MEG 0523'  'MEG 0521'  'MEG 0522+0523'
                'MEG 0532'  'MEG 0533'  'MEG 0531'  'MEG 0532+0533'
                'MEG 0542'  'MEG 0543'  'MEG 0541'  'MEG 0542+0543'
                'MEG 0612'  'MEG 0613'  'MEG 0611'  'MEG 0612+0613'
                'MEG 0622'  'MEG 0623'  'MEG 0621'  'MEG 0622+0623'
                'MEG 0632'  'MEG 0633'  'MEG 0631'  'MEG 0632+0633'
                'MEG 0642'  'MEG 0643'  'MEG 0641'  'MEG 0642+0643'
                'MEG 0712'  'MEG 0713'  'MEG 0711'  'MEG 0712+0713'
                'MEG 0722'  'MEG 0723'  'MEG 0721'  'MEG 0722+0723'
                'MEG 0732'  'MEG 0733'  'MEG 0731'  'MEG 0732+0733'
                'MEG 0742'  'MEG 0743'  'MEG 0741'  'MEG 0742+0743'
                'MEG 0812'  'MEG 0813'  'MEG 0811'  'MEG 0812+0813'
                'MEG 0822'  'MEG 0823'  'MEG 0821'  'MEG 0822+0823'
                'MEG 0912'  'MEG 0913'  'MEG 0911'  'MEG 0912+0913'
                'MEG 0922'  'MEG 0923'  'MEG 0921'  'MEG 0922+0923'
                'MEG 0932'  'MEG 0933'  'MEG 0931'  'MEG 0932+0933'
                'MEG 0942'  'MEG 0943'  'MEG 0941'  'MEG 0942+0943'
                'MEG 1012'  'MEG 1013'  'MEG 1011'  'MEG 1012+1013'
                'MEG 1022'  'MEG 1023'  'MEG 1021'  'MEG 1022+1023'
                'MEG 1032'  'MEG 1033'  'MEG 1031'  'MEG 1032+1033'
                'MEG 1042'  'MEG 1043'  'MEG 1041'  'MEG 1042+1043'
                'MEG 1112'  'MEG 1113'  'MEG 1111'  'MEG 1112+1113'
                'MEG 1122'  'MEG 1123'  'MEG 1121'  'MEG 1122+1123'
                'MEG 1132'  'MEG 1133'  'MEG 1131'  'MEG 1132+1133'
                'MEG 1142'  'MEG 1143'  'MEG 1141'  'MEG 1142+1143'
                'MEG 1212'  'MEG 1213'  'MEG 1211'  'MEG 1212+1213'
                'MEG 1222'  'MEG 1223'  'MEG 1221'  'MEG 1222+1223'
                'MEG 1232'  'MEG 1233'  'MEG 1231'  'MEG 1232+1233'
                'MEG 1242'  'MEG 1243'  'MEG 1241'  'MEG 1242+1243'
                'MEG 1312'  'MEG 1313'  'MEG 1311'  'MEG 1312+1313'
                'MEG 1322'  'MEG 1323'  'MEG 1321'  'MEG 1322+1323'
                'MEG 1332'  'MEG 1333'  'MEG 1331'  'MEG 1332+1333'
                'MEG 1342'  'MEG 1343'  'MEG 1341'  'MEG 1342+1343'
                'MEG 1412'  'MEG 1413'  'MEG 1411'  'MEG 1412+1413'
                'MEG 1422'  'MEG 1423'  'MEG 1421'  'MEG 1422+1423'
                'MEG 1432'  'MEG 1433'  'MEG 1431'  'MEG 1432+1433'
                'MEG 1442'  'MEG 1443'  'MEG 1441'  'MEG 1442+1443'
                'MEG 1512'  'MEG 1513'  'MEG 1511'  'MEG 1512+1513'
                'MEG 1522'  'MEG 1523'  'MEG 1521'  'MEG 1522+1523'
                'MEG 1532'  'MEG 1533'  'MEG 1531'  'MEG 1532+1533'
                'MEG 1542'  'MEG 1543'  'MEG 1541'  'MEG 1542+1543'
                'MEG 1612'  'MEG 1613'  'MEG 1611'  'MEG 1612+1613'
                'MEG 1622'  'MEG 1623'  'MEG 1621'  'MEG 1622+1623'
                'MEG 1632'  'MEG 1633'  'MEG 1631'  'MEG 1632+1633'
                'MEG 1642'  'MEG 1643'  'MEG 1641'  'MEG 1642+1643'
                'MEG 1712'  'MEG 1713'  'MEG 1711'  'MEG 1712+1713'
                'MEG 1722'  'MEG 1723'  'MEG 1721'  'MEG 1722+1723'
                'MEG 1732'  'MEG 1733'  'MEG 1731'  'MEG 1732+1733'
                'MEG 1742'  'MEG 1743'  'MEG 1741'  'MEG 1742+1743'
                'MEG 1812'  'MEG 1813'  'MEG 1811'  'MEG 1812+1813'
                'MEG 1822'  'MEG 1823'  'MEG 1821'  'MEG 1822+1823'
                'MEG 1832'  'MEG 1833'  'MEG 1831'  'MEG 1832+1833'
                'MEG 1842'  'MEG 1843'  'MEG 1841'  'MEG 1842+1843'
                'MEG 1912'  'MEG 1913'  'MEG 1911'  'MEG 1912+1913'
                'MEG 1922'  'MEG 1923'  'MEG 1921'  'MEG 1922+1923'
                'MEG 1932'  'MEG 1933'  'MEG 1931'  'MEG 1932+1933'
                'MEG 1942'  'MEG 1943'  'MEG 1941'  'MEG 1942+1943'
                'MEG 2012'  'MEG 2013'  'MEG 2011'  'MEG 2012+2013'
                'MEG 2022'  'MEG 2023'  'MEG 2021'  'MEG 2022+2023'
                'MEG 2032'  'MEG 2033'  'MEG 2031'  'MEG 2032+2033'
                'MEG 2042'  'MEG 2043'  'MEG 2041'  'MEG 2042+2043'
                'MEG 2112'  'MEG 2113'  'MEG 2111'  'MEG 2112+2113'
                'MEG 2122'  'MEG 2123'  'MEG 2121'  'MEG 2122+2123'
                'MEG 2132'  'MEG 2133'  'MEG 2131'  'MEG 2132+2133'
                'MEG 2142'  'MEG 2143'  'MEG 2141'  'MEG 2142+2143'
                'MEG 2212'  'MEG 2213'  'MEG 2211'  'MEG 2212+2213'
                'MEG 2222'  'MEG 2223'  'MEG 2221'  'MEG 2222+2223'
                'MEG 2232'  'MEG 2233'  'MEG 2231'  'MEG 2232+2233'
                'MEG 2242'  'MEG 2243'  'MEG 2241'  'MEG 2242+2243'
                'MEG 2312'  'MEG 2313'  'MEG 2311'  'MEG 2312+2313'
                'MEG 2322'  'MEG 2323'  'MEG 2321'  'MEG 2322+2323'
                'MEG 2332'  'MEG 2333'  'MEG 2331'  'MEG 2332+2333'
                'MEG 2342'  'MEG 2343'  'MEG 2341'  'MEG 2342+2343'
                'MEG 2412'  'MEG 2413'  'MEG 2411'  'MEG 2412+2413'
                'MEG 2422'  'MEG 2423'  'MEG 2421'  'MEG 2422+2423'
                'MEG 2432'  'MEG 2433'  'MEG 2431'  'MEG 2432+2433'
                'MEG 2442'  'MEG 2443'  'MEG 2441'  'MEG 2442+2443'
                'MEG 2512'  'MEG 2513'  'MEG 2511'  'MEG 2512+2513'
                'MEG 2522'  'MEG 2523'  'MEG 2521'  'MEG 2522+2523'
                'MEG 2532'  'MEG 2533'  'MEG 2531'  'MEG 2532+2533'
                'MEG 2542'  'MEG 2543'  'MEG 2541'  'MEG 2542+2543'
                'MEG 2612'  'MEG 2613'  'MEG 2611'  'MEG 2612+2613'
                'MEG 2622'  'MEG 2623'  'MEG 2621'  'MEG 2622+2623'
                'MEG 2632'  'MEG 2633'  'MEG 2631'  'MEG 2632+2633'
                'MEG 2642'  'MEG 2643'  'MEG 2641'  'MEG 2642+2643'
                % this is an alternative set of labels without a space in them
                'MEG0112'  'MEG0113'  'MEG0111'  'MEG0112+0113'
                'MEG0122'  'MEG0123'  'MEG0121'  'MEG0122+0123'
                'MEG0132'  'MEG0133'  'MEG0131'  'MEG0132+0133'
                'MEG0142'  'MEG0143'  'MEG0141'  'MEG0142+0143'
                'MEG0212'  'MEG0213'  'MEG0211'  'MEG0212+0213'
                'MEG0222'  'MEG0223'  'MEG0221'  'MEG0222+0223'
                'MEG0232'  'MEG0233'  'MEG0231'  'MEG0232+0233'
                'MEG0242'  'MEG0243'  'MEG0241'  'MEG0242+0243'
                'MEG0312'  'MEG0313'  'MEG0311'  'MEG0312+0313'
                'MEG0322'  'MEG0323'  'MEG0321'  'MEG0322+0323'
                'MEG0332'  'MEG0333'  'MEG0331'  'MEG0332+0333'
                'MEG0342'  'MEG0343'  'MEG0341'  'MEG0342+0343'
                'MEG0412'  'MEG0413'  'MEG0411'  'MEG0412+0413'
                'MEG0422'  'MEG0423'  'MEG0421'  'MEG0422+0423'
                'MEG0432'  'MEG0433'  'MEG0431'  'MEG0432+0433'
                'MEG0442'  'MEG0443'  'MEG0441'  'MEG0442+0443'
                'MEG0512'  'MEG0513'  'MEG0511'  'MEG0512+0513'
                'MEG0522'  'MEG0523'  'MEG0521'  'MEG0522+0523'
                'MEG0532'  'MEG0533'  'MEG0531'  'MEG0532+0533'
                'MEG0542'  'MEG0543'  'MEG0541'  'MEG0542+0543'
                'MEG0612'  'MEG0613'  'MEG0611'  'MEG0612+0613'
                'MEG0622'  'MEG0623'  'MEG0621'  'MEG0622+0623'
                'MEG0632'  'MEG0633'  'MEG0631'  'MEG0632+0633'
                'MEG0642'  'MEG0643'  'MEG0641'  'MEG0642+0643'
                'MEG0712'  'MEG0713'  'MEG0711'  'MEG0712+0713'
                'MEG0722'  'MEG0723'  'MEG0721'  'MEG0722+0723'
                'MEG0732'  'MEG0733'  'MEG0731'  'MEG0732+0733'
                'MEG0742'  'MEG0743'  'MEG0741'  'MEG0742+0743'
                'MEG0812'  'MEG0813'  'MEG0811'  'MEG0812+0813'
                'MEG0822'  'MEG0823'  'MEG0821'  'MEG0822+0823'
                'MEG0912'  'MEG0913'  'MEG0911'  'MEG0912+0913'
                'MEG0922'  'MEG0923'  'MEG0921'  'MEG0922+0923'
                'MEG0932'  'MEG0933'  'MEG0931'  'MEG0932+0933'
                'MEG0942'  'MEG0943'  'MEG0941'  'MEG0942+0943'
                'MEG1012'  'MEG1013'  'MEG1011'  'MEG1012+1013'
                'MEG1022'  'MEG1023'  'MEG1021'  'MEG1022+1023'
                'MEG1032'  'MEG1033'  'MEG1031'  'MEG1032+1033'
                'MEG1042'  'MEG1043'  'MEG1041'  'MEG1042+1043'
                'MEG1112'  'MEG1113'  'MEG1111'  'MEG1112+1113'
                'MEG1122'  'MEG1123'  'MEG1121'  'MEG1122+1123'
                'MEG1132'  'MEG1133'  'MEG1131'  'MEG1132+1133'
                'MEG1142'  'MEG1143'  'MEG1141'  'MEG1142+1143'
                'MEG1212'  'MEG1213'  'MEG1211'  'MEG1212+1213'
                'MEG1222'  'MEG1223'  'MEG1221'  'MEG1222+1223'
                'MEG1232'  'MEG1233'  'MEG1231'  'MEG1232+1233'
                'MEG1242'  'MEG1243'  'MEG1241'  'MEG1242+1243'
                'MEG1312'  'MEG1313'  'MEG1311'  'MEG1312+1313'
                'MEG1322'  'MEG1323'  'MEG1321'  'MEG1322+1323'
                'MEG1332'  'MEG1333'  'MEG1331'  'MEG1332+1333'
                'MEG1342'  'MEG1343'  'MEG1341'  'MEG1342+1343'
                'MEG1412'  'MEG1413'  'MEG1411'  'MEG1412+1413'
                'MEG1422'  'MEG1423'  'MEG1421'  'MEG1422+1423'
                'MEG1432'  'MEG1433'  'MEG1431'  'MEG1432+1433'
                'MEG1442'  'MEG1443'  'MEG1441'  'MEG1442+1443'
                'MEG1512'  'MEG1513'  'MEG1511'  'MEG1512+1513'
                'MEG1522'  'MEG1523'  'MEG1521'  'MEG1522+1523'
                'MEG1532'  'MEG1533'  'MEG1531'  'MEG1532+1533'
                'MEG1542'  'MEG1543'  'MEG1541'  'MEG1542+1543'
                'MEG1612'  'MEG1613'  'MEG1611'  'MEG1612+1613'
                'MEG1622'  'MEG1623'  'MEG1621'  'MEG1622+1623'
                'MEG1632'  'MEG1633'  'MEG1631'  'MEG1632+1633'
                'MEG1642'  'MEG1643'  'MEG1641'  'MEG1642+1643'
                'MEG1712'  'MEG1713'  'MEG1711'  'MEG1712+1713'
                'MEG1722'  'MEG1723'  'MEG1721'  'MEG1722+1723'
                'MEG1732'  'MEG1733'  'MEG1731'  'MEG1732+1733'
                'MEG1742'  'MEG1743'  'MEG1741'  'MEG1742+1743'
                'MEG1812'  'MEG1813'  'MEG1811'  'MEG1812+1813'
                'MEG1822'  'MEG1823'  'MEG1821'  'MEG1822+1823'
                'MEG1832'  'MEG1833'  'MEG1831'  'MEG1832+1833'
                'MEG1842'  'MEG1843'  'MEG1841'  'MEG1842+1843'
                'MEG1912'  'MEG1913'  'MEG1911'  'MEG1912+1913'
                'MEG1922'  'MEG1923'  'MEG1921'  'MEG1922+1923'
                'MEG1932'  'MEG1933'  'MEG1931'  'MEG1932+1933'
                'MEG1942'  'MEG1943'  'MEG1941'  'MEG1942+1943'
                'MEG2012'  'MEG2013'  'MEG2011'  'MEG2012+2013'
                'MEG2022'  'MEG2023'  'MEG2021'  'MEG2022+2023'
                'MEG2032'  'MEG2033'  'MEG2031'  'MEG2032+2033'
                'MEG2042'  'MEG2043'  'MEG2041'  'MEG2042+2043'
                'MEG2112'  'MEG2113'  'MEG2111'  'MEG2112+2113'
                'MEG2122'  'MEG2123'  'MEG2121'  'MEG2122+2123'
                'MEG2132'  'MEG2133'  'MEG2131'  'MEG2132+2133'
                'MEG2142'  'MEG2143'  'MEG2141'  'MEG2142+2143'
                'MEG2212'  'MEG2213'  'MEG2211'  'MEG2212+2213'
                'MEG2222'  'MEG2223'  'MEG2221'  'MEG2222+2223'
                'MEG2232'  'MEG2233'  'MEG2231'  'MEG2232+2233'
                'MEG2242'  'MEG2243'  'MEG2241'  'MEG2242+2243'
                'MEG2312'  'MEG2313'  'MEG2311'  'MEG2312+2313'
                'MEG2322'  'MEG2323'  'MEG2321'  'MEG2322+2323'
                'MEG2332'  'MEG2333'  'MEG2331'  'MEG2332+2333'
                'MEG2342'  'MEG2343'  'MEG2341'  'MEG2342+2343'
                'MEG2412'  'MEG2413'  'MEG2411'  'MEG2412+2413'
                'MEG2422'  'MEG2423'  'MEG2421'  'MEG2422+2423'
                'MEG2432'  'MEG2433'  'MEG2431'  'MEG2432+2433'
                'MEG2442'  'MEG2443'  'MEG2441'  'MEG2442+2443'
                'MEG2512'  'MEG2513'  'MEG2511'  'MEG2512+2513'
                'MEG2522'  'MEG2523'  'MEG2521'  'MEG2522+2523'
                'MEG2532'  'MEG2533'  'MEG2531'  'MEG2532+2533'
                'MEG2542'  'MEG2543'  'MEG2541'  'MEG2542+2543'
                'MEG2612'  'MEG2613'  'MEG2611'  'MEG2612+2613'
                'MEG2622'  'MEG2623'  'MEG2621'  'MEG2622+2623'
                'MEG2632'  'MEG2633'  'MEG2631'  'MEG2632+2633'
                'MEG2642'  'MEG2643'  'MEG2641'  'MEG2642+2643'
                };
            neuromag306_combined = label(:,4);
            neuromag306alt_combined = label(:,4);
            label = label(:,1:3);

        case 'eeg1020'
            label = {
                'Fp1'
                'Fpz'
                'Fp2'
                'F7'
                'F3'
                'Fz'
                'F4'
                'F8'
                'T7'
                'C3'
                'Cz'
                'C4'
                'T8'
                'P7'
                'P3'
                'Pz'
                'P4'
                'P8'
                'O1'
                'Oz'
                'O2'};

            % Add also reference and some alternative labels that might be used
            label = cat(1, label, {'A1' 'A2' 'M1' 'M2' 'T3' 'T4' 'T5' 'T6'}');

        case 'eeg1010'
            label = {
                'Fp1'
                'Fpz'
                'Fp2'
                'AF9'
                'AF7'
                'AF5'
                'AF3'
                'AF1'
                'AFz'
                'AF2'
                'AF4'
                'AF6'
                'AF8'
                'AF10'
                'F9'
                'F7'
                'F5'
                'F3'
                'F1'
                'Fz'
                'F2'
                'F4'
                'F6'
                'F8'
                'F10'
                'FT9'
                'FT7'
                'FC5'
                'FC3'
                'FC1'
                'FCz'
                'FC2'
                'FC4'
                'FC6'
                'FT8'
                'FT10'
                'T9'
                'T7'
                'C5'
                'C3'
                'C1'
                'Cz'
                'C2'
                'C4'
                'C6'
                'T8'
                'T10'
                'TP9'
                'TP7'
                'CP5'
                'CP3'
                'CP1'
                'CPz'
                'CP2'
                'CP4'
                'CP6'
                'TP8'
                'TP10'
                'P9'
                'P7'
                'P5'
                'P3'
                'P1'
                'Pz'
                'P2'
                'P4'
                'P6'
                'P8'
                'P10'
                'PO9'
                'PO7'
                'PO5'
                'PO3'
                'PO1'
                'POz'
                'PO2'
                'PO4'
                'PO6'
                'PO8'
                'PO10'
                'O1'
                'Oz'
                'O2'
                'I1'
                'Iz'
                'I2'
                };

            % Add also reference and some alternative labels that might be used
            label = cat(1, label, {'A1' 'A2' 'M1' 'M2' 'T3' 'T4' 'T5' 'T6'}');

        case 'eeg1005'
            label = {
                'Fp1'
                'Fpz'
                'Fp2'
                'AF9'
                'AF7'
                'AF5'
                'AF3'
                'AF1'
                'AFz'
                'AF2'
                'AF4'
                'AF6'
                'AF8'
                'AF10'
                'F9'
                'F7'
                'F5'
                'F3'
                'F1'
                'Fz'
                'F2'
                'F4'
                'F6'
                'F8'
                'F10'
                'FT9'
                'FT7'
                'FC5'
                'FC3'
                'FC1'
                'FCz'
                'FC2'
                'FC4'
                'FC6'
                'FT8'
                'FT10'
                'T9'
                'T7'
                'C5'
                'C3'
                'C1'
                'Cz'
                'C2'
                'C4'
                'C6'
                'T8'
                'T10'
                'TP9'
                'TP7'
                'CP5'
                'CP3'
                'CP1'
                'CPz'
                'CP2'
                'CP4'
                'CP6'
                'TP8'
                'TP10'
                'P9'
                'P7'
                'P5'
                'P3'
                'P1'
                'Pz'
                'P2'
                'P4'
                'P6'
                'P8'
                'P10'
                'PO9'
                'PO7'
                'PO5'
                'PO3'
                'PO1'
                'POz'
                'PO2'
                'PO4'
                'PO6'
                'PO8'
                'PO10'
                'O1'
                'Oz'
                'O2'
                'I1'
                'Iz'
                'I2'
                'AFp9h'
                'AFp7h'
                'AFp5h'
                'AFp3h'
                'AFp1h'
                'AFp2h'
                'AFp4h'
                'AFp6h'
                'AFp8h'
                'AFp10h'
                'AFF9h'
                'AFF7h'
                'AFF5h'
                'AFF3h'
                'AFF1h'
                'AFF2h'
                'AFF4h'
                'AFF6h'
                'AFF8h'
                'AFF10h'
                'FFT9h'
                'FFT7h'
                'FFC5h'
                'FFC3h'
                'FFC1h'
                'FFC2h'
                'FFC4h'
                'FFC6h'
                'FFT8h'
                'FFT10h'
                'FTT9h'
                'FTT7h'
                'FCC5h'
                'FCC3h'
                'FCC1h'
                'FCC2h'
                'FCC4h'
                'FCC6h'
                'FTT8h'
                'FTT10h'
                'TTP9h'
                'TTP7h'
                'CCP5h'
                'CCP3h'
                'CCP1h'
                'CCP2h'
                'CCP4h'
                'CCP6h'
                'TTP8h'
                'TTP10h'
                'TPP9h'
                'TPP7h'
                'CPP5h'
                'CPP3h'
                'CPP1h'
                'CPP2h'
                'CPP4h'
                'CPP6h'
                'TPP8h'
                'TPP10h'
                'PPO9h'
                'PPO7h'
                'PPO5h'
                'PPO3h'
                'PPO1h'
                'PPO2h'
                'PPO4h'
                'PPO6h'
                'PPO8h'
                'PPO10h'
                'POO9h'
                'POO7h'
                'POO5h'
                'POO3h'
                'POO1h'
                'POO2h'
                'POO4h'
                'POO6h'
                'POO8h'
                'POO10h'
                'OI1h'
                'OI2h'
                'Fp1h'
                'Fp2h'
                'AF9h'
                'AF7h'
                'AF5h'
                'AF3h'
                'AF1h'
                'AF2h'
                'AF4h'
                'AF6h'
                'AF8h'
                'AF10h'
                'F9h'
                'F7h'
                'F5h'
                'F3h'
                'F1h'
                'F2h'
                'F4h'
                'F6h'
                'F8h'
                'F10h'
                'FT9h'
                'FT7h'
                'FC5h'
                'FC3h'
                'FC1h'
                'FC2h'
                'FC4h'
                'FC6h'
                'FT8h'
                'FT10h'
                'T9h'
                'T7h'
                'C5h'
                'C3h'
                'C1h'
                'C2h'
                'C4h'
                'C6h'
                'T8h'
                'T10h'
                'TP9h'
                'TP7h'
                'CP5h'
                'CP3h'
                'CP1h'
                'CP2h'
                'CP4h'
                'CP6h'
                'TP8h'
                'TP10h'
                'P9h'
                'P7h'
                'P5h'
                'P3h'
                'P1h'
                'P2h'
                'P4h'
                'P6h'
                'P8h'
                'P10h'
                'PO9h'
                'PO7h'
                'PO5h'
                'PO3h'
                'PO1h'
                'PO2h'
                'PO4h'
                'PO6h'
                'PO8h'
                'PO10h'
                'O1h'
                'O2h'
                'I1h'
                'I2h'
                'AFp9'
                'AFp7'
                'AFp5'
                'AFp3'
                'AFp1'
                'AFpz'
                'AFp2'
                'AFp4'
                'AFp6'
                'AFp8'
                'AFp10'
                'AFF9'
                'AFF7'
                'AFF5'
                'AFF3'
                'AFF1'
                'AFFz'
                'AFF2'
                'AFF4'
                'AFF6'
                'AFF8'
                'AFF10'
                'FFT9'
                'FFT7'
                'FFC5'
                'FFC3'
                'FFC1'
                'FFCz'
                'FFC2'
                'FFC4'
                'FFC6'
                'FFT8'
                'FFT10'
                'FTT9'
                'FTT7'
                'FCC5'
                'FCC3'
                'FCC1'
                'FCCz'
                'FCC2'
                'FCC4'
                'FCC6'
                'FTT8'
                'FTT10'
                'TTP9'
                'TTP7'
                'CCP5'
                'CCP3'
                'CCP1'
                'CCPz'
                'CCP2'
                'CCP4'
                'CCP6'
                'TTP8'
                'TTP10'
                'TPP9'
                'TPP7'
                'CPP5'
                'CPP3'
                'CPP1'
                'CPPz'
                'CPP2'
                'CPP4'
                'CPP6'
                'TPP8'
                'TPP10'
                'PPO9'
                'PPO7'
                'PPO5'
                'PPO3'
                'PPO1'
                'PPOz'
                'PPO2'
                'PPO4'
                'PPO6'
                'PPO8'
                'PPO10'
                'POO9'
                'POO7'
                'POO5'
                'POO3'
                'POO1'
                'POOz'
                'POO2'
                'POO4'
                'POO6'
                'POO8'
                'POO10'
                'OI1'
                'OIz'
                'OI2'
                };

            % Add also reference and some alternative labels that might be used
            label = cat(1, label, {'A1' 'A2' 'M1' 'M2' 'T3' 'T4' 'T5' 'T6'}');

        case 'ext1020'
            % start with the eeg1005 list
            label = {
                'Fp1'
                'Fpz'
                'Fp2'
                'AF9'
                'AF7'
                'AF5'
                'AF3'
                'AF1'
                'AFz'
                'AF2'
                'AF4'
                'AF6'
                'AF8'
                'AF10'
                'F9'
                'F7'
                'F5'
                'F3'
                'F1'
                'Fz'
                'F2'
                'F4'
                'F6'
                'F8'
                'F10'
                'FT9'
                'FT7'
                'FC5'
                'FC3'
                'FC1'
                'FCz'
                'FC2'
                'FC4'
                'FC6'
                'FT8'
                'FT10'
                'T9'
                'T7'
                'C5'
                'C3'
                'C1'
                'Cz'
                'C2'
                'C4'
                'C6'
                'T8'
                'T10'
                'TP9'
                'TP7'
                'CP5'
                'CP3'
                'CP1'
                'CPz'
                'CP2'
                'CP4'
                'CP6'
                'TP8'
                'TP10'
                'P9'
                'P7'
                'P5'
                'P3'
                'P1'
                'Pz'
                'P2'
                'P4'
                'P6'
                'P8'
                'P10'
                'PO9'
                'PO7'
                'PO5'
                'PO3'
                'PO1'
                'POz'
                'PO2'
                'PO4'
                'PO6'
                'PO8'
                'PO10'
                'O1'
                'Oz'
                'O2'
                'I1'
                'Iz'
                'I2'
                'AFp9h'
                'AFp7h'
                'AFp5h'
                'AFp3h'
                'AFp1h'
                'AFp2h'
                'AFp4h'
                'AFp6h'
                'AFp8h'
                'AFp10h'
                'AFF9h'
                'AFF7h'
                'AFF5h'
                'AFF3h'
                'AFF1h'
                'AFF2h'
                'AFF4h'
                'AFF6h'
                'AFF8h'
                'AFF10h'
                'FFT9h'
                'FFT7h'
                'FFC5h'
                'FFC3h'
                'FFC1h'
                'FFC2h'
                'FFC4h'
                'FFC6h'
                'FFT8h'
                'FFT10h'
                'FTT9h'
                'FTT7h'
                'FCC5h'
                'FCC3h'
                'FCC1h'
                'FCC2h'
                'FCC4h'
                'FCC6h'
                'FTT8h'
                'FTT10h'
                'TTP9h'
                'TTP7h'
                'CCP5h'
                'CCP3h'
                'CCP1h'
                'CCP2h'
                'CCP4h'
                'CCP6h'
                'TTP8h'
                'TTP10h'
                'TPP9h'
                'TPP7h'
                'CPP5h'
                'CPP3h'
                'CPP1h'
                'CPP2h'
                'CPP4h'
                'CPP6h'
                'TPP8h'
                'TPP10h'
                'PPO9h'
                'PPO7h'
                'PPO5h'
                'PPO3h'
                'PPO1h'
                'PPO2h'
                'PPO4h'
                'PPO6h'
                'PPO8h'
                'PPO10h'
                'POO9h'
                'POO7h'
                'POO5h'
                'POO3h'
                'POO1h'
                'POO2h'
                'POO4h'
                'POO6h'
                'POO8h'
                'POO10h'
                'OI1h'
                'OI2h'
                'Fp1h'
                'Fp2h'
                'AF9h'
                'AF7h'
                'AF5h'
                'AF3h'
                'AF1h'
                'AF2h'
                'AF4h'
                'AF6h'
                'AF8h'
                'AF10h'
                'F9h'
                'F7h'
                'F5h'
                'F3h'
                'F1h'
                'F2h'
                'F4h'
                'F6h'
                'F8h'
                'F10h'
                'FT9h'
                'FT7h'
                'FC5h'
                'FC3h'
                'FC1h'
                'FC2h'
                'FC4h'
                'FC6h'
                'FT8h'
                'FT10h'
                'T9h'
                'T7h'
                'C5h'
                'C3h'
                'C1h'
                'C2h'
                'C4h'
                'C6h'
                'T8h'
                'T10h'
                'TP9h'
                'TP7h'
                'CP5h'
                'CP3h'
                'CP1h'
                'CP2h'
                'CP4h'
                'CP6h'
                'TP8h'
                'TP10h'
                'P9h'
                'P7h'
                'P5h'
                'P3h'
                'P1h'
                'P2h'
                'P4h'
                'P6h'
                'P8h'
                'P10h'
                'PO9h'
                'PO7h'
                'PO5h'
                'PO3h'
                'PO1h'
                'PO2h'
                'PO4h'
                'PO6h'
                'PO8h'
                'PO10h'
                'O1h'
                'O2h'
                'I1h'
                'I2h'
                'AFp9'
                'AFp7'
                'AFp5'
                'AFp3'
                'AFp1'
                'AFpz'
                'AFp2'
                'AFp4'
                'AFp6'
                'AFp8'
                'AFp10'
                'AFF9'
                'AFF7'
                'AFF5'
                'AFF3'
                'AFF1'
                'AFFz'
                'AFF2'
                'AFF4'
                'AFF6'
                'AFF8'
                'AFF10'
                'FFT9'
                'FFT7'
                'FFC5'
                'FFC3'
                'FFC1'
                'FFCz'
                'FFC2'
                'FFC4'
                'FFC6'
                'FFT8'
                'FFT10'
                'FTT9'
                'FTT7'
                'FCC5'
                'FCC3'
                'FCC1'
                'FCCz'
                'FCC2'
                'FCC4'
                'FCC6'
                'FTT8'
                'FTT10'
                'TTP9'
                'TTP7'
                'CCP5'
                'CCP3'
                'CCP1'
                'CCPz'
                'CCP2'
                'CCP4'
                'CCP6'
                'TTP8'
                'TTP10'
                'TPP9'
                'TPP7'
                'CPP5'
                'CPP3'
                'CPP1'
                'CPPz'
                'CPP2'
                'CPP4'
                'CPP6'
                'TPP8'
                'TPP10'
                'PPO9'
                'PPO7'
                'PPO5'
                'PPO3'
                'PPO1'
                'PPOz'
                'PPO2'
                'PPO4'
                'PPO6'
                'PPO8'
                'PPO10'
                'POO9'
                'POO7'
                'POO5'
                'POO3'
                'POO1'
                'POOz'
                'POO2'
                'POO4'
                'POO6'
                'POO8'
                'POO10'
                'OI1'
                'OIz'
                'OI2'
                };

            % Add also reference and some alternative labels that might be used
            label = cat(1, label, {'A1' 'A2' 'M1' 'M2' 'T3' 'T4' 'T5' 'T6'}');

            % This is to account for all variants of case in 1020 systems
            label = unique(cat(1, label, upper(label), lower(label)));

        case 'biosemi64'
            label = {
                'A1'
                'A2'
                'A3'
                'A4'
                'A5'
                'A6'
                'A7'
                'A8'
                'A9'
                'A10'
                'A11'
                'A12'
                'A13'
                'A14'
                'A15'
                'A16'
                'A17'
                'A18'
                'A19'
                'A20'
                'A21'
                'A22'
                'A23'
                'A24'
                'A25'
                'A26'
                'A27'
                'A28'
                'A29'
                'A30'
                'A31'
                'A32'
                'B1'
                'B2'
                'B3'
                'B4'
                'B5'
                'B6'
                'B7'
                'B8'
                'B9'
                'B10'
                'B11'
                'B12'
                'B13'
                'B14'
                'B15'
                'B16'
                'B17'
                'B18'
                'B19'
                'B20'
                'B21'
                'B22'
                'B23'
                'B24'
                'B25'
                'B26'
                'B27'
                'B28'
                'B29'
                'B30'
                'B31'
                'B32'
                };

        case 'biosemi128'
            label = {
                'A1'
                'A2'
                'A3'
                'A4'
                'A5'
                'A6'
                'A7'
                'A8'
                'A9'
                'A10'
                'A11'
                'A12'
                'A13'
                'A14'
                'A15'
                'A16'
                'A17'
                'A18'
                'A19'
                'A20'
                'A21'
                'A22'
                'A23'
                'A24'
                'A25'
                'A26'
                'A27'
                'A28'
                'A29'
                'A30'
                'A31'
                'A32'
                'B1'
                'B2'
                'B3'
                'B4'
                'B5'
                'B6'
                'B7'
                'B8'
                'B9'
                'B10'
                'B11'
                'B12'
                'B13'
                'B14'
                'B15'
                'B16'
                'B17'
                'B18'
                'B19'
                'B20'
                'B21'
                'B22'
                'B23'
                'B24'
                'B25'
                'B26'
                'B27'
                'B28'
                'B29'
                'B30'
                'B31'
                'B32'
                'C1'
                'C2'
                'C3'
                'C4'
                'C5'
                'C6'
                'C7'
                'C8'
                'C9'
                'C10'
                'C11'
                'C12'
                'C13'
                'C14'
                'C15'
                'C16'
                'C17'
                'C18'
                'C19'
                'C20'
                'C21'
                'C22'
                'C23'
                'C24'
                'C25'
                'C26'
                'C27'
                'C28'
                'C29'
                'C30'
                'C31'
                'C32'
                'D1'
                'D2'
                'D3'
                'D4'
                'D5'
                'D6'
                'D7'
                'D8'
                'D9'
                'D10'
                'D11'
                'D12'
                'D13'
                'D14'
                'D15'
                'D16'
                'D17'
                'D18'
                'D19'
                'D20'
                'D21'
                'D22'
                'D23'
                'D24'
                'D25'
                'D26'
                'D27'
                'D28'
                'D29'
                'D30'
                'D31'
                'D32'
                };

        case 'biosemi256'
            label = {
                'A1'
                'A2'
                'A3'
                'A4'
                'A5'
                'A6'
                'A7'
                'A8'
                'A9'
                'A10'
                'A11'
                'A12'
                'A13'
                'A14'
                'A15'
                'A16'
                'A17'
                'A18'
                'A19'
                'A20'
                'A21'
                'A22'
                'A23'
                'A24'
                'A25'
                'A26'
                'A27'
                'A28'
                'A29'
                'A30'
                'A31'
                'A32'
                'B1'
                'B2'
                'B3'
                'B4'
                'B5'
                'B6'
                'B7'
                'B8'
                'B9'
                'B10'
                'B11'
                'B12'
                'B13'
                'B14'
                'B15'
                'B16'
                'B17'
                'B18'
                'B19'
                'B20'
                'B21'
                'B22'
                'B23'
                'B24'
                'B25'
                'B26'
                'B27'
                'B28'
                'B29'
                'B30'
                'B31'
                'B32'
                'C1'
                'C2'
                'C3'
                'C4'
                'C5'
                'C6'
                'C7'
                'C8'
                'C9'
                'C10'
                'C11'
                'C12'
                'C13'
                'C14'
                'C15'
                'C16'
                'C17'
                'C18'
                'C19'
                'C20'
                'C21'
                'C22'
                'C23'
                'C24'
                'C25'
                'C26'
                'C27'
                'C28'
                'C29'
                'C30'
                'C31'
                'C32'
                'D1'
                'D2'
                'D3'
                'D4'
                'D5'
                'D6'
                'D7'
                'D8'
                'D9'
                'D10'
                'D11'
                'D12'
                'D13'
                'D14'
                'D15'
                'D16'
                'D17'
                'D18'
                'D19'
                'D20'
                'D21'
                'D22'
                'D23'
                'D24'
                'D25'
                'D26'
                'D27'
                'D28'
                'D29'
                'D30'
                'D31'
                'D32'
                'E1'
                'E2'
                'E3'
                'E4'
                'E5'
                'E6'
                'E7'
                'E8'
                'E9'
                'E10'
                'E11'
                'E12'
                'E13'
                'E14'
                'E15'
                'E16'
                'E17'
                'E18'
                'E19'
                'E20'
                'E21'
                'E22'
                'E23'
                'E24'
                'E25'
                'E26'
                'E27'
                'E28'
                'E29'
                'E30'
                'E31'
                'E32'
                'F1'
                'F2'
                'F3'
                'F4'
                'F5'
                'F6'
                'F7'
                'F8'
                'F9'
                'F10'
                'F11'
                'F12'
                'F13'
                'F14'
                'F15'
                'F16'
                'F17'
                'F18'
                'F19'
                'F20'
                'F21'
                'F22'
                'F23'
                'F24'
                'F25'
                'F26'
                'F27'
                'F28'
                'F29'
                'F30'
                'F31'
                'F32'
                'G1'
                'G2'
                'G3'
                'G4'
                'G5'
                'G6'
                'G7'
                'G8'
                'G9'
                'G10'
                'G11'
                'G12'
                'G13'
                'G14'
                'G15'
                'G16'
                'G17'
                'G18'
                'G19'
                'G20'
                'G21'
                'G22'
                'G23'
                'G24'
                'G25'
                'G26'
                'G27'
                'G28'
                'G29'
                'G30'
                'G31'
                'G32'
                'H1'
                'H2'
                'H3'
                'H4'
                'H5'
                'H6'
                'H7'
                'H8'
                'H9'
                'H10'
                'H11'
                'H12'
                'H13'
                'H14'
                'H15'
                'H16'
                'H17'
                'H18'
                'H19'
                'H20'
                'H21'
                'H22'
                'H23'
                'H24'
                'H25'
                'H26'
                'H27'
                'H28'
                'H29'
                'H30'
                'H31'
                'H32'
                };

        case 'egi32'
            % this should be  uppercase for consistency with ft_read_header
            label = cell(33, 1);
            for i = 1:33
                label{i} = sprintf('E%d', i);
            end
            % there might also be a reference channel, but its name is inconsistent
            % it might be Cz, REF, VREF or 'vertex reference'

        case 'egi64'
            % this should be  uppercase for consistency with ft_read_header
            label = cell(65, 1);
            for i = 1:65
                label{i} = sprintf('E%d', i);
            end
            % there might also be a reference channel, but its name is inconsistent
            % it might be Cz, REF, VREF or 'vertex reference'

        case 'egi128'
            % this should be  uppercase for consistency with ft_read_header
            label = cell(129, 1);
            for i = 1:129
                label{i} = sprintf('E%d', i);
            end
            % there might also be a reference channel, but its name is inconsistent
            % it might be Cz, REF, VREF or 'vertex reference'

        case 'egi256'
            % this should be  uppercase for consistency with ft_read_header
            label = cell(257, 1);
            for i = 1:257
                label{i} = sprintf('E%d', i);
            end
            % there might also be a reference channel, but its name is inconsistent
            % it might be Cz, REF, VREF or 'vertex reference'

        case 'itab28'
            label = {
                'MAG_1'
                'MAG_2'
                'MAG_3'
                'MAG_5'
                'MAG_7'
                'MAG_8'
                'MAG_9'
                'MAG_11'
                'MAG_12'
                'MAG_13'
                'MAG_15'
                'MAG_17'
                'MAG_18'
                'MAG_21'
                'MAG_22'
                'MAG_23'
                'MAG_25'
                'MAG_26'
                'MAG_27'
                'MAG_28'
                };

        case 'itab153'
            label = cell(153,1);
            for i=1:153
                % channel names start counting at zero
                label{i} = sprintf('MAG_%03d',  i-1);
            end

        case 'itab153_planar'
            label = cell(153,3);
            for i=1:153
                % channel names start counting at zero
                label{i,1} = sprintf('MAG_%03d_dH', i-1);
                label{i,2} = sprintf('MAG_%03d_dV', i-1);
                label{i,3} = sprintf('MAG_%03d',  i-1);
            end
            itab153_planar_combined = label(:,3);
            label = label(:,1:2);

        case 'yokogawa9'
            % note that this uses MATLAB style 1-offset indexing and not C style 0-offset indexing
            % this should be consistent with: read_yokogawa_header, ft_channelselection, yokogawa2grad
            label = cell(9,1);
            for i=1:9
                label{i} = sprintf('M%03d',  i);
            end

        case 'yokogawa64'
            % note that this uses MATLAB style 1-offset indexing and not C style 0-offset indexing
            % this should be consistent with: read_yokogawa_header, ft_channelselection, yokogawa2grad
            label = cell(64,1);
            for i=1:64
                label{i} = sprintf('AG%03d', i);
            end

        case 'yokogawa64_planar'
            % note that this uses MATLAB style 1-offset indexing and not C style 0-offset indexing
            % this should be consistent with: read_yokogawa_header, ft_channelselection, yokogawa2grad
            label = cell(64,3);
            for i=1:64
                label{i,1} = sprintf('AG%03d_dH', i);
                label{i,2} = sprintf('AG%03d_dV', i);
                label{i,3} = sprintf('AG%03d', i);
            end
            yokogawa64_planar_combined = label(:,3);
            label = label(:,1:2);

        case 'yokogawa160'
            % note that this uses MATLAB style 1-offset indexing and not C style 0-offset indexing
            % this should be consistent with: read_yokogawa_header, ft_channelselection, yokogawa2grad
            label = cell(160,1);
            for i=1:160
                label{i} = sprintf('AG%03d', i);
            end

        case 'yokogawa160_planar'
            % note that this uses MATLAB style 1-offset indexing and not C style 0-offset indexing
            % this should be consistent with: read_yokogawa_header, ft_channelselection, yokogawa2grad
            label = cell(160,2);
            for i=1:160
                label{i,1} = sprintf('AG%03d_dH', i);
                label{i,2} = sprintf('AG%03d_dV', i);
                label{i,3} = sprintf('AG%03d', i);
            end
            yokogawa160_planar_combined = label(:,3);
            label = label(:,1:2);

        case 'yokogawa440'
            % this should be consistent with read_yokogawa_header, with ft_channelselection and with yokogawa2grad
            label = {
                'AG001'
                'AG002'
                'AG003'
                'AG004'
                'AG005'
                'AG006'
                'AG007'
                'AG008'
                'AG009'
                'AG010'
                'AG011'
                'AG012'
                'AG013'
                'AG014'
                'AG015'
                'AG016'
                'AG017'
                'AG018'
                'AG019'
                'AG020'
                'AG021'
                'AG022'
                'AG023'
                'AG024'
                'AG025'
                'AG026'
                'AG027'
                'AG028'
                'AG029'
                'AG030'
                'AG031'
                'AG032'
                'PG033'
                'PG034'
                'PG035'
                'PG036'
                'PG037'
                'PG038'
                'PG039'
                'PG040'
                'PG041'
                'PG042'
                'PG043'
                'PG044'
                'PG045'
                'PG046'
                'PG047'
                'PG048'
                'PG049'
                'PG050'
                'PG051'
                'PG052'
                'PG053'
                'PG054'
                'PG055'
                'PG056'
                'PG057'
                'PG058'
                'PG059'
                'PG060'
                'PG061'
                'PG062'
                'PG063'
                'PG064'
                'AG065'
                'AG066'
                'AG067'
                'AG068'
                'AG069'
                'AG070'
                'AG071'
                'AG072'
                'AG073'
                'AG074'
                'AG075'
                'AG076'
                'AG077'
                'AG078'
                'AG079'
                'AG080'
                'AG081'
                'AG082'
                'AG083'
                'AG084'
                'AG085'
                'AG086'
                'AG087'
                'AG088'
                'AG089'
                'AG090'
                'AG091'
                'AG092'
                'AG093'
                'AG094'
                'AG095'
                'AG096'
                'PG097'
                'PG098'
                'PG099'
                'PG100'
                'PG101'
                'PG102'
                'PG103'
                'PG104'
                'PG105'
                'PG106'
                'PG107'
                'PG108'
                'PG109'
                'PG110'
                'PG111'
                'PG112'
                'PG113'
                'PG114'
                'PG115'
                'PG116'
                'PG117'
                'PG118'
                'PG119'
                'PG120'
                'PG121'
                'AG122'
                'PG123'
                'PG124'
                'PG125'
                'PG126'
                'PG127'
                'PG128'
                'AG129'
                'AG130'
                'AG131'
                'AG132'
                'AG133'
                'AG134'
                'AG135'
                'AG136'
                'AG137'
                'AG138'
                'AG139'
                'AG140'
                'AG141'
                'AG142'
                'AG143'
                'AG144'
                'AG145'
                'AG146'
                'AG147'
                'AG148'
                'AG149'
                'AG150'
                'AG151'
                'AG152'
                'AG153'
                'AG154'
                'AG155'
                'AG156'
                'AG157'
                'AG158'
                'AG159'
                'AG160'
                'AG161'
                'AG162'
                'AG163'
                'AG164'
                'AG165'
                'PG166'
                'PG167'
                'PG168'
                'PG169'
                'PG170'
                'PG171'
                'PG172'
                'PG173'
                'PG174'
                'PG175'
                'PG176'
                'PG177'
                'PG178'
                'PG179'
                'PG180'
                'PG181'
                'PG182'
                'PG183'
                'PG184'
                'PG185'
                'PG186'
                'PG187'
                'PG188'
                'PG189'
                'PG190'
                'PG191'
                'PG192'
                'AG193'
                'AG194'
                'AG195'
                'AG196'
                'AG197'
                'AG198'
                'AG199'
                'AG200'
                'AG201'
                'AG202'
                'AG203'
                'AG204'
                'AG205'
                'AG206'
                'AG207'
                'AG208'
                'AG209'
                'AG210'
                'AG211'
                'AG212'
                'AG213'
                'AG214'
                'AG215'
                'AG216'
                'AG217'
                'AG218'
                'AG219'
                'AG220'
                'AG221'
                'AG222'
                'AG223'
                'AG224'
                'AG225'
                'AG226'
                'AG227'
                'PG228'
                'PG229'
                'PG230'
                'PG231'
                'PG232'
                'PG233'
                'PG234'
                'PG235'
                'PG236'
                'PG237'
                'PG238'
                'PG239'
                'PG240'
                'PG241'
                'PG242'
                'PG243'
                'PG244'
                'PG245'
                'PG246'
                'PG247'
                'PG248'
                'PG249'
                'PG250'
                'PG251'
                'PG252'
                'PG253'
                'PG254'
                'PG255'
                'PG256'
                'AG257'
                'AG258'
                'AG259'
                'AG260'
                'AG261'
                'AG262'
                'AG263'
                'AG264'
                'AG265'
                'AG266'
                'AG267'
                'AG268'
                'AG269'
                'AG270'
                'AG271'
                'AG272'
                'AG273'
                'AG274'
                'AG275'
                'AG276'
                'AG277'
                'AG278'
                'AG279'
                'AG280'
                'AG281'
                'AG282'
                'AG283'
                'AG284'
                'AG285'
                'AG286'
                'AG287'
                'AG288'
                'PG289'
                'PG290'
                'PG291'
                'PG292'
                'PG293'
                'PG294'
                'PG295'
                'PG296'
                'PG297'
                'PG298'
                'PG299'
                'PG300'
                'PG301'
                'PG302'
                'PG303'
                'PG304'
                'PG305'
                'PG306'
                'PG307'
                'PG308'
                'PG309'
                'PG310'
                'PG311'
                'PG312'
                'PG313'
                'PG314'
                'PG315'
                'PG316'
                'PG317'
                'PG318'
                'PG319'
                'PG320'
                'AG321'
                'AG322'
                'AG323'
                'AG324'
                'AG325'
                'AG326'
                'AG327'
                'AG328'
                'AG329'
                'AG330'
                'AG331'
                'AG332'
                'AG333'
                'AG334'
                'AG335'
                'AG336'
                'AG337'
                'AG338'
                'AG339'
                'AG340'
                'AG341'
                'AG342'
                'AG343'
                'AG344'
                'AG345'
                'AG346'
                'AG347'
                'AG348'
                'AG349'
                'AG350'
                'AG351'
                'AG352'
                'PG353'
                'PG354'
                'PG355'
                'PG356'
                'PG357'
                'PG358'
                'PG359'
                'PG360'
                'PG361'
                'PG362'
                'PG363'
                'PG364'
                'PG365'
                'PG366'
                'PG367'
                'PG368'
                'PG369'
                'PG370'
                'PG371'
                'PG372'
                'PG373'
                'PG374'
                'PG375'
                'PG376'
                'PG377'
                'AG378'
                'PG379'
                'PG380'
                'PG381'
                'PG382'
                'PG383'
                'PG384'
                'AG385'
                'AG386'
                'AG387'
                'AG388'
                'AG389'
                'AG390'
                'AG391'
                'AG392'
                'PG393'
                'PG394'
                'PG395'
                'PG396'
                'PG397'
                'PG398'
                'PG399'
                'PG400'
                'RM401'
                'RM402'
                'RM403'
                'RM404'
                'RM405'
                'RM406'
                'RM407'
                'RM408'
                'RM409'
                'RM410'
                'RM411'
                'RM412'
                };

        case 'yokogawa440_planar'
            % this should be consistent with read_yokogawa_header, with
            % ft_channelselection and with yokogawa2grad
            label = {
                'AG001_dH'  'AG001_dV'  'AG001'
                'AG002_dH'  'AG002_dV'  'AG002'
                'AG003_dH'  'AG003_dV'  'AG003'
                'AG004_dH'  'AG004_dV'  'AG004'
                'AG005_dH'  'AG005_dV'  'AG005'
                'AG006_dH'  'AG006_dV'  'AG006'
                'AG007_dH'  'AG007_dV'  'AG007'
                'AG008_dH'  'AG008_dV'  'AG008'
                'AG009_dH'  'AG009_dV'  'AG009'
                'AG010_dH'  'AG010_dV'  'AG010'
                'AG011_dH'  'AG011_dV'  'AG011'
                'AG012_dH'  'AG012_dV'  'AG012'
                'AG013_dH'  'AG013_dV'  'AG013'
                'AG014_dH'  'AG014_dV'  'AG014'
                'AG015_dH'  'AG015_dV'  'AG015'
                'AG016_dH'  'AG016_dV'  'AG016'
                'AG017_dH'  'AG017_dV'  'AG017'
                'AG018_dH'  'AG018_dV'  'AG018'
                'AG019_dH'  'AG019_dV'  'AG019'
                'AG020_dH'  'AG020_dV'  'AG020'
                'AG021_dH'  'AG021_dV'  'AG021'
                'AG022_dH'  'AG022_dV'  'AG022'
                'AG023_dH'  'AG023_dV'  'AG023'
                'AG024_dH'  'AG024_dV'  'AG024'
                'AG025_dH'  'AG025_dV'  'AG025'
                'AG026_dH'  'AG026_dV'  'AG026'
                'AG027_dH'  'AG027_dV'  'AG027'
                'AG028_dH'  'AG028_dV'  'AG028'
                'AG029_dH'  'AG029_dV'  'AG029'
                'AG030_dH'  'AG030_dV'  'AG030'
                'AG031_dH'  'AG031_dV'  'AG031'
                'AG032_dH'  'AG032_dV'  'AG032'
                'AG065_dH'  'AG065_dV'  'AG065'
                'AG066_dH'  'AG066_dV'  'AG066'
                'AG067_dH'  'AG067_dV'  'AG067'
                'AG068_dH'  'AG068_dV'  'AG068'
                'AG069_dH'  'AG069_dV'  'AG069'
                'AG070_dH'  'AG070_dV'  'AG070'
                'AG071_dH'  'AG071_dV'  'AG071'
                'AG072_dH'  'AG072_dV'  'AG072'
                'AG073_dH'  'AG073_dV'  'AG073'
                'AG074_dH'  'AG074_dV'  'AG074'
                'AG075_dH'  'AG075_dV'  'AG075'
                'AG076_dH'  'AG076_dV'  'AG076'
                'AG077_dH'  'AG077_dV'  'AG077'
                'AG078_dH'  'AG078_dV'  'AG078'
                'AG079_dH'  'AG079_dV'  'AG079'
                'AG080_dH'  'AG080_dV'  'AG080'
                'AG081_dH'  'AG081_dV'  'AG081'
                'AG082_dH'  'AG082_dV'  'AG082'
                'AG083_dH'  'AG083_dV'  'AG083'
                'AG084_dH'  'AG084_dV'  'AG084'
                'AG085_dH'  'AG085_dV'  'AG085'
                'AG086_dH'  'AG086_dV'  'AG086'
                'AG087_dH'  'AG087_dV'  'AG087'
                'AG088_dH'  'AG088_dV'  'AG088'
                'AG089_dH'  'AG089_dV'  'AG089'
                'AG090_dH'  'AG090_dV'  'AG090'
                'AG091_dH'  'AG091_dV'  'AG091'
                'AG092_dH'  'AG092_dV'  'AG092'
                'AG093_dH'  'AG093_dV'  'AG093'
                'AG094_dH'  'AG094_dV'  'AG094'
                'AG095_dH'  'AG095_dV'  'AG095'
                'AG096_dH'  'AG096_dV'  'AG096'
                'AG122_dH'  'AG122_dV'  'AG122'
                'AG129_dH'  'AG129_dV'  'AG129'
                'AG130_dH'  'AG130_dV'  'AG130'
                'AG131_dH'  'AG131_dV'  'AG131'
                'AG132_dH'  'AG132_dV'  'AG132'
                'AG133_dH'  'AG133_dV'  'AG133'
                'AG134_dH'  'AG134_dV'  'AG134'
                'AG135_dH'  'AG135_dV'  'AG135'
                'AG136_dH'  'AG136_dV'  'AG136'
                'AG137_dH'  'AG137_dV'  'AG137'
                'AG138_dH'  'AG138_dV'  'AG138'
                'AG139_dH'  'AG139_dV'  'AG139'
                'AG140_dH'  'AG140_dV'  'AG140'
                'AG141_dH'  'AG141_dV'  'AG141'
                'AG142_dH'  'AG142_dV'  'AG142'
                'AG143_dH'  'AG143_dV'  'AG143'
                'AG144_dH'  'AG144_dV'  'AG144'
                'AG145_dH'  'AG145_dV'  'AG145'
                'AG146_dH'  'AG146_dV'  'AG146'
                'AG147_dH'  'AG147_dV'  'AG147'
                'AG148_dH'  'AG148_dV'  'AG148'
                'AG149_dH'  'AG149_dV'  'AG149'
                'AG150_dH'  'AG150_dV'  'AG150'
                'AG151_dH'  'AG151_dV'  'AG151'
                'AG152_dH'  'AG152_dV'  'AG152'
                'AG153_dH'  'AG153_dV'  'AG153'
                'AG154_dH'  'AG154_dV'  'AG154'
                'AG155_dH'  'AG155_dV'  'AG155'
                'AG156_dH'  'AG156_dV'  'AG156'
                'AG157_dH'  'AG157_dV'  'AG157'
                'AG158_dH'  'AG158_dV'  'AG158'
                'AG159_dH'  'AG159_dV'  'AG159'
                'AG160_dH'  'AG160_dV'  'AG160'
                'AG161_dH'  'AG161_dV'  'AG161'
                'AG162_dH'  'AG162_dV'  'AG162'
                'AG163_dH'  'AG163_dV'  'AG163'
                'AG164_dH'  'AG164_dV'  'AG164'
                'AG165_dH'  'AG165_dV'  'AG165'
                'AG193_dH'  'AG193_dV'  'AG193'
                'AG194_dH'  'AG194_dV'  'AG194'
                'AG195_dH'  'AG195_dV'  'AG195'
                'AG196_dH'  'AG196_dV'  'AG196'
                'AG197_dH'  'AG197_dV'  'AG197'
                'AG198_dH'  'AG198_dV'  'AG198'
                'AG199_dH'  'AG199_dV'  'AG199'
                'AG200_dH'  'AG200_dV'  'AG200'
                'AG201_dH'  'AG201_dV'  'AG201'
                'AG202_dH'  'AG202_dV'  'AG202'
                'AG203_dH'  'AG203_dV'  'AG203'
                'AG204_dH'  'AG204_dV'  'AG204'
                'AG205_dH'  'AG205_dV'  'AG205'
                'AG206_dH'  'AG206_dV'  'AG206'
                'AG207_dH'  'AG207_dV'  'AG207'
                'AG208_dH'  'AG208_dV'  'AG208'
                'AG209_dH'  'AG209_dV'  'AG209'
                'AG210_dH'  'AG210_dV'  'AG210'
                'AG211_dH'  'AG211_dV'  'AG211'
                'AG212_dH'  'AG212_dV'  'AG212'
                'AG213_dH'  'AG213_dV'  'AG213'
                'AG214_dH'  'AG214_dV'  'AG214'
                'AG215_dH'  'AG215_dV'  'AG215'
                'AG216_dH'  'AG216_dV'  'AG216'
                'AG217_dH'  'AG217_dV'  'AG217'
                'AG218_dH'  'AG218_dV'  'AG218'
                'AG219_dH'  'AG219_dV'  'AG219'
                'AG220_dH'  'AG220_dV'  'AG220'
                'AG221_dH'  'AG221_dV'  'AG221'
                'AG222_dH'  'AG222_dV'  'AG222'
                'AG223_dH'  'AG223_dV'  'AG223'
                'AG224_dH'  'AG224_dV'  'AG224'
                'AG225_dH'  'AG225_dV'  'AG225'
                'AG226_dH'  'AG226_dV'  'AG226'
                'AG227_dH'  'AG227_dV'  'AG227'
                'AG257_dH'  'AG257_dV'  'AG257'
                'AG258_dH'  'AG258_dV'  'AG258'
                'AG259_dH'  'AG259_dV'  'AG259'
                'AG260_dH'  'AG260_dV'  'AG260'
                'AG261_dH'  'AG261_dV'  'AG261'
                'AG262_dH'  'AG262_dV'  'AG262'
                'AG263_dH'  'AG263_dV'  'AG263'
                'AG264_dH'  'AG264_dV'  'AG264'
                'AG265_dH'  'AG265_dV'  'AG265'
                'AG266_dH'  'AG266_dV'  'AG266'
                'AG267_dH'  'AG267_dV'  'AG267'
                'AG268_dH'  'AG268_dV'  'AG268'
                'AG269_dH'  'AG269_dV'  'AG269'
                'AG270_dH'  'AG270_dV'  'AG270'
                'AG271_dH'  'AG271_dV'  'AG271'
                'AG272_dH'  'AG272_dV'  'AG272'
                'AG273_dH'  'AG273_dV'  'AG273'
                'AG274_dH'  'AG274_dV'  'AG274'
                'AG275_dH'  'AG275_dV'  'AG275'
                'AG276_dH'  'AG276_dV'  'AG276'
                'AG277_dH'  'AG277_dV'  'AG277'
                'AG278_dH'  'AG278_dV'  'AG278'
                'AG279_dH'  'AG279_dV'  'AG279'
                'AG280_dH'  'AG280_dV'  'AG280'
                'AG281_dH'  'AG281_dV'  'AG281'
                'AG282_dH'  'AG282_dV'  'AG282'
                'AG283_dH'  'AG283_dV'  'AG283'
                'AG284_dH'  'AG284_dV'  'AG284'
                'AG285_dH'  'AG285_dV'  'AG285'
                'AG286_dH'  'AG286_dV'  'AG286'
                'AG287_dH'  'AG287_dV'  'AG287'
                'AG288_dH'  'AG288_dV'  'AG288'
                'AG321_dH'  'AG321_dV'  'AG321'
                'AG322_dH'  'AG322_dV'  'AG322'
                'AG323_dH'  'AG323_dV'  'AG323'
                'AG324_dH'  'AG324_dV'  'AG324'
                'AG325_dH'  'AG325_dV'  'AG325'
                'AG326_dH'  'AG326_dV'  'AG326'
                'AG327_dH'  'AG327_dV'  'AG327'
                'AG328_dH'  'AG328_dV'  'AG328'
                'AG329_dH'  'AG329_dV'  'AG329'
                'AG330_dH'  'AG330_dV'  'AG330'
                'AG331_dH'  'AG331_dV'  'AG331'
                'AG332_dH'  'AG332_dV'  'AG332'
                'AG333_dH'  'AG333_dV'  'AG333'
                'AG334_dH'  'AG334_dV'  'AG334'
                'AG335_dH'  'AG335_dV'  'AG335'
                'AG336_dH'  'AG336_dV'  'AG336'
                'AG337_dH'  'AG337_dV'  'AG337'
                'AG338_dH'  'AG338_dV'  'AG338'
                'AG339_dH'  'AG339_dV'  'AG339'
                'AG340_dH'  'AG340_dV'  'AG340'
                'AG341_dH'  'AG341_dV'  'AG341'
                'AG342_dH'  'AG342_dV'  'AG342'
                'AG343_dH'  'AG343_dV'  'AG343'
                'AG344_dH'  'AG344_dV'  'AG344'
                'AG345_dH'  'AG345_dV'  'AG345'
                'AG346_dH'  'AG346_dV'  'AG346'
                'AG347_dH'  'AG347_dV'  'AG347'
                'AG348_dH'  'AG348_dV'  'AG348'
                'AG349_dH'  'AG349_dV'  'AG349'
                'AG350_dH'  'AG350_dV'  'AG350'
                'AG351_dH'  'AG351_dV'  'AG351'
                'AG352_dH'  'AG352_dV'  'AG352'
                'AG378_dH'  'AG378_dV'  'AG378'
                'AG385_dH'  'AG385_dV'  'AG385'
                'AG386_dH'  'AG386_dV'  'AG386'
                'AG387_dH'  'AG387_dV'  'AG387'
                'AG388_dH'  'AG388_dV'  'AG388'
                'AG389_dH'  'AG389_dV'  'AG389'
                'AG390_dH'  'AG390_dV'  'AG390'
                'AG391_dH'  'AG391_dV'  'AG391'
                'AG392_dH'  'AG392_dV'  'AG392'
                };
            yokogawa440_planar_combined = label(:,3);
            label = label(:,1:2);

        case {'eeg' 'electrode'}
            % there is no default set of electrode labels for all possible EEG systems
            % but nevertheless the requested input type should not result in an error
            label = {};

        otherwise
            error('the requested sensor type "%s" is not supported', type);

    end % switch

    % remember this set of labels to speed up subsequent function calls
    eval(sprintf('%s = label;', type));
    clear label

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch output
    case 'normal'
        % return labels as 2*Nx1 cell-array for planar systems or 3*Nx1 for neuromag306
        % return labels as   Nx1 cell-array for non-planar systems
        label = eval(type);

    case 'planarcombined'
        % return labels as Nx3 cell-array for the planar channels, 3rd column contains the combination
        planar    = eval(type);
        combined  = eval([type '_combined']);
        label     = [planar(:,1:2) combined]; % magnetometers are in the 3rd column for neuromag306

    otherwise
        error('unsupported output "%s"', output);

end

function vol = ea_ft_headmodel_simbio(geom, varargin)

% FT_HEADMODEL_SIMBIO creates a volume conduction model of the head
% using the finite element method (FEM) for EEG. This function takes
% as input a volumetric mesh (hexahedral or tetrahedral) and
% returns as output a volume conduction model which can be used to
% compute leadfields.
%
% This implements
%       ...
%
% Use as
%   vol = ft_headmodel_simbio(geom,'conductivity', conductivities, ...)
%
% The geom is given as a volumetric mesh, using ft_datatype_parcellation
%   geom.pos = vertex positions
%   geom.tet/geom.hex = list of volume elements
%   geom.tissue = tissue assignment for elements
%   geom.tissuelabel = labels correspondig to tissues
%
% Required input arguments should be specified in key-value pairs and have
% to include
%   conductivity   = vector containing tissue conductivities using ordered
%                    corresponding to geom.tissuelabel
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% To run this on Windows the following packages are necessary:
%
% Microsoft Visual C++ 2008 Redistributable
%
% Intel Visual Fortran Redistributables
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

% $Id: ft_headmodel_simbio.m 8445 2013-09-03 10:01:42Z johvor $


% get the optional arguments
conductivity    = ea_ft_getopt(varargin, 'conductivity');

% start with an empty volume conductor
geom = ea_ft_datatype_parcellation(geom);
vol = [];
if isfield(geom,'pos')
    vol.pos = geom.pos;
else
    error('Vertex field is required!')
end

if isfield(geom,'tet')
    vol.tet = geom.tet;
elseif isfield(geom,'hex')
    vol.hex = geom.hex;
else
    error('Connectivity information is required!')
end

if isfield(geom,'tissue')
    vol.tissue = geom.tissue;
else
    error('No element indices declared!')
end

if isempty(conductivity)
    error('No conductivity information!')
end

if length(conductivity) >= length(unique(vol.tissue))
    vol.cond = conductivity;
else
    keyboard
    error('Wrong conductivity information!')
end

if ~isfield(geom,'tissuelabel')
    numlabels = size(unique(geom.tissue),1);
    vol.tissuelabel = {};
    ulabel = unique(labels);
    for i = 1:numlabels
        vol.tissuelabel{i} = num2str(ulabel(i));
    end
else
    vol.tissuelabel = geom.tissuelabel;
end

vol.stiff = ea_sb_calc_stiff(vol);
vol.type = 'simbio';


function parcellation = ea_ft_datatype_parcellation(parcellation, varargin)

% FT_DATATYPE_PARCELLATION describes the FieldTrip MATLAB structure for parcellated
% cortex-based data and atlases. A parcellation can either be indexed or probabilistic
% (see below).
%
% A parcellation describes the tissue types for each of the surface elements.
% Parcellations are often, but not always labeled. A parcellatoin can be used to
% estimate the activity from MEG data in a known region of interest. A surface-based
% atlas is basically a very detailled parcellation with an anatomical label for each
% vertex.
%
% An example of a surface based Brodmann parcellation looks like this
%
%              pos: [8192x3]         positions of the vertices forming the cortical sheet
%              tri: [16382x3]        triangles of the cortical sheet
%         coordsys: 'ctf'            the (head) coordinate system in which the vertex positions are expressed
%             unit: 'mm'             the units in which the coordinate system is expressed
%         brodmann: [8192x1 uint8]   values from 1 to N, the value 0 means unknown
%    brodmannlabel: {Nx1 cell}
%
% An alternative representation of this parcellation is
%
%              pos: [8192x3]           positions of the vertices forming the cortical sheet
%              tri: [16382x3]          triangles of the cortical sheet
%         coordsys: 'ctf'              the (head) coordinate system in which the vertex positions are expressed
%             unit: 'mm'               the units in which the coordinate system is expressed
%  Brodmann_Area_1: [8192x1 logical]   binary map representing the voxels belonging to the specific area
%  Brodmann_Area_2: [8192x1 logical]   binary map representing the voxels belonging to the specific area
%  Brodmann_Area_3: [8192x1 logical]   binary map representing the voxels belonging to the specific area
%  ...
%
% The examples above demonstrate that a parcellation can be indexed, i.e. consisting of
% subsequent integer numbers (1, 2, ...) or probabilistic, consisting of real numbers
% ranging from 0 to 1 that represent probabilities between 0% and 100%. An extreme case
% is one where the probability is either 0 or 1, in which case the probability can be
% represented as a binary or logical array.
%
% The only difference to the source data structure is that the parcellation structure
% contains the additional fields xxx and xxxlabel. See FT_DATATYPE_SOURCE for further
% details.
%
% Required fields:
%   - pos
%
% Optional fields:
%   - tri, coordsys, unit
%
% Deprecated fields:
%   - none
%
% Obsoleted fields:
%   - none
%
% Revision history:
% (2012/latest) The initial version was defined in accordance with the representation of
% a voxel-based segmentation.
%
% See also ea_ft_datatype, FT_DATATYPE_SOURCE, ea_ft_datatype_segmentation

% Copyright (C) 2012, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ft_datatype_parcellation.m 10213 2015-02-11 19:38:33Z roboos $

% get the optional input arguments, which should be specified as key-value pairs
version           = ea_ft_getopt(varargin, 'version', 'latest');
parcellationstyle = ea_ft_getopt(varargin, 'parcellationstyle');  % can be indexed or probabilistic

if strcmp(version, 'latest')
    parcelversion = '2012';
    sourceversion = 'latest';
    clear version
else
    parcelversion = version;
    sourceversion = version;
    clear version
end

if isempty(parcellation)
    return;
end

switch parcelversion
    case '2012'

        if isfield(parcellation, 'pnt')
            parcellation.pos = parcellation.pnt;
            parcellation = rmfield(parcellation, 'pnt');
        end

        % convert the inside/outside fields, they should be logical rather than an index
        if isfield(parcellation, 'inside')
            parcellation = ea_fixinside(parcellation, 'logical');
        end

        dim = size(parcellation.pos,1);

        % make a list of fields that represent a parcellation
        fn = fieldnames(parcellation);
        fn = setdiff(fn, 'inside'); % exclude the inside field from any conversions
        sel = false(size(fn));
        for i=1:numel(fn)
            sel(i) = isnumeric(parcellation.(fn{i})) && numel(parcellation.(fn{i}))==dim;
        end
        % only consider numeric fields of the correct size
        fn = fn(sel);

        % determine whether the style of the input fields is probabilistic or indexed
        [indexed, probabilistic] = ea_determine_segmentationstyle(parcellation, fn, dim);

        % ignore the fields that do not contain a parcellation
        sel = indexed | probabilistic;
        fn            = fn(sel);
        indexed       = indexed(sel);
        probabilistic = probabilistic(sel);

        if ~any(probabilistic) && ~any(indexed)
            % rather than being described with a tissue label for each vertex
            % it can also be described with a tissue label for each surface or volme element
            for i = 1:length(fn)
                fname = fn{i};
                switch fname
                    case 'tri'
                        dim = size(parcellation.tri,1);
                    case 'hex'
                        dim = size(parcellation.hex,1);
                    case 'tet'
                        dim = size(parcellation.tet,1);
                end
            end
            [indexed, probabilistic] = ea_determine_segmentationstyle(parcellation, fn, dim);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ensure that the parcellation is internally consistent
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if any(probabilistic)
            parcellation = ea_fixsegmentation(parcellation, fn(probabilistic), 'probabilistic');
        end

        if any(indexed)
            parcellation = ea_fixsegmentation(parcellation, fn(indexed), 'indexed');
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % convert the parcellation to the desired style
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if strcmp(parcellationstyle, 'indexed') && any(probabilistic)
            parcellation  = convert_segmentationstyle(parcellation, fn(probabilistic), [dim 1], 'indexed');
        elseif strcmp(parcellationstyle, 'probabilistic') && any(indexed)
            parcellation  = convert_segmentationstyle(parcellation, fn(indexed), [dim 1], 'probabilistic');
        end % converting converting to desired style

    otherwise
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        error('unsupported version "%s" for parcellation datatype', parcelversion);
end

% the parcellation is a speciat type of volume structure, so ensure that it also fulfills the requirements for that
parcellation = ea_ft_datatype_source(parcellation, 'version', sourceversion);


function source = ea_ft_datatype_source(source, varargin)

% FT_DATATYPE_SOURCE describes the FieldTrip MATLAB structure for data that is
% represented at the source level. This is typically obtained with a beamformer of
% minimum-norm source reconstruction using FT_SOURCEANALYSIS.
%
% An example of a source structure obtained after performing DICS (a frequency
% domain beamformer scanning method) is shown here
%
%           pos: [6732x3 double]       positions at which the source activity could have been estimated
%        inside: [6732x1 logical]      boolean vector that indicates at which positions the source activity was estimated
%           dim: [xdim ydim zdim]      if the positions can be described as a 3D regular grid, this contains the
%                                       dimensionality of the 3D volume
%     cumtapcnt: [120x1 double]        information about the number of tapers per original trial
%          time: 0.100                 the latency at which the activity is estimated (in seconds)
%          freq: 30                    the frequency at which the activity is estimated (in Hz)
%           pow: [6732x120 double]     the estimated power at each source position
%     powdimord: 'pos_rpt'             defines how the numeric data has to be interpreted,
%                                       in this case 6732 dipole positions x 120 repetitions (i.e. trials)
%           cfg: [1x1 struct]          the configuration used by the function that generated this data structure
%
% Required fields:
%   - pos
%
% Optional fields:
%   - time, freq, pow, coh, eta, mom, ori, cumtapcnt, dim, transform, inside, cfg, dimord, other fields with a dimord
%
% Deprecated fields:
%   - method, outside
%
% Obsoleted fields:
%   - xgrid, ygrid, zgrid, transform, latency, frequency
%
% Historical fields:
%   - avg, cfg, cumtapcnt, df, dim, freq, frequency, inside, method,
%   outside, pos, time, trial, vol, see bug2513
%
% Revision history:
%
% (2014) The subfields in the avg and trial fields are now present in the
% main structure, e.g. source.avg.pow is now source.pow. Furthermore, the
% inside is always represented as logical vector.
%
% (2011) The source representation should always be irregular, i.e. not
% a 3-D volume, contain a "pos" field and not contain a "transform".
%
% (2010) The source structure should contain a general "dimord" or specific
% dimords for each of the fields. The source reconstruction in the avg and
% trial substructures has been moved to the toplevel.
%
% (2007) The xgrid/ygrid/zgrid fields have been removed, because they are
% redundant.
%
% (2003) The initial version was defined
%
% See also ea_ft_datatype, FT_DATATYPE_COMP, FT_DATATYPE_DIP, FT_DATATYPE_FREQ,
% FT_DATATYPE_MVAR, FT_DATATYPE_RAW, FT_DATATYPE_SOURCE, FT_DATATYPE_SPIKE,
% FT_DATATYPE_TIMELOCK, FT_DATATYPE_VOLUME

% Copyright (C) 2013-2014, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ft_datatype_source.m 10265 2015-03-04 12:20:41Z jansch $

% FIXME: I am not sure whether the removal of the xgrid/ygrid/zgrid fields
% was really in 2007

% get the optional input arguments, which should be specified as key-value pairs
version = ea_ft_getopt(varargin, 'version', 'latest');

if strcmp(version, 'latest') || strcmp(version, 'upcoming')
    version = '2014';
end

if isempty(source)
    return;
end

% old data structures may use latency/frequency instead of time/freq. It is
% unclear when these were introduced and removed again, but they were never
% used by any fieldtrip function itself
if isfield(source, 'frequency')
    source.freq = source.frequency;
    source      = rmfield(source, 'frequency');
end
if isfield(source, 'latency')
    source.time = source.latency;
    source      = rmfield(source, 'latency');
end

switch version
    case '2014'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ensure that it has individual source positions
        source = ea_fixpos(source);

        % ensure that it is always logical
        source = ea_fixinside(source, 'logical');

        % remove obsolete fields
        if isfield(source, 'method')
            source = rmfield(source, 'method');
        end
        if isfield(source, 'transform')
            source = rmfield(source, 'transform');
        end
        if isfield(source, 'xgrid')
            source = rmfield(source, 'xgrid');
        end
        if isfield(source, 'ygrid')
            source = rmfield(source, 'ygrid');
        end
        if isfield(source, 'zgrid')
            source = rmfield(source, 'zgrid');
        end

        if isfield(source, 'avg') && isstruct(source.avg)
            % move the average fields to the main structure
            fn = fieldnames(source.avg);
            for i=1:length(fn)
                dat = source.avg.(fn{i});
                if isequal(size(dat), [1 size(source.pos,1)])
                    source.(fn{i}) = dat';
                else
                    source.(fn{i}) = dat;
                end
                clear dat
            end % j
            source = rmfield(source, 'avg');
        end

        if isfield(source, 'inside')
            % the inside is by definition logically indexed
            probe = find(source.inside, 1, 'first');
        else
            % just take the first source position
            probe = 1;
        end

        if isfield(source, 'trial') && isstruct(source.trial)
            npos = size(source.pos,1);

            % concatenate the fields for each trial and move them to the main structure
            fn = fieldnames(source.trial);

            for i=1:length(fn)
                % some fields are descriptive and hence identical over trials
                if strcmp(fn{i}, 'csdlabel')
                    source.csdlabel = dat;
                    continue
                end

                % start with the first trial
                dat    = source.trial(1).(fn{i});
                datsiz = ea_getdimsiz(source, fn{i});
                nrpt   = datsiz(1);
                datsiz = datsiz(2:end);


                if iscell(dat)
                    datsiz(1) = nrpt; % swap the size of pos with the size of rpt
                    val  = cell(npos,1);
                    indx = find(source.inside);
                    for k=1:length(indx)
                        val{indx(k)}          = nan(datsiz);
                        val{indx(k)}(1,:,:,:) = dat{indx(k)};
                    end
                    % concatenate all data as {pos}_rpt_etc
                    for j=2:nrpt
                        dat = source.trial(j).(fn{i});
                        for k=1:length(indx)
                            val{indx(k)}(j,:,:,:) = dat{indx(k)};
                        end

                    end % for all trials
                    source.(fn{i}) = val;

                else
                    % concatenate all data as pos_rpt_etc
                    val = nan([datsiz(1) nrpt datsiz(2:end)]);
                    val(:,1,:,:,:) = dat(:,:,:,:);
                    for j=2:length(source.trial)
                        dat = source.trial(j).(fn{i});
                        val(:,j,:,:,:) = dat(:,:,:,:);
                    end % for all trials
                    source.(fn{i}) = val;

                    %         else
                    %           siz = size(dat);
                    %           if prod(siz)==npos
                    %             siz = [npos nrpt];
                    %           elseif siz(1)==npos
                    %             siz = [npos nrpt siz(2:end)];
                    %           end
                    %           val = nan(siz);
                    %           % concatenate all data as pos_rpt_etc
                    %           val(:,1,:,:,:) = dat(:);
                    %           for j=2:length(source.trial)
                    %             dat = source.trial(j).(fn{i});
                    %             val(:,j,:,:,:) = dat(:);
                    %           end % for all trials
                    %           source.(fn{i}) = val;

                end
            end % for each field

            source = rmfield(source, 'trial');

        end % if trial

        % ensure that it has a dimord (or multiple for the different fields)
        source = ea_fixdimord(source);


    case '2011'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ensure that it has individual source positions
        source = ea_fixpos(source);

        % remove obsolete fields
        if isfield(source, 'xgrid')
            source = rmfield(source, 'xgrid');
        end
        if isfield(source, 'ygrid')
            source = rmfield(source, 'ygrid');
        end
        if isfield(source, 'zgrid')
            source = rmfield(source, 'zgrid');
        end
        if isfield(source, 'transform')
            source = rmfield(source, 'transform');
        end

        % ensure that it has a dimord (or multiple for the different fields)
        source = ea_fixdimord(source);

    case '2010'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ensure that it has individual source positions
        source = ea_fixpos(source);

        % remove obsolete fields
        if isfield(source, 'xgrid')
            source = rmfield(source, 'xgrid');
        end
        if isfield(source, 'ygrid')
            source = rmfield(source, 'ygrid');
        end
        if isfield(source, 'zgrid')
            source = rmfield(source, 'zgrid');
        end

        % ensure that it has a dimord (or multiple for the different fields)
        source = ea_fixdimord(source);

    case '2007'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ensure that it has individual source positions
        source = ea_fixpos(source);

        % remove obsolete fields
        if isfield(source, 'dimord')
            source = rmfield(source, 'dimord');
        end
        if isfield(source, 'xgrid')
            source = rmfield(source, 'xgrid');
        end
        if isfield(source, 'ygrid')
            source = rmfield(source, 'ygrid');
        end
        if isfield(source, 'zgrid')
            source = rmfield(source, 'zgrid');
        end

    case '2003'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isfield(source, 'dimord')
            source = rmfield(source, 'dimord');
        end

        if ~isfield(source, 'xgrid') || ~isfield(source, 'ygrid') || ~isfield(source, 'zgrid')
            if isfield(source, 'dim')
                minx = min(source.pos(:,1));
                maxx = max(source.pos(:,1));
                miny = min(source.pos(:,2));
                maxy = max(source.pos(:,2));
                minz = min(source.pos(:,3));
                maxz = max(source.pos(:,3));
                source.xgrid = linspace(minx, maxx, source.dim(1));
                source.ygrid = linspace(miny, maxy, source.dim(2));
                source.zgrid = linspace(minz, maxz, source.dim(3));
            end
        end

    otherwise
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        error('unsupported version "%s" for source datatype', version);
end

function dimsiz = ea_getdimsiz(data, field)

% GETDIMSIZ
%
% Use as
%   dimsiz = getdimsiz(data, field)
%
% See also GETDIMORD

if ~isfield(data, field) && isfield(data, 'avg') && isfield(data.avg, field)
    field = ['avg.' field];
elseif ~isfield(data, field) && isfield(data, 'trial') && isfield(data.trial, field)
    field = ['trial.' field];
elseif ~isfield(data, field)
    error('field "%s" not present in data', field);
end

if strncmp(field, 'avg.', 4)
    prefix = [];
    field = field(5:end); % strip the avg
    data.(field) = data.avg.(field); % move the avg into the main structure
    data = rmfield(data, 'avg');
elseif strncmp(field, 'trial.', 6)
    prefix = numel(data.trial);
    field = field(7:end); % strip the trial
    data.(field) = data.trial(1).(field); % move the first trial into the main structure
    data = rmfield(data, 'trial');
else
    prefix = [];
end

dimsiz = ea_cellmatsize(data.(field));

% add nrpt in case of source.trial
dimsiz = [prefix dimsiz];


function siz = ea_cellmatsize(x)
if iscell(x)
    cellsize = numel(x);          % the number of elements in the cell-array
    [dum, indx] = max(cellfun(@numel, x));
    matsize = size(x{indx});      % the size of the content of the cell-array
    siz = [cellsize matsize];     % concatenate the two
else
    siz = size(x);
end

function [data] = ea_fixdimord(data)

% FIXDIMORD ensures consistency between the dimord string and the axes
% that describe the data dimensions. The main purpose of this function
% is to ensure backward compatibility of all functions with data that has
% been processed by older FieldTrip versions
%
% Use as
%   [data] = fixdimord(data)
% This will modify the data.dimord field to ensure consistency.
% The name of the axis is the same as the name of the dimord, i.e. if
% dimord='freq_time', then data.freq and data.time should be present.
%
% The default dimensions in the data are described by
%  'time'
%  'freq'
%  'chan'
%  'chancmb'
%  'refchan'
%  'subj'
%  'rpt'
%  'rpttap'
%  'pos'
%  'ori'
%  'rgb'
%  'comp'
%  'voxel'

% Copyright (C) 2009-2014, Robert Oostenveld, Jan-Mathijs Schoffelen
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: fixdimord.m 9972 2014-11-19 08:09:34Z roboos $

% if nargin<2, keepsourcedimord = 0; end
%
% if any(ea_ft_datatype(data, {'source', 'volume'})) && isfield(data, 'dimord') && ~keepsourcedimord
%   % the old source data representation does not have a dimord, whereas the new source data representation does have a dimord
%   warning(sprintf('removing dimord "%s" from source representation data', data.dimord));
%   data = rmfield(data, 'dimord');
%   return
% else
%   % it is ok
%   return
% end

if ~isfield(data, 'dimord')
    if ea_ft_datatype(data, 'raw')
        % it is raw data, which does not have a dimord -> this is ok
        return
    elseif ea_ft_datatype(data, 'comp')
        % it is component data, which resembles raw data -> this is ok
        return
    elseif ea_ft_datatype(data, 'volume')
        % it is volume data, which does not have a dimord -> this is ok
        return
    else
        fn = fieldnames(data);
        sel = true(size(fn));
        for i=1:length(fn)
            sel(i) = contains(fn{i}, 'dimord');
        end
        df = fn(sel);

        if isempty(df)
            if ea_ft_datatype(data, 'source') || ea_ft_datatype(data, 'parcellation')
                % it is old-style source data -> this is ok
                % ft_checkdata will convert it to new-style
                return
            else
                error('the data does not contain a dimord, but it also does not resemble raw or component data');
            end
        end

        % use this function recursively on the XXXdimord fields
        for i=1:length(df)
            data.dimord = data.(df{i});
            data = fixdimord(data);
            data.(df{i}) = data.dimord;
            data = rmfield(data, 'dimord');
        end
        % after the recursive call it should be ok
        return
    end
end

if strcmp(data.dimord, 'voxel')
    % this means that it is position
    data.dimord = 'pos';
end

dimtok = tokenize(data.dimord, '_');
if strncmp('{pos_pos}', data.dimord, 9)
    % keep these together for bivariate source structures
    dimtok = {'{pos_pos}', dimtok{3:end}};
end

for i=1:length(dimtok)
    switch dimtok{i}
        case {'tim' 'time' 'toi' 'latency'}
            dimtok{i} = 'time';

        case {'frq' 'freq' 'foi' 'frequency'}
            dimtok{i} = 'freq';

        case {'sgn' 'label' 'chan'}
            dimtok{i} = 'chan';

        case {'rpt' 'trial'}
            dimtok{i} = 'rpt';

        case {'subj' 'subject'}
            dimtok{i} = 'subj';

        case {'comp'}
            % don't change, it is ok

        case {'sgncmb' 'labelcmb' 'chancmb'}
            dimtok{i} = 'chancmb';

        case {'rpttap'}
            % this is a 2-D field, coding trials and tapers along the same dimension
            % don't change, it is ok

        case {'refchan'}
            % don't change, it is ok

        case {'ori'}
            % don't change, it is ok

        case {'rgb'}
            % don't change, it is ok

        case {'voxel' 'vox' 'repl' 'wcond'}
            % these are used in some fieldtrip functions, but are not considered standard
            warning_once('unexpected dimord "%s"', data.dimord);

        case {'pos'}
            % this is for source data on a 3-d grid, a cortical sheet, or unstructured positions

        case {'{pos}' '{pos}_rpt' '{pos}_rpttap'}
            % this is for source data on a 3-d grid, a cortical sheet, or unstructured positions
            % the data itself is represented in a cell-array, e.g. source.mom or source.leadfield

        case {'{pos_pos}'}
            % this is for bivariate source data on a 3-d grid, a cortical sheet, or unstructured positions

        otherwise
            error(sprintf('unexpected dimord "%s"', data.dimord));

    end % switch dimtok
end % for length dimtok

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(data, 'tim'),         data.time      = data.tim         ; data = rmfield(data, 'tim')        ; end
if isfield(data, 'toi'),         data.time      = data.toi         ; data = rmfield(data, 'toi')        ; end
if isfield(data, 'latency'),     data.time      = data.latency     ; data = rmfield(data, 'latency')    ; end
if isfield(data, 'frq'),         data.freq      = data.frq         ; data = rmfield(data, 'frq')        ; end
if isfield(data, 'foi'),         data.freq      = data.foi         ; data = rmfield(data, 'foi')        ; end
if isfield(data, 'frequency'),   data.freq      = data.frequency   ; data = rmfield(data, 'frequency')  ; end
if isfield(data, 'sgn'),         data.label     = data.sgn         ; data = rmfield(data, 'sgn')        ; end
if isfield(data, 'chan'),        data.label     = data.chan        ; data = rmfield(data, 'chan')       ; end
% if isfield(data, 'trial'),         data.rpt     = data.trial         ; data = rmfield(data, 'trial')        ; end  % DO NOT CONVERT -> this is an exception
if isfield(data, 'subject'),     data.subj      = data.subject     ; data = rmfield(data, 'subject')    ; end
if isfield(data, 'sgncmb'),      data.labelcmb  = data.sgncmb      ; data = rmfield(data, 'sgncmb')     ; end
if isfield(data, 'chancmb'),     data.labelcmb  = data.chancmb     ; data = rmfield(data, 'chancmb')    ; end

% ensure that it is a column
if isfield(data, 'label')
    data.label = data.label(:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if isfield(data, 'trial')
%   mat = data.trial;
% elseif isfield(data, 'individual')
%   mat = data.individual;
% elseif isfield(data, 'avg')
%   mat = data.avg;
% elseif isfield(data, 'crsspctrm')
%   mat = data.crsspctrm;
% elseif isfield(data, 'powspctrm')
%   mat = data.powspctrm;
% elseif isfield(data, 'fourierspctrm')
%   mat = data.fourierspctrm;
% end
%
% add the descriptive axis for each dimension
% for i=1:length(dimtok)
%   if isfield(data, dimtok{i})
%     % the dimension is already described with its own axis
%     % data = setfield(data, dimtok{i}, getfield(data, dimtok{i}));
%   else
%     % add an axis to the output data
%     data = setfield(data, dimtok{i}, 1:size(mat,i));
%   end
% end

% undo the tokenization
data.dimord = dimtok{1};
for i=2:length(dimtok)
    data.dimord = [data.dimord '_' dimtok{i}];
end

function [type, dimord] = ea_ft_datatype(data, desired)

% ea_ft_datatype determines the type of data represented in a FieldTrip data
% structure and returns a string with raw, freq, timelock source, comp,
% spike, source, volume, dip.
%
% Use as
%   [type, dimord] = ea_ft_datatype(data)
%   [status]       = ea_ft_datatype(data, desired)
%
% See also FT_DATATYPE_COMP FT_DATATYPE_FREQ FT_DATATYPE_MVAR
% ea_ft_datatype_segmentation FT_DATATYPE_PARCELLATION FT_DATATYPE_SOURCE
% FT_DATATYPE_TIMELOCK FT_DATATYPE_DIP FT_DATATYPE_HEADMODEL
% FT_DATATYPE_RAW FT_DATATYPE_SENS FT_DATATYPE_SPIKE FT_DATATYPE_VOLUME

% Copyright (C) 2008-2012, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ea_ft_datatype.m 10064 2014-12-22 14:30:50Z roboos $

if nargin<2
    desired = [];
end

% determine the type of input data, this can be raw, freq, timelock, comp, spike, source, volume, dip, segmentation, parcellation
israw          =  isfield(data, 'label') && isfield(data, 'time') && isa(data.time, 'cell') && isfield(data, 'trial') && isa(data.trial, 'cell') && ~isfield(data,'trialtime');
isfreq         = (isfield(data, 'label') || isfield(data, 'labelcmb')) && isfield(data, 'freq') && ~isfield(data,'trialtime') && ~isfield(data,'origtrial'); %&& (isfield(data, 'powspctrm') || isfield(data, 'crsspctrm') || isfield(data, 'cohspctrm') || isfield(data, 'fourierspctrm') || isfield(data, 'powcovspctrm'));
istimelock     =  isfield(data, 'label') && isfield(data, 'time') && ~isfield(data, 'freq') && ~isfield(data,'timestamp') && ~isfield(data,'trialtime') && ~(isfield(data, 'trial') && iscell(data.trial)); %&& ((isfield(data, 'avg') && isnumeric(data.avg)) || (isfield(data, 'trial') && isnumeric(data.trial) || (isfield(data, 'cov') && isnumeric(data.cov))));
iscomp         =  isfield(data, 'label') && isfield(data, 'topo') || isfield(data, 'topolabel');
isvolume       =  isfield(data, 'transform') && isfield(data, 'dim') && ~isfield(data, 'pos');
issource       =  isfield(data, 'pos');
isdip          =  isfield(data, 'dip');
ismvar         =  isfield(data, 'dimord') && contains(data.dimord, 'lag');
isfreqmvar     =  isfield(data, 'freq') && isfield(data, 'transfer');
ischan         = ea_check_chan(data);
issegmentation = ea_check_segmentation(data);
isparcellation = ea_check_parcellation(data);

if ~isfreq
    % this applies to a freq structure from 2003 up to early 2006
    isfreq = all(isfield(data, {'foi', 'label', 'dimord'})) && contains(data.dimord, 'frq');
end

% check if it is a spike structure
spk_hastimestamp  = isfield(data,'label') && isfield(data, 'timestamp') && isa(data.timestamp, 'cell');
spk_hastrials     = isfield(data,'label') && isfield(data, 'time') && isa(data.time, 'cell') && isfield(data, 'trial') && isa(data.trial, 'cell') && isfield(data, 'trialtime') && isa(data.trialtime, 'numeric');
spk_hasorig       = isfield(data,'origtrial') && isfield(data,'origtime'); % for compatibility
isspike           = isfield(data, 'label') && (spk_hastimestamp || spk_hastrials || spk_hasorig);

% check if it is a sensor array
isgrad = isfield(data, 'label') && isfield(data, 'coilpos') && isfield(data, 'coilori');
iselec = isfield(data, 'label') && isfield(data, 'elecpos');

if isspike
    type = 'spike';
elseif israw && iscomp
    type = 'raw+comp';
elseif istimelock && iscomp
    type = 'timelock+comp';
elseif isfreq && iscomp
    type = 'freq+comp';
elseif israw
    type = 'raw';
elseif iscomp
    type = 'comp';
elseif isfreqmvar
    % freqmvar should conditionally go before freq, otherwise the returned ea_ft_datatype will be freq in the case of frequency mvar data
    type = 'freqmvar';
elseif isfreq
    type = 'freq';
elseif ismvar
    type = 'mvar';
elseif isdip
    % dip should conditionally go before timelock, otherwise the ea_ft_datatype will be timelock
    type = 'dip';
elseif istimelock
    type = 'timelock';
elseif issegmentation
    % a segmentation data structure is a volume data structure, but in general not vice versa
    % segmentation should conditionally go before volume, otherwise the returned ea_ft_datatype will be volume
    type = 'segmentation';
elseif isvolume
    type = 'volume';
elseif isparcellation
    % a parcellation data structure is a source data structure, but in general not vice versa
    % parcellation should conditionally go before source, otherwise the returned ea_ft_datatype will be source
    type = 'parcellation';
elseif issource
    type = 'source';
elseif ischan
    % this results from avgovertime/avgoverfreq after timelockstatistics or freqstatistics
    type = 'chan';
elseif iselec
    type = 'elec';
elseif isgrad
    type = 'grad';
else
    type = 'unknown';
end

if nargin>1
    % return a boolean value
    switch desired
        case 'raw'
            type = any(strcmp(type, {'raw', 'raw+comp'}));
        case 'timelock'
            type = any(strcmp(type, {'timelock', 'timelock+comp'}));
        case 'freq'
            type = any(strcmp(type, {'freq', 'freq+comp'}));
        case 'comp'
            type = any(strcmp(type, {'comp', 'raw+comp', 'timelock+comp', 'freq+comp'}));
        case 'volume'
            type = any(strcmp(type, {'volume', 'segmentation'}));
        case 'source'
            type = any(strcmp(type, {'source', 'parcellation'}));
        case 'sens'
            type = any(strcmp(type, {'elec', 'grad'}));
        otherwise
            type = strcmp(type, desired);
    end % switch
end

if nargout>1
    % FIXME this should be replaced with getdimord in the calling code
    % also return the dimord of the input data
    if isfield(data, 'dimord')
        dimord = data.dimord;
    else
        dimord = 'unknown';
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res] = ea_check_chan(data)

if ~isstruct(data) || any(isfield(data, {'time', 'freq', 'pos', 'dim', 'transform'}))
    res = false;
elseif isfield(data, 'dimord') && any(strcmp(data.dimord, {'chan', 'chan_chan'}))
    res = true;
else
    res = false;
    fn = fieldnames(data);
    for i=1:numel(fn)
        if isfield(data, [fn{i} 'dimord']) && any(strcmp(data.([fn{i} 'dimord']), {'chan', 'chan_chan'}))
            res = true;
            break;
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res] = ea_check_segmentation(volume)
res = false;

if ~isfield(volume, 'dim')
    return
end

if isfield(volume, 'pos')
    return
end

if any(isfield(volume, {'seg', 'csf', 'white', 'gray', 'skull', 'scalp', 'brain'}))
    res = true;
    return
end

fn = fieldnames(volume);
isboolean = [];
cnt = 0;
for i=1:length(fn)
    if isfield(volume, [fn{i} 'label'])
        res = true;
        return
    else
        if (islogical(volume.(fn{i})) || isnumeric(volume.(fn{i}))) && isequal(size(volume.(fn{i})),volume.dim)
            cnt = cnt+1;
            if islogical(volume.(fn{i}))
                isboolean(cnt) = true;
            else
                isboolean(cnt) = false;
            end
        end
    end
end
if ~isempty(isboolean)
    res = all(isboolean);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res] = ea_check_parcellation(source)
res = false;

if ~isfield(source, 'pos')
    return
end

fn = fieldnames(source);
fb = false(size(fn));
npos = size(source.pos,1);
for i=1:numel(fn)
    % for each of the fields check whether it might be a logical array with the size of the number of sources
    tmp = source.(fn{i});
    fb(i) = numel(tmp)==npos && islogical(tmp);
end
if sum(fb)>1
    % the presence of multiple logical arrays suggests it is a parcellation
    res = true;
end

if res == false      % check if source has more D elements
    check = 0;
    for i = 1: length(fn)
        fname = fn{i};
        switch fname
            case 'tri'
                npos = size(source.tri,1);
                check = 1;
            case 'hex'
                npos = size(source.hex,1);
                check = 1;
            case 'tet'
                npos = size(source.tet,1);
                check = 1;
        end
    end
    if check == 1   % check if elements are labelled
        for i=1:numel(fn)
            tmp = source.(fn{i});
            fb(i) = numel(tmp)==npos && islogical(tmp);
        end
        if sum(fb)>1
            res = true;
        end
    end
end

fn = fieldnames(source);
for i=1:length(fn)
    if isfield(source, [fn{i} 'label']) && isnumeric(source.(fn{i}))
        res = true;
        return
    end
end


function source = ea_fixpos(source)
if ~isfield(source, 'pos')
    if isfield(source, 'xgrid') && isfield(source, 'ygrid') && isfield(source, 'zgrid')
        source.pos = ea_grid2pos(source.xgrid, source.ygrid, source.zgrid);
    elseif isfield(source, 'dim') && isfield(source, 'transform')
        source.pos = ea_dim2pos(source.dim, source.transform);
    else
        error('cannot reconstruct individual source positions');
    end
end

function pos = ea_grid2pos(xgrid, ygrid, zgrid)
[X, Y, Z] = ndgrid(xgrid, ygrid, zgrid);
pos = [X(:) Y(:) Z(:)];


function pos = ea_dim2pos(dim, transform)
[X, Y, Z] = ndgrid(1:dim(1), 1:dim(2), 1:dim(3));
pos = [X(:) Y(:) Z(:)];
pos = ea_ft_warp_apply(transform, pos, 'homogenous');


function [indexed, probabilistic] = ea_determine_segmentationstyle(segmentation, fn, dim)

% DETERMINE_SEGMENTATIONSTYLE is a helper function that determines the type of segmentation
% contained in each of the fields. It is used by ea_ft_datatype_segmentation and
% ft_datatype_parcellation.
%
% See also FIXSEGMENTATION, CONVERT_SEGMENTATIONSTYLE

indexed       = false(size(fn));
probabilistic = false(size(fn));

% determine for each of the fields whether it is probabilistic, indexed or something else
for i=1:numel(fn)
    if numel(segmentation.(fn{i}))~=prod(dim)
        % this does not look like a segmentation
        continue
    elseif strcmp(fn{i}, 'anatomy')
        % this should not be interpreted as segmentation, also not when it is a uint8 or uint16 representation
        continue
    else
        if isfield(segmentation, [fn{i} 'label'])
            % the xxxlabel field exists, which only makes sense for an indexed representation
            probabilistic(i) = false;
            indexed(i)       = true;
        else
            % this looks like a segmentation
            tmp = segmentation.(fn{i});
            tmp = tmp(:);       % convert to vector
            sel = isnan(tmp);   % find NaN values
            if any(sel)
                % note that the the finding and removing of NaNs have been separated to speed up the code
                tmp = tmp(~sel);  % remove NaN values
            end
            clear sel
            probabilistic(i) =  islogical(tmp) || all(tmp>=-0.001 & tmp<=1.001); % allow some roundoff error
            indexed(i)       = ~islogical(tmp) && all(abs(tmp - round(tmp))<1000*eps);

            if probabilistic(i) && indexed(i)
                % the xxxlabel does not exist, so treat it as a probabilistic representation
                probabilistic(i) = true;
                indexed(i)       = false;
            end
        end
    end
end % for each of the fields

function rows = ea_sb_sparse_to_mat(diinsy)

% SB_SPARSE_TO_MAT
%
% $Id: sb_sparse_to_mat.m 8776 2013-11-14 09:04:48Z roboos $

rows = zeros(max(diinsy),1);
rows(diinsy) = 1;
rows = [1;rows];
rows = cumsum(rows);
rows = rows(1:end-1);


function [stiff, diinsy, cols, sysmat] = ea_sb_calc_stiff(vol)

% SB_CALC_STIFF
%
% $Id: sb_calc_stiff.m 8776 2013-11-14 09:04:48Z roboos $

if(~(size(vol.pos,2)==3))
    if(size(vol.pos,1)==3)
        node = vol.pos';
        warning('Dimension of vol.pos should be #nodes x 3!')
    else
        error('vol.pos has wrong dimension!')
    end
else
    node = vol.pos;
end
npnt = size(node,1);
npnt = int32(npnt);

if isfield(vol,'tet')
    if size(vol.tet,1) == 4
        mele = size(vol.tet,1);
        elem = vol.tet;
    elseif size(vol.tet,2) == 4
        mele = size(vol.tet,2);
        elem = vol.tet';
    else
        error('vol.tet has wrong dimensions!')
    end
    elem = [elem; zeros(4,size(elem,2))];
elseif isfield(vol,'hex')
    if size(vol.hex,1) == 8
        mele = size(vol.hex,1);
        elem = vol.hex;
    elseif size(vol.hex,2) == 8
        mele = size(vol.hex,2);
        elem = vol.hex';
    else
        error('vol.hex has wrong dimensions!')
    end
else
    error('Could not find connectivity information!')
end

if min(min(elem(1:mele,:))) == 0
    elem = elem + 1;
    warning('Numbering of nodes in vol.tet/vol.hex must start at 1 (Fortran numbering)!')
elseif min(min(elem(1:mele,:))) < 0
    error('No negative indices for conectivity information allowed!')
end

if isfield(vol,'cond') && isfield(vol,'tissue') && isfield(vol,'tissuelabel')
    if length(vol.tissuelabel) == length(vol.cond)
        if length(vol.tissue) == size(elem,2)
            cond = zeros(size(elem,2),1);
            numlabels = length(vol.tissuelabel);
            for i=1:numlabels
                cond(vol.tissue == i) = vol.cond(i);
            end
        else
            error('Dimensions of vol.tet or vol.hex and vol.tissue do not fit!');
        end
    else
        error('Dimensions of vol.cond and entries of vol.tissuelabel do not fit!');
    end
end

mele = int32(mele);
elem = int32(elem);

% check whether the nodes have right orientation

if isfield(vol,'tet')
    if ~ea_sb_test_ori(node,elem(1:4,:)')
        error('Elements have wrong orientation, consider exchanging node 3 and 4');
        return;
    end
elseif isfield(vol,'hex')
    if ~ea_sb_test_ori(node,elem')
        error('Elements have wrong orientation or are degenerated');
        return
    end
end

try
    [diinsy,cols,sysmat] = ea_calc_stiff_matrix_val(node,elem,cond,mele);
catch err
    if ispc && strcmp(err.identifier,'MATLAB:invalidMEXFile')
        error('Error executing mex-file. Microsoft Visual C++ 2008 Redistributables and Intel Visual Fortran Redistributables are required.')
    else
        rethrow(err)
    end
end
npnt = double(npnt);
diinsy = double(diinsy);
cols = double(cols);
rows = ea_sb_sparse_to_mat(diinsy);
stiff = sparse(rows,cols,sysmat,npnt,npnt,length(sysmat));


function err = ea_sb_test_ori(pnt,elem)
err = 1;
if(size(elem,2) == 4)
    det = sum(cross(pnt(elem(:,2),:)-pnt(elem(:,1),:),pnt(elem(:,4),:)-pnt(elem(:,1),:),2).*(pnt(elem(:,3),:)-pnt(elem(:,1),:)),2);
    if length(find(det <= 0)) > 0
        err = 0;
    end
elseif(size(elem,2) == 8)
    det1 = sum(cross(pnt(elem(:,6),:)-pnt(elem(:,1),:),pnt(elem(:,8),:)-pnt(elem(:,1),:),2).*(pnt(elem(:,5),:)-pnt(elem(:,1),:)),2);
    det2 = sum(cross(pnt(elem(:,3),:)-pnt(elem(:,1),:),pnt(elem(:,6),:)-pnt(elem(:,1),:),2).*(pnt(elem(:,2),:)-pnt(elem(:,1),:)),2);
    det3 = sum(cross(pnt(elem(:,8),:)-pnt(elem(:,1),:),pnt(elem(:,3),:)-pnt(elem(:,1),:),2).*(pnt(elem(:,4),:)-pnt(elem(:,1),:)),2);
    det4 = sum(cross(pnt(elem(:,8),:)-pnt(elem(:,3),:),pnt(elem(:,6),:)-pnt(elem(:,3),:),2).*(pnt(elem(:,7),:)-pnt(elem(:,3),:)),2);
    if (length(find(det1 <= 0)) > 0 || length(find(det2 <= 0)) > 0 || length(find(det3 <= 0)) > 0 || length(find(det4 <= 0)) > 0)
        err = 0;
    end
else
    error('Invalid number of nodes per element!');
end
