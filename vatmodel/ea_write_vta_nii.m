function varargout=ea_write_vta_nii(S,stimname,midpts,indices,elspec,dpvx,voltix,constvol,thresh,mesh,gradient,side,resultfig,options)
stack = dbstack;
if any(ismember({'ea_generate_optim_vat', 'ea_generate_base_vats','ea_generate_vat'}, {stack.name}))
    isClearTuneRun = 1;
else
   isClearTuneRun = 0;
end
if ~isempty(resultfig)
    vatgrad=getappdata(resultfig,'vatgrad');
    if isempty(vatgrad)
        clear('vatgrad');
    end
end

% define midpoints of quiver field
vatgrad(side).x=midpts(indices,1);
vatgrad(side).y=midpts(indices,2);
vatgrad(side).z=midpts(indices,3);

vizz=0;
gradvis=gradient(indices,:);
mag_gradvis=sqrt(sum(gradvis'.^2,1))';
nmag_gradvis=mag_gradvis; % copy to get normalized version
nmag_gradvis(nmag_gradvis>thresh)=thresh;
nmag_gradvis=(nmag_gradvis-min(nmag_gradvis(:)))/max(nmag_gradvis-min(nmag_gradvis(:))); % norm from 1 - 0
nmag_gradvis=nmag_gradvis/5; % largest grad vector will be 1/50 mm long

% now apply scaling to gradvis:
gradvis=gradvis.*repmat(nmag_gradvis,1,3);
gradvis=gradvis./repmat(mag_gradvis,1,3);

vatgrad(side).qx=gradvis(:,1);
vatgrad(side).qy=gradvis(:,2);
vatgrad(side).qz=gradvis(:,3);

%figure, quiver3(midpts(:,1),midpts(:,2),midpts(:,3),gradient(:,1),gradient(:,2),gradient(:,3))

% calculate electric field ET by calculating midpoints of each
% mesh-connection and setting difference of voltage to these points.

vat.pos=midpts;

%plot3(midpts(:,1),midpts(:,2),midpts(:,3),'g.');

if ~isempty(resultfig)
    setappdata(resultfig,'vatgrad',vatgrad);
end

ngrad=sqrt(sum(gradient'.^2,1));
vat.ET=ngrad; % vol.cond(vol.tissue).*ngrad; would be stromstaerke.

% reload elstruct to make sure to take correct one (native vs. template)
[coords_mm,trajectory,markers]=ea_load_reconstruction(options);
elstruct(1).coords_mm=coords_mm;
elstruct(1).trajectory=trajectory;
elstruct(1).name=options.patientname;
elstruct(1).markers=markers;
if options.prefs.machine.vatsettings.horn_removeElectrode
    vat = jr_remove_electrode(vat,elstruct,mesh,side,elspec);
end

ea_dispt('Preparing VAT...');

vat.tET=vat.ET>thresh;
vat.tpos=vat.pos(vat.tET,:);
outliers=ea_removeoutliers(vat.tpos,mean(dpvx,1),voltix,constvol);
vat.tpos(outliers,:)=[];
if vizz
    figure, plot3(vat.tpos(:,1),vat.tpos(:,2),vat.tpos(:,3),'r.');
    [coords_mm,trajectory,markers,elmodel,manually_corrected,coords_acpc]=ea_load_reconstruction(options);
    hold on
    plot3(trajectory{side}(:,1),trajectory{side}(:,2),trajectory{side}(:,3),'k*');
    plot3(nvat.tpos(:,1),nvat.tpos(:,2),nvat.tpos(:,3),'m.');
    toptions=options;
    toptions.native=1;
    [coords_mm,trajectory,markers,elmodel,manually_corrected,coords_acpc]=ea_load_reconstruction(toptions);
    plot3(trajectory{side}(:,1),trajectory{side}(:,2),trajectory{side}(:,3),'g*');
end

% the following will be used for volume 2 isosurf creation as well as
% volumetrics of the vat in mm^3.
ea_dispt('Calculating interpolant on scattered FEM mesh data...');
F=scatteredInterpolant(vat.pos(:,1),vat.pos(:,2),vat.pos(:,3),vat.ET','linear','none');

ea_dispt('Converting to equispaced image data...');
res=100;
gv=cell(3,1); spacing=zeros(3,1);
try
    for dim=1:3
        gv{dim}=linspace(min(round(vat.tpos(:,dim)))-5,max(round(vat.tpos(:,dim)))+5,res);
        spacing(dim)=abs(gv{dim}(1)-gv{dim}(2));
    end
catch
    vatfv=nan;
    vatvolume=nan;
    radius=nan;
    return
end

ea_dispt('Creating nifti header for export...');
% create nifti
[~, ~, endian] = computer;
switch endian
    case 'L'
        endian = 0;
    case 'B'
        endian = 1;
end

chun1=randperm(res); chun2=randperm(res); chun3=randperm(res);
Vvat.mat=mldivide([(chun1);(chun2);(chun3);ones(1,res)]',[gv{1}(chun1);gv{2}(chun2);gv{3}(chun3);ones(1,res)]')';
Vvat.dim=[res,res,res];
Vvat.dt = [4, endian];
Vvat.n=[1 1];
Vvat.descrip='lead dbs - vat';

ea_dispt('Filling data with values from interpolant...');
eeg = F(gv);
eeg(isnan(eeg))=0;
eeg(eeg>options.prefs.vat.efieldmax)=options.prefs.vat.efieldmax; % upperlimit files to 10000.

% figure, plot3(F.Points(:,1),F.Points(:,2),F.Points(:,3),'r.')
% hold on
% plot3(vat.pos(:,1),vat.pos(:,2),vat.pos(:,3),'b.')

ea_dispt('Calculating output file data...');

% binary e-field - "vat"
eg=eeg;
eg=eg>thresh;

% normalized e-field (zscored).
neeg=eeg;
neeg(~eg)=nan;
neeg(neeg>0)=ea_normal(neeg(neeg>0),1,0,' ',0,1,'TRUE');%
neeg(~isnan(neeg))=neeg(~isnan(neeg))-min(neeg(~isnan(neeg)));
neeg(~isnan(neeg))=neeg(~isnan(neeg))/sum(neeg(~isnan(neeg))); % 0-1 distributed.

[xg,yg,zg] = meshgrid(gv{1},gv{2},gv{3});

XYZmax=[max(yg(eg>0)),max(xg(eg>0)),max(zg(eg>0))]; % x and y need to be permuted here (should be correct but wouldnt matter anyways since only serves to calc radius)
try
    radius=ea_pdist([XYZmax;dpvx]);
catch
    keyboard
end

ea_dispt('Calculating volume...');

vatvolume=sum(eg(:))*spacing(1)*spacing(2)*spacing(3); % returns volume of vat in mm^3
S.volume(side)=vatvolume;

ea_dispt('Writing files...');

% determine stimulation name
stimDir = fullfile(options.subj.stimDir, ea_nt(options), stimname);
ea_mkdir(stimDir);
filePrefix = ['sub-', options.subj.subjId, '_sim-'];

modelLabel = ea_simModel2Label(S.model);

switch side
    case 1
        Vvat.fname = [stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-R.nii'];
        Vvate=Vvat;
        Vvatne=Vvat;
        Vvate.fname = [stimDir, filesep, filePrefix, 'efield_model-', modelLabel, '_hemi-R.nii'];
        Vvatne.fname = [stimDir, filesep, filePrefix, 'efieldgauss_model-', modelLabel, '_hemi-R.nii'];
    case 2
        Vvat.fname = [stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-L.nii'];
        Vvate = Vvat;
        Vvatne = Vvat;
        Vvate.fname = [stimDir, filesep, filePrefix, 'efield_model-', modelLabel, '_hemi-L.nii'];
        Vvatne.fname = [stimDir, filesep, filePrefix, 'efieldgauss_model-', modelLabel, '_hemi-L.nii'];
end

ea_savestimulation(S,options);

Vvate.img=eeg; %permute(eeg,[2,1,3]);
Vvate.dt = [16, endian];
if ~isClearTuneRun || options.writeVTA
    ea_write_nii(Vvate);
end
%ea_write_nii(Vvate);
Vvatne.img=neeg; %permute(neeg,[2,1,3]);
if ~isClearTuneRun || options.writeVTA
    ea_write_nii(Vvatne);
end
%ea_write_nii(Vvatne);

Vvat.img=eg; %permute(eg,[1,2,3]);
if ~isClearTuneRun || options.writeVTA
    ea_write_nii(Vvat);
end
%ea_write_nii(Vvat);

ea_dispt('Calculating isosurface to display...');
vatfv=isosurface(xg,yg,zg,permute(Vvat.img,[2,1,3]),0.75);

caps=isocaps(xg,yg,zg,permute(Vvat.img,[2,1,3]),0.5);

vatfv.faces=[vatfv.faces;caps.faces+size(vatfv.vertices,1)];
vatfv.vertices=[vatfv.vertices;caps.vertices];

try
    vatfv=ea_smoothpatch(vatfv,1,35);
catch
    try
        cd([ea_getearoot,'ext_libs',filesep,'smoothpatch']);
        mex ea_smoothpatch_curvature_double.c -v
        mex ea_smoothpatch_inversedistance_double.c -v
        mex ea_vertex_neighbours_double.c -v
        vatfv=ea_smoothpatch(vatfv);
    catch
        warndlg('Patch could not be smoothed. Please supply a compatible Matlab compiler to smooth VTAs.');
    end
end
% new save by Till to save VAT and quiver in seperate .mat-file for quick
% visualization
switch side
    case 1
        vatfvname = [stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-R.mat'];
    case 2
        vatfvname = [stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-R.mat'];
end

vatgrad = vatgrad(side);
save(vatfvname,'vatfv','vatgrad','vatvolume');

%% new vta.nii save, filled and eroded/dilated by 3 voxels.
Vvat.img=imfill(Vvat.img,'holes');
SE = strel('sphere',3);
Vvat.img = imerode(Vvat.img,SE);
Vvat.img = imdilate(Vvat.img,SE);
if ~isClearTuneRun || options.writeVTA
    ea_write_nii(Vvat);
end

%ea_write_nii(Vvat);
if isClearTuneRun
    varargout{1} = Vvate;
else
    varargout{1} = vatfv;
    varargout{2} = vatvolume;
    varargout{3} = radius;
end 

function vat = jr_remove_electrode(vat,elstruct,mesh,side,elspec)

% anonymous functions for rotation matrices
rotx = @(t) [1 0 0; 0 cosd(t) -sind(t) ; 0 sind(t) cosd(t)] ;
roty = @(t) [cosd(t) 0 sind(t) ; 0 1 0 ; -sind(t) 0  cosd(t)] ;
rotz = @(t) [cosd(t) -sind(t) 0 ; sind(t) cosd(t) 0; 0 0 1] ;

% Remove electrode
vat.pos(mesh.tissue>2,:) = [];
vat.ET(mesh.tissue>2) = [];

oldlocas = vat.pos;

% Assign tip of electrode as origin
org = elstruct.trajectory{1,side}(1,:);
tra = elstruct.trajectory{1,side} - org;
pos = vat.pos - org;

% Rotate coordinate system so that y-axis aligns with electrode
elvec = (tra(end,:))';
yvec = [0 1 0]';

th1 = atan2d(elvec(2),elvec(1));
M1z = rotz(-th1);
th2 = atan2d(yvec(2),yvec(1));
M2z = rotz(-th2);
v1 = M1z*elvec;
v2 = M2z*yvec;
b = atan2d(v2(1),v2(3));
a = atan2d(v1(1),v1(3));
My = roty(b-a);
R = M2z'*My*M1z;
r_elvec = R*elvec;

if r_elvec(1) + r_elvec(3) > 0.01
    disp('Error in electrode removal: Rotation did not align along y-axis')
end

r_pos = R*pos';

cr_pos = r_pos(:,r_pos(2,:)>0); % Determine all points at electrode level
fact = (vecnorm(cr_pos([1 3],:),2)-elspec.lead_diameter/2)./vecnorm(cr_pos([1 3],:),2); % Determine shifting factors
cr_pos([1 3],:) = cr_pos([1 3],:).*fact;
cr_pos(:,fact<0) = nan;

r_pos(:,r_pos(2,:)>0) = cr_pos;
artpts = find(isnan(r_pos(1,:)));

% Rotate back
yvec = (tra(end,:))';
elvec = [0 1 0]';
th1 = atan2d(elvec(2),elvec(1));
M1z = rotz(-th1);
th2 = atan2d(yvec(2),yvec(1));
M2z = rotz(-th2);
v1 = M1z*elvec;
v2 = M2z*yvec;
b = atan2d(v2(1),v2(3));
a = atan2d(v1(1),v1(3));
My = roty(b-a);
R = M2z'*My*M1z;

vat.pos = R*r_pos+org';

% Remove potential artefacts
vat.pos(:,artpts) = [];
vat.ET(artpts) = [];

oldlocas(artpts,:) = [];

moved = vecnorm(vat.pos-oldlocas',2);
moved(moved<0.2) = [];
moved = abs(moved-elspec.lead_diameter/2);

if max(moved)>0.2
    disp('something went wrong during electrode removal... point movement too great or too small');
end

vat.pos = vat.pos';


function outliers=ea_removeoutliers(pointcloud,dpvx,voltix,constvol)
% using the heuristic Maedler/Coenen model to detect outliers
vmax=max(abs(voltix));

if ~constvol
   vmax=vmax*1000;
end
r=maedler12_eq3(vmax,1000);
r=r*4.5;
mp=dpvx;

D=pointcloud-repmat(mp,size(pointcloud,1),1);
S=[r,r,r];
outliers=D>repmat(S,size(D,1),1);
outliers=any(outliers,2);


function r=maedler12_eq3(U,Im)
% This function calculates the  radius of Volume of Activated Tissue for
% stimulation settings U (Maedler 2012). Clinical measurements of DBS
% electrode impedance typically range from 500-1500 Ohm (Butson 2006).

r=0; %
if U %(U>0)
    k1=-1.0473;
    k3=0.2786;
    k4=0.0009856;
    r=-(k4*Im-sqrt(k4^2*Im^2 + 2*k1*k4*Im + k1^2 + 4*k3*U) + k1)/(2*k3);
end
