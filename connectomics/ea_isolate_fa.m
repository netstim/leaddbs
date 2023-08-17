function ea_isolate_fa(options)

disp('Export FA...');
directory=[options.root,options.patientname,filesep];

Vdti=spm_vol([directory,options.prefs.dti]);
Xdti=spm_read_vols(Vdti);

bval=load([directory,options.prefs.bval]);
bvec=load([directory,options.prefs.bvec]);

if size(bvec,1)~=3
    bvec=bvec';
end
if size(bval,1)~=1
    bval=bval';
end

for i=1:size(Xdti,4)
   DTIdata(i).VoxelData=single(squeeze(Xdti(:,:,:,i)));
   DTIdata(i).Gradient=bvec(:,i);
   DTIdata(i).Bvalue=bval(:,i);
   if ~DTIdata(i).Bvalue
       DTIdata(i).Gradient=zeros(3,1);
   end
end

% Constants DTI
parametersDTI=[];
parametersDTI.BackgroundThreshold=50;
parametersDTI.WhiteMatterExtractionThreshold=0.10;
parametersDTI.textdisplay=true;

[ADC,FA,VectorF,DifT,Bo]=ea_DTI(DTIdata,parametersDTI);

V=spm_vol([directory,options.prefs.dti,',1']);
V.fname=[directory,options.prefs.fa];
V.dt(1) = 16;
spm_write_vol(V,FA);

V.fname=[directory,options.prefs.b0];
V.dt(1) = 16;
spm_write_vol(V,Bo);
