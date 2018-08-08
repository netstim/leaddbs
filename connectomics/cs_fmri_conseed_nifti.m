function cs_fmri_conseed_nifti(restfile,seedfile,options)

directory=[fileparts(restfile),filesep];
options=ea_getptopts(directory,options);
options=ea_assignpretra(options);

[~,restfname]=fileparts(restfile);
[~,anatfname]=fileparts(options.prefs.prenii_unnormalized);

if ~exist([directory,'sr',restfname,'.nii'],'file') ...
    || ~exist([directory,'r',restfname,'_c1',anatfname,'.nii'],'file') % preproecessing needs to be performed
    disp('No preprocessed fMRI-images found, processing...');
    options.prefs.rest=[restfname,'.nii'];
    options.overwriteapproved = 0;
    ea_preprocess_fmri(options);
    disp('Done preprocessing fMRI data.');
end

%% set some initial parameters here:
usesmooth=1;
if usesmooth
    spfx='sr';
else
    spfx='r';
end
rest=ea_load_nii([directory,spfx,restfname,'.nii']);
signallength=size(rest.img,4);
interpol_tc=nan(numel(rest.img(:,:,:,1)),size(rest.img,4));
for tmpt = 1:signallength
    thisvol=rest.img(:,:,:,tmpt);
    interpol_tc(:,tmpt)=thisvol(:);
end
%% Extract timecourses of complete volume for signal regression..
%alltc=spm_read_vols(spm_vol(restfilename));

load([directory,'TR.mat']);
%% Data corrections steps


disp('Calculating C2 and CSF-signals for signal regression...');

% regression steps

c2=ea_load_nii([directory,'r',restfname,'_c2',options.prefs.prenii_unnormalized]);
c3=ea_load_nii([directory,'r',restfname,'_c3',options.prefs.prenii_unnormalized]);

ec2map=c2.img(:); ec2map(ec2map<0.6)=0; ec2map=logical(ec2map);
ec3map=c3.img(:); ec3map(ec3map<0.6)=0; ec3map=logical(ec3map);


%% regress out WM- and CSF-Timecourses
WMTimecourse=zeros(signallength,1);
CSFTimecourse=zeros(signallength,1);
for tmpt = 1:signallength
    OneTimePoint=interpol_tc(:,tmpt);
    WMTimecourse(tmpt)=squeeze(nanmean(nanmean(nanmean(OneTimePoint(ec2map)))));
    CSFTimecourse(tmpt)=squeeze(nanmean(nanmean(nanmean(OneTimePoint(ec3map)))));
end

disp('Done. Regressing out nuisance variables...');

% %% actual regression:
% for voxx=1:size(interpol_tc,1)
%
%     beta_hat        = (X'*X)\X'*squeeze(interpol_tc(voxx,:))';
%     if ~isnan(beta_hat)
%     interpol_tc(voxx,:)=squeeze(interpol_tc(voxx,:))'-X*beta_hat;
%     else
%         warning('Regression of WM-/CSF-Signals could not be performed.');
%     end
% end
%
% clear X

%% regress out movement parameters

load([directory,'rp_',restfname,'.txt']); % rigid body motion parameters.
rp_rest=eval(['rp_',restfname]);
X(:,1)=ones(signallength,1);
X(:,2)=WMTimecourse;
X(:,3)=CSFTimecourse;
X(:,4)=rp_rest(1:signallength,1);
X(:,5)=rp_rest(1:signallength,2);
X(:,6)=rp_rest(1:signallength,3);
X(:,7)=rp_rest(1:signallength,4);
X(:,8)=rp_rest(1:signallength,5);
X(:,9)=rp_rest(1:signallength,6);

for voxx=1:size(interpol_tc,1)

    beta_hat        = (X'*X)\X'*squeeze(interpol_tc(voxx,:))';
    if ~isnan(beta_hat)
    interpol_tc(voxx,:)=squeeze(interpol_tc(voxx,:))'-X*beta_hat;
    else
        warning('Regression of Motion parameters could not be performed.');
    end

end


%% begin rest bandpass

lp_HighCutoff=0.08;
hp_LowCutoff=0.009;


disp('Done. Bandpass-filtering...');
sampleFreq 	 = 1/TR;
sampleLength = signallength;
paddedLength = rest_nextpow2_one35(sampleLength); %2^nextpow2(sampleLength);

voxelmask.locsvx=ones(size(interpol_tc,1),1);
voxelmask.locsvx(:)=1;
maskLowPass =	repmat(voxelmask.locsvx, [1, paddedLength]);
maskHighPass=	maskLowPass;
clear mask.locsvx;
%% GENERATE LOW PASS WINDOW	20070514, reference: fourior_filter.c in AFNI
%Revised by YAN Chao-Gan, 100420. Fixed the bug in calculating the frequency band.

% Low pass, such as freq < 0.08 Hz
idxCutoff	=round(lp_HighCutoff *paddedLength *TR + 1); % Revised by YAN Chao-Gan, 100420. Fixed the bug in calculating the frequency band. %idxCutoff	=round(ALowPass_HighCutoff *paddedLength *TR);
idxCutoff2	=paddedLength+2 -idxCutoff;				%Center Index =(paddedLength/2 +1)
maskLowPass(:,idxCutoff+1:idxCutoff2-1)=0; %High eliminate

%%GENERATE HIGH PASS WINDOW

% high pass, such as freq > 0.01 Hz
idxCutoff	=round(hp_LowCutoff *paddedLength *TR + 1); % Revised by YAN Chao-Gan, 100420. Fixed the bug in calculating the frequency band. %idxCutoff	=round(AHighPass_LowCutoff *paddedLength *TR);
idxCutoff2	=paddedLength+2 -idxCutoff;				%Center Index =(paddedLength/2 +1)
maskHighPass(:,1:idxCutoff-1)=0;	%Low eliminate
maskHighPass(:,idxCutoff2+1:paddedLength)=0;	%Low eliminate


% 	%20070513	remove trend --> FFT --> filter --> inverse FFT --> retrend
% YAN Chao-Gan, 100401. remove the mean --> FFT --> filter --> inverse FFT --> add mean back
fftw('dwisdom');

theMean=mean(interpol_tc,2);
interpol_tc=interpol_tc-repmat(theMean,[1, sampleLength]);
interpol_tc=cat(2,interpol_tc,zeros(size(interpol_tc,1),paddedLength-sampleLength));

%FFT
interpol_tc =fft(interpol_tc, [], 2);

%Apply the filter Low Pass
interpol_tc(~maskLowPass)=0;


%Apply the filter High Pass
interpol_tc(~maskHighPass)=0;

%inverse FFT
interpol_tc =ifft(interpol_tc, [], 2);
interpol_tc =interpol_tc(:, 1:sampleLength);%remove the padded parts

% Add the mean back after filter.
interpol_tc=interpol_tc+repmat(theMean,[1, sampleLength]);

%% end  bandpass
disp('Done.');


%% export ROI map:
for s=1:length(seedfile)
    seed{s}=ea_load_nii(seedfile{s});
    seed{s}.img(seed{s}.img==0)=nan;
    seed_tc=interpol_tc.*repmat(seed{s}.img(:),1,signallength);
    seed_tc=ea_nanmean(seed_tc,1);
    R=corr(seed_tc',interpol_tc','rows','pairwise');

    seed{s}.img(:)=R;
    [pth,sf]=fileparts(seed{s}.fname);
    outputfolder = options.lcm.func.connectome;
    if strfind(sf,'rest')
        sf(strfind(sf,'rest')-1:end)=[];
    end
    seed{s}.fname=fullfile(pth,outputfolder,[sf,'_AvgR_native_unsmoothed.nii']);

    if ~exist(fullfile(pth,outputfolder),'dir')
        mkdir(fullfile(pth,outputfolder))
    end
    ea_write_nii(seed{s});
    seed{s}.img(:)=atanh(seed{s}.img(:));
    seed{s}.fname=fullfile(pth,outputfolder,[sf,'_AvgR_Fz_native_unsmoothed.nii']);
    ea_write_nii(seed{s});

    matlabbatch{1}.spm.spatial.smooth.data = {fullfile(pth,outputfolder,[sf,'_AvgR_native_unsmoothed.nii'])
        fullfile(pth,outputfolder,[sf,'_AvgR_Fz_native_unsmoothed.nii'])};
    matlabbatch{1}.spm.spatial.smooth.fwhm = [8 8 8];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    spm_jobman('run',{matlabbatch}); clear matlabbatch
    movefile(fullfile(pth,outputfolder,['s',sf,'_AvgR_native_unsmoothed.nii']),...
             fullfile(pth,outputfolder,[sf,'_AvgR_native.nii']));
    movefile(fullfile(pth,outputfolder,['s',sf,'_AvgR_Fz_native_unsmoothed.nii']),...
             fullfile(pth,outputfolder,[sf,'_AvgR_Fz_native.nii']));

    % warp back to MNI:

    options.coregmr.method='SPM'; % hard code for now

    copyfile(fullfile(pth,outputfolder,[sf,'_AvgR_Fz_native_unsmoothed.nii']),...
             fullfile(pth,outputfolder,[sf,'_AvgR_Fz.nii']));

         if strcmp(options.coregmr.method,'SPM') && exist([directory,'r',restfname,'2',ea_stripex(options.prefs.prenii_unnormalized),'_spm.mat'],'file')
             load([directory,'r',restfname,'2',ea_stripex(options.prefs.prenii_unnormalized),'_spm.mat']);
             nii=ea_load_nii(fullfile(pth,outputfolder,[sf,'_AvgR_Fz.nii']));
             nii.mat=spmaffine;
             ea_write_nii(nii);

             matlabbatch{1}.spm.spatial.coreg.write.ref = {[directory,options.prefs.prenii_unnormalized,',1']};
             matlabbatch{1}.spm.spatial.coreg.write.source = {nii.fname};
             matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0;
             matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
             matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
             matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
             spm_jobman('run',{matlabbatch});
             clear matlabbatch
             [pth,fn,ext]=fileparts(nii.fname);
             movefile(fullfile(pth,['r',fn,ext]),fullfile(pth,[fn,ext]));
         else
             % STILL NEED TO WRITE THIS:
             V = ea_load_nii([directory,restfname,'.nii,1']);
             V.fname=[directory,restfname,'_first_TR.nii'];
             ea_write_nii(V);
             ea_backuprestore([directory,restfname,'_first_TR.nii']);
             ea_coreg2images(options,[directory,restfname,'_first_TR.nii'],...
                 [directory,options.prefs.prenii_unnormalized],...
                 [directory,restfname,'_first_TR.nii'],...
                 {fullfile(pth,outputfolder,[sf,'_AvgR_Fz.nii'])});
             delete([directory,restfname,'_first_TR.nii']);
             delete([directory,'raw_',restfname,'_first_TR.nii']);

         end


    ea_apply_normalization_tofile(options,...
        {fullfile(pth,[sf,'_AvgR_Fz.nii'])},...
        {fullfile(pth,[sf,'_AvgR_Fz.nii'])},...
        directory,0,1,ea_niigz([ea_getearoot,'templates',filesep,'spacedefinitions',filesep,'222.nii']));

    nii=ea_load_nii(fullfile(pth,[sf,'_AvgR_Fz.nii']));
    nii.img(nii.img==0)=nan;
    nii.dt(2)=1;
    ea_write_nii(nii);
    matlabbatch{1}.spm.spatial.smooth.data = {fullfile(pth,[sf,'_AvgR_Fz.nii'])};
    matlabbatch{1}.spm.spatial.smooth.fwhm = [8 8 8];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 1;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    spm_jobman('run',{matlabbatch}); clear matlabbatch
    movefile(fullfile(pth,[sf,'_AvgR_Fz.nii']),fullfile(pth,[sf,'_AvgR_Fz_unsmoothed.nii']));
    movefile(fullfile(pth,['s',sf,'_AvgR_Fz.nii']),fullfile(pth,[sf,'_AvgR_Fz.nii']));
end


function sl=ea_detsiglength(fname)

V=spm_vol(fname);
sl=length(V);



function Result = rest_nextpow2_one35(n)
%Compute the min length for FFT according to AFNI's algorithm, By Xiao-Wei Song
%------------------------------------------------------------------------------------------------------------------------------
%	Copyright(c) 2007~2010
%	State Key Laboratory of Cognitive Neuroscience and Learning in Beijing Normal University
%	Written by Xiao-Wei Song
%	http://resting-fmri.sourceforge.net
% 	<a href="Dawnwei.Song@gmail.com">Mail to Author</a>: Xiaowei Song
%	Version=1.0;
%	Release=20070903;

if length(n)>1
    n = cast(length(n),class(n));
end
if n<16
    Result =2^nextpow2(n);
    return;
end

limit =nextpow2(n);             %n=134, limit=8
tbl=[2^(limit-1):2^limit];      %tbl =128, 129, ... , 256
tbl =tbl(find(tbl>=n));          %tbl =134, 135, ... , 256
for x=1:length(tbl)
    Result =tbl(x);
    [f,p]=log2(Result);
    if ~isempty(f) & f == 0.5   %Copy from nextpow2.m
        return;
    end
    if mod(Result,3*5)==0
        y= Result /(3*5);
        [f,p]=log2(y);
        if ~isempty(f) & f == 0.5   %Copy from nextpow2.m
            return;
        end
    end
    if mod(Result,3)==0
        y= Result /3;
        [f,p]=log2(y);
        if ~isempty(f) & f == 0.5   %Copy from nextpow2.m
            return;
        end
    end
    if mod(Result,5)==0
        y= Result /5;
        [f,p]=log2(y);
        if ~isempty(f) & f == 0.5   %Copy from nextpow2.m
            return;
        end
    end
end
Result =NaN;    % Should not reach, except when n=1

% csfft_nextup35 in AFNI list 1~1024, 20070516, dawnsong
% 2
% 4
% 6
% 8
% 10
% 12
% 16
% 20
% 24
% 30
% 32
% 40
% 48
% 60
% 64
% 80
% 96
% 120
% 128
% 160
% 192
% 240
% 256
% 320
% 384
% 480
% 512
% 640
% 768
% 960
% 1024

function preparecombinedvat(directory,stim)
% merge niftifiles i.e. with imcalc here.
