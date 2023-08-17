function seed_tc=cs_fmri_conseed_nifti(restfile,seedfile,options)

directory=[fileparts(restfile),filesep];
options=ea_getptopts(directory,options);
options=ea_assignpretra(options);

[~,restfname]=fileparts(restfile);
[~,anatfname]=fileparts(options.prefs.prenii_unnormalized);

options.prefs.rest=[restfname,'.nii'];

if ~exist([directory,'sr',restfname,'.nii'],'file') ...
        || ~exist([directory,'r',restfname,'_c1',anatfname,'.nii'],'file') % preproecessing needs to be performed
    disp('No preprocessed fMRI-images found, processing...');
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

%% Generate brain mask in fMRI space ('meanrest.nii')
if ~exist([directory,'mean',restfname,'_mask.nii'], 'file')
    copyfile([ea_getearoot,'templates',filesep,'spacedefinitions',filesep,'222.nii.gz'],...
        [directory,'temp.nii.gz']);
    gunzip([directory,'temp.nii.gz']);
    ea_delete([directory,'temp.nii.gz']);

    % Warp mask from MNI space to patient T1 space
    ea_apply_normalization_tofile(ea_getptopts(directory),...
        {[directory,'temp.nii']},...
        {[directory,'temp.nii']}, 1, 0);

    % Check coregistration method
    coregmethodsused = load([directory,'ea_coregmrmethod_applied.mat']);
    coregPrefix = ['r',restfname,'_',anatfname];
    if isfield(coregmethodsused, coregPrefix) && ~isempty(coregmethodsused.(coregPrefix))
        % Disable Hybrid coregistration
        coregmethod = strrep(coregmethodsused.(coregPrefix), 'Hybrid SPM & ', '');
        fprintf(['For this pair of coregistrations, the user specifically approved the ',coregmethod,' method.\n',...
            'Will overwrite the current global options and try to use this transform.\n']);
    else
        coregmethod = 'SPM'; % fallback to SPM coregistration
    end
    options.coregmr.method = coregmethod;

    % Check if the corresponding transform already exists
    xfm = [anatfname, '2r', restfname, '_', lower(coregmethod), '\d*\.(mat|h5)$'];
    transform = ea_regexpdir(directory, xfm, 0);

    if numel(transform) == 0
        warning('Transformation not found! Running coregistration now!');
        transform = ea_coregimages(options,[directory,options.prefs.prenii_unnormalized],...
            [directory,'mean',restfname,'.nii'],...
            [directory,'r',restfname,'_',options.prefs.prenii_unnormalized],...
            [],1,[],1);
        % Fix transformation names, replace 'mean' by 'r' for fMRI
        cellfun(@(f) movefile(f, strrep(f, 'mean', 'r')), transform);
        transform = strrep(transform, 'mean', 'r');
        transform = transform{1}; % Forward transformation
    else
        if numel(transform) > 1
            warning(['Multiple transformations of the same type found! ' ...
                'Will use the last one:\n%s'], transform{end});
        end
        transform = transform{end};
    end

    % Transform mask from patient T1 space to fMRI space ('meanrest.nii')
    ea_apply_coregistration([directory,'mean',restfname,'.nii'],...
        [directory,'temp.nii'],...
        [directory,'mean',restfname,'_mask.nii'],...
        transform,'nn');
    ea_delete([directory,'temp.nii']);
end

%% Extract timecourses of complete volume for signal regression
brainmask = ea_load_nii([directory,'mean',restfname,'_mask.nii']);
rest = ea_load_nii([directory,spfx,restfname,'.nii']);
signallength = size(rest.img,4);
interpol_tc = nan(numel(rest.img(:,:,:,1)),size(rest.img,4));
for tmpt = 1:signallength
    thisvol = rest.img(:,:,:,tmpt);
    thisvol = thisvol .* brainmask.img; % Mask out voxels outside the brain
    interpol_tc(:,tmpt) = thisvol(:);
end

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

signallength=size(rest.img,4);
for s=1:length(seedfile)
    seed{s}=ea_load_nii([seedfile{s}]);
    [xx,yy,zz]=ind2sub(size(seed{s}.img),1:numel(seed{s}.img));
    stringnum=cell(signallength,1);

    for i=1:signallength
        stringnum{i}=num2str(i);
    end
    single_s_files=cellfun(@(x) [directory,spfx,restfname,'.nii',',',x],stringnum,'Uniformoutput',false);
    single_s_files=single_s_files';
    V=spm_vol(single_s_files);

    nonzeros=find(seed{s}.img(:));
    vv=seed{s}.img(nonzeros);

    [xx,yy,zz]=ind2sub(size(seed{s}.img),nonzeros);

    voxelmask.locsvx=[xx,yy,zz,ones(size(xx,1),1)]';
    voxelmask.locsmm=[seed{s}.mat*voxelmask.locsvx]'; % get from voxels in parcellations to mm
    voxelmask.locsvx=[V{1}.mat\voxelmask.locsmm']'; % get from mm to voxels in restfile
    voxelmask.locsvx=voxelmask.locsvx(:,1:3);
    voxelmask.locsmm=voxelmask.locsmm(:,1:3);

    [allxx,allyy,allzz]=ind2sub(size(seed{s}.img),1:numel(seed{s}.img));

    allvoxelmask.locsvx=[allxx',allyy',allzz',ones(size(allxx,2),1)]';
    allvoxelmask.locsmm=[seed{s}.mat*allvoxelmask.locsvx]'; % get from voxels in parcellations to mm
    allvoxelmask.locsvx=[V{1}.mat\allvoxelmask.locsmm']'; % get from mm to voxels in restfile
    allvoxelmask.locsvx=allvoxelmask.locsvx(:,1:3);
    allvoxelmask.locsmm=allvoxelmask.locsmm(:,1:3);

    weights=vv./sum(vv);
    ea_dispercent(0,'Extracting time courses');
    clear seed_tc_all
    for i=1:signallength
        seed_tc_all(i,:)=spm_sample_vol(V{i},double(voxelmask.locsvx(:,1)),double(voxelmask.locsvx(:,2)),double(voxelmask.locsvx(:,3)),1);
        seed_tc{s}(i)=sum(seed_tc_all(i,:).*weights');
        ea_dispercent(i/signallength);
    end
    ea_dispercent(1,'end');

    disp('Done. Regressing out nuisance variables...');
end
interpol_tc=[cell2mat(seed_tc');interpol_tc];

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

voxelmask.locsvxi=ones(size(interpol_tc,1),1);
voxelmask.locsvxi(:)=1;
maskLowPass =	repmat(voxelmask.locsvxi, [1, paddedLength]);
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

% 20070513 remove trend --> FFT --> filter --> inverse FFT --> retrend
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

% cut seed_tc from voxel tc again:
for s=1:length(seedfile)
   seed_tc{s}=interpol_tc(s,:);
end
interpol_tc(1:length(seedfile),:)=[];


%% end  bandpass
disp('Done.');
for s=1:length(seedfile)
    %% export ROI map:
    if s>1
        if ~isequal(size(seed{s}),size(seed{s-1})) % this should not happen from within lead-dbs but potentially if people use the low-level function directly.
            ea_error('Seed-files have different size. Please supply these files separately!');
        end
    end
    %seed{s}.img(seed{s}.img==0)=nan;

    if ~isfield(options,'csfMRInowriteout')
        R=corr(seed_tc{s}',interpol_tc','rows','pairwise');
        expvol=rest.img(:,:,:,1);
        expvol(:)=R;
        interpvol=interp3(expvol,...
            allvoxelmask.locsvx(:,2),allvoxelmask.locsvx(:,1),allvoxelmask.locsvx(:,3));

        seed{s}.img(:)=interpvol;
        [pth,seedfn]=fileparts(seed{s}.fname);
        outputfolder = options.lcm.func.connectome;
        if strfind(seedfn,'rest')
            seedfn(strfind(seedfn,'rest')-1:end)=[];
        end
        seed{s}.fname=fullfile(pth,outputfolder,[seedfn,'_AvgR_native_unsmoothed.nii']);

        if ~exist(fullfile(pth,outputfolder),'dir')
            mkdir(fullfile(pth,outputfolder))
        end
        ea_write_nii(seed{s});
        seed{s}.dt(1)=16;
        seed{s}.img(:)=atanh(seed{s}.img(:));
        seed{s}.fname=fullfile(pth,outputfolder,[seedfn,'_AvgR_Fz_native_unsmoothed.nii']);
        ea_write_nii(seed{s});

        matlabbatch{1}.spm.spatial.smooth.data = {fullfile(pth,outputfolder,[seedfn,'_AvgR_native_unsmoothed.nii'])
            fullfile(pth,outputfolder,[seedfn,'_AvgR_Fz_native_unsmoothed.nii'])};
        matlabbatch{1}.spm.spatial.smooth.fwhm = [8 8 8];
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = 's';
        spm_jobman('run',{matlabbatch}); clear matlabbatch
        movefile(fullfile(pth,outputfolder,['s',seedfn,'_AvgR_native_unsmoothed.nii']),...
            fullfile(pth,outputfolder,[seedfn,'_AvgR_native.nii']));
        movefile(fullfile(pth,outputfolder,['s',seedfn,'_AvgR_Fz_native_unsmoothed.nii']),...
            fullfile(pth,outputfolder,[seedfn,'_AvgR_Fz_native.nii']));

        % warp back to MNI:

        copyfile(fullfile(pth,outputfolder,[seedfn,'_AvgR_Fz_native_unsmoothed.nii']),...
            fullfile(pth,outputfolder,[seedfn,'_AvgR_Fz.nii']));

        % Check coregistration method
        coregmethodsused = load([directory,'ea_coregmrmethod_applied.mat']);
        coregPrefix = ['r',restfname,'_',anatfname];
        if isfield(coregmethodsused, coregPrefix) && ~isempty(coregmethodsused.(coregPrefix))
            % Disable Hybrid coregistration
            coregmethod = strrep(coregmethodsused.(coregPrefix), 'Hybrid SPM & ', '');
            fprintf(['For this pair of coregistrations, the user specifically approved the ',coregmethod,' method.\n',...
                'Will overwrite the current global options and use this method.\n']);
        else
            coregmethod = 'SPM'; % fallback to SPM coregistration
        end
        options.coregmr.method = coregmethod;

        % Check if the transformation already exists
        xfm = ['hdmean', restfname, '2', anatfname, '_', lower(coregmethod), '\d*\.(mat|h5)$'];
        transform = ea_regexpdir(directory, xfm, 0);

        if numel(transform) == 0
            warning('Transformation not found! Running coregistration now!');
            transform = ea_coregimages(options,[directory,options.prefs.prenii_unnormalized],...
                [directory, 'hdmean', restfname, '.nii'],...
                [directory,'tmp.nii'],...
                [],1,[],1);
            ea_delete([directory,'tmp.nii']);
            transform = transform{2}; % Inverse transformation
        else
            if numel(transform) > 1
                warning(['Multiple transformations of the same type found! ' ...
                    'Will use the last one:\n%s'], transform{end});
            end
            transform = transform{end};
        end

        % Apply coregistration
        ea_apply_coregistration([directory,options.prefs.prenii_unnormalized], ...
            fullfile(pth,outputfolder,[seedfn,'_AvgR_Fz.nii']), ...
            fullfile(pth,outputfolder,[seedfn,'_AvgR_Fz.nii']), ...
            transform, 'linear');

        % Apply normalization
        ea_apply_normalization_tofile(options,...
            {fullfile(pth,outputfolder,[seedfn,'_AvgR_Fz.nii'])},...
            {fullfile(pth,outputfolder,[seedfn,'_AvgR_Fz.nii'])},...
            0,1,ea_niigz([ea_getearoot,'templates',filesep,'spacedefinitions',filesep,'222.nii']));

        nii=ea_load_nii(fullfile(pth,outputfolder,[seedfn,'_AvgR_Fz.nii']));
        nii.img(nii.img==0)=nan;
        nii.dt(2)=1;
        ea_write_nii(nii);
        matlabbatch{1}.spm.spatial.smooth.data = {fullfile(pth,outputfolder,[seedfn,'_AvgR_Fz.nii'])};
        matlabbatch{1}.spm.spatial.smooth.fwhm = [8 8 8];
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 1;
        matlabbatch{1}.spm.spatial.smooth.prefix = 's';
        spm_jobman('run',{matlabbatch}); clear matlabbatch
        movefile(fullfile(pth,outputfolder,[seedfn,'_AvgR_Fz.nii']),fullfile(pth,outputfolder,[seedfn,'_AvgR_Fz_unsmoothed.nii']));
        movefile(fullfile(pth,outputfolder,['s',seedfn,'_AvgR_Fz.nii']),fullfile(pth,outputfolder,[seedfn,'_AvgR_Fz.nii']));
    end
end


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
