function [gmtc,interpol_tc,voxelmask]=ea_extract_timecourses(options)
% Extract timecourses based on template

%% prepare voxelmask:
directory=[options.root,options.patientname,filesep];
ea_warp_parcellation(['r', options.prefs.rest], options);
vizz=0;

%% create voxelmask
[~,rrf]=fileparts(options.prefs.rest);
Vatl=spm_vol([directory,'templates',filesep,'labeling',filesep,'r',rrf,'w',options.lc.general.parcellation,'.nii']);
Xatl=spm_read_vols(Vatl);
Xatl(isnan(Xatl))=0;
nonzeros=find(Xatl(:));
vv=Xatl(nonzeros);

[xx,yy,zz]=ind2sub(size(Xatl),nonzeros);
restfilename=options.prefs.pprest;
signallength=ea_detsiglength([directory,restfilename]);
stringnum=cell(signallength,1);

for i=1:signallength
    stringnum{i}=num2str(i);
end
single_s_files=cellfun(@(x) [directory,restfilename,',',x],stringnum,'Uniformoutput',false);
single_s_files=single_s_files';
V=spm_vol(single_s_files);
voxelmask.locsvx=[xx,yy,zz,ones(size(xx,1),1)]';
voxelmask.locsmm=[Vatl.mat*voxelmask.locsvx]'; % get from voxels in parcellations to mm
voxelmask.locsvx=[V{1}.mat\voxelmask.locsmm']'; % get from mm to voxels in restfile
voxelmask.locsvx=voxelmask.locsvx(:,1:3);
voxelmask.locsmm=voxelmask.locsmm(:,1:3);
voxelmask.vals=round(vv);
todel=[];
vmaskentries=unique(voxelmask.vals)';
for entry=vmaskentries
    if sum(voxelmask.vals==entry)>10000 % normally sized ROI, can downsample for speed
    thisdelete=find(voxelmask.vals==entry);
    thisdelete=thisdelete(1:2:end);
    todel=[todel;thisdelete];
    end
end
voxelmask.locsmm(todel,:)=[];
voxelmask.locsvx(todel,:)=[];
voxelmask.vals(todel)=[];

%% set some initial parameters here:
TR=options.lc.func.prefs.TR;
save([directory,'TR.mat'],'TR');


%% Extract timecourses of specified ROI
% ea_dispercent(0,'Extracting time courses');
parfor i=1:signallength
    disp(['Extracting time courses: ', num2str(i,'%03d'),'/',num2str(signallength)]);
    interpol_tc(i,:)=spm_sample_vol(V{i},double(voxelmask.locsvx(:,1)),double(voxelmask.locsvx(:,2)),double(voxelmask.locsvx(:,3)),-2);
%     ea_dispercent(i/signallength);
end
% ea_dispercent(1,'end');
interpol_tc=interpol_tc';

disp('Done.');

%% Load complete volume for signal regression..
rfile=[directory,restfilename];
alltc=spm_read_vols(spm_vol(rfile));

% % detrend global signal
% nDim1 = size(alltc,1); nDim2 = size(alltc,2); nDim3 = size(alltc,3); nDim4 =size(alltc,4);
%
% for x=1:nDim1
%     oneslice =double(alltc(x, :, :, :));
%     oneslice =reshape(oneslice, 1*nDim2*nDim3, nDim4)';
%     oneslice =detrend(oneslice) +repmat(mean(oneslice), [size(oneslice,1), 1]);
%     oneslice =reshape(oneslice', 1,nDim2,nDim3, nDim4);
%     alltc(x, :, :, :) =(oneslice);
% end;

[~,restbase]=fileparts(restfilename);
if exist([directory,restbase(3:end),'_sessvec.mat'],'file')
    multsess=1;
    load([directory,restbase(3:end),'_sessvec.mat']);
else
    multsess=0;
    sessvec=ones(size(interpol_tc,2),1);
end

sessvec=ea_gensessvec(sessvec);

if vizz
    figure
    pcnt=1;
    subplot(3,2,pcnt)
    pcnt=pcnt+1;
    plot(interpol_tc(round(1:size(interpol_tc,1)/1000:size(interpol_tc,1)),:)');
    title('Raw timeseries')
end

%% Data corrections steps
disp('Regressing out nuisance variables...');
if multsess
    if vizz
        subplot(3,2,pcnt)
        pcnt=pcnt+1;
        plot(sessvec);
        title('Session vector')
    end

    % regress out sessions
    X0(:,1)=ones(signallength,1);
    X0=[X0,sessvec+1];

    % actual regression:
    X0reg=(X0'*X0)\X0';
    for voxx=1:size(interpol_tc,1)
        beta_hat  = X0reg*squeeze(interpol_tc(voxx,:))';
        if ~any(isnan(beta_hat))
            interpol_tc(voxx,:)=squeeze(interpol_tc(voxx,:))'-X0*beta_hat;
        else

            warning('Regression of WM-/CSF-Signals could not be performed.');
        end
    end

    if vizz
        subplot(3,2,pcnt)
        pcnt=pcnt+1;
        plot(interpol_tc(round(1:size(interpol_tc,1)/1000:size(interpol_tc,1)),:)');
        title('Time series cleaned from session vector.')
    end

    % do the same on whole brain tc to get WM/GM/Global signal:
    for xx=1:size(alltc,1)
        for yy=1:size(alltc,2)
            for zz=1:size(alltc,3)
                beta_hat  = X0reg*squeeze(alltc(xx,yy,zz,:));
                if ~any(isnan(beta_hat))
                    alltc(xx,yy,zz,:)=squeeze(alltc(xx,yy,zz,:))-X0*beta_hat;
                else
                    warning('Regression of WM-/CSF-Signals could not be performed.');
                end
            end
        end
    end
end

% average Glob-, WM- and CSF-Timecourses
disp('Calculating Global, WM and CSF-signals for signal regression...');

[~,rf]=fileparts(options.prefs.rest);
% regression steps
c1=ea_load_nii([directory,'r',rf,'_c1',options.prefs.prenii_unnormalized]);
c2=ea_load_nii([directory,'r',rf,'_c2',options.prefs.prenii_unnormalized]);
c3=ea_load_nii([directory,'r',rf,'_c3',options.prefs.prenii_unnormalized]);

globmap=logical((c1.img>0.5)+(c2.img)>0.5+(c3.img>0.5));
ec2map=c2.img; ec2map(ec2map<0.6)=0; ec2map=logical(ec2map);
ec3map=c3.img; ec3map(ec3map<0.6)=0; ec3map=logical(ec3map);

WMTimecourse=zeros(signallength,1);
CSFTimecourse=zeros(signallength,1);
GlobTimecourse=zeros(signallength,1);
parfor tmpt = 1:signallength
    OneTimePoint=alltc(:,:,:,tmpt);
    try
        GlobTimecourse(tmpt)=squeeze(nanmean(OneTimePoint(globmap(:))));
        WMTimecourse(tmpt)=squeeze(nanmean(OneTimePoint(ec2map(:))));
        CSFTimecourse(tmpt)=squeeze(nanmean(OneTimePoint(ec3map(:))));
    catch
        keyboard
    end
end

if vizz
    subplot(3,2,pcnt)
    pcnt=pcnt+1;
    plot([GlobTimecourse,WMTimecourse,CSFTimecourse]);
    title('Global/WM/CSF Timecourses (cleaned from session).');
end

% regress out movement parameters
rp_rest = load([directory,'rp_',rf,'.txt']); % rigid body motion parameters.
X1(:,1)=ones(signallength,1);
X1(:,2)=rp_rest(1:signallength,1);
X1(:,3)=rp_rest(1:signallength,2);
X1(:,4)=rp_rest(1:signallength,3);
X1(:,5)=rp_rest(1:signallength,4);
X1(:,6)=rp_rest(1:signallength,5);
X1(:,7)=rp_rest(1:signallength,6);

if multsess
    % regress sessions from movement parameters
    for mov=2:size(X1,2)
        beta_hat  = X0reg*squeeze(X1(:,mov));
        if ~any(isnan(beta_hat))
            X1(:,mov)=squeeze(X1(:,mov))-X0*beta_hat;
        else
            warning('Regression of motion parameters could not be performed.');
        end
    end
end

if vizz
    subplot(3,2,pcnt)
    pcnt=pcnt+1;
    plot(X1);
    title('Motion parameters (cleaned from session vector).');
end

% actual regression:
X1reg=(X1'*X1)\X1';
for voxx=1:size(interpol_tc,1)
    beta_hat  = X1reg*squeeze(interpol_tc(voxx,:))';
    if ~any(isnan(beta_hat))
        interpol_tc(voxx,:)=squeeze(interpol_tc(voxx,:))'-X1*beta_hat;
    else
        warning('Regression of motion parameters could not be performed.');
    end
end

if vizz
    subplot(3,2,pcnt)
    pcnt=pcnt+1;
    plot(interpol_tc(round(1:size(interpol_tc,1)/1000:size(interpol_tc,1)),:)');
    title('Time series cleaned from motion parameters.');
end

X2(:,1)=ones(signallength,1);

if options.prefs.lc.func.regress_wmcsf
    X2(:,2) = WMTimecourse;
    X2(:,3) = CSFTimecourse;
end

if options.prefs.lc.func.regress_global
    X2(:,4) = GlobTimecourse;
end

% actual regression of cleaned X2 (WM/Global) from time courses:
X2reg=(X2'*X2)\X2';
for voxx=1:size(interpol_tc,1)

    beta_hat = X2reg*squeeze(interpol_tc(voxx,:))';
    if ~any(isnan(beta_hat))
        interpol_tc(voxx,:)=squeeze(interpol_tc(voxx,:))'-X2*beta_hat;
    else
        warning('Regression of WM-/CSF-Signals could not be performed.');
    end
end

clear X X2
disp('Done.')

if vizz
    subplot(3,2,pcnt)
    pcnt=pcnt+1;
    plot(interpol_tc(round(1:size(interpol_tc,1)/1000:size(interpol_tc,1)),:)');
    title('Time series (cleaned from wm&csf).');
    % title('Time series (cleaned from wm&csf&global).');
end

%% Bandpass filtering
lp_HighCutoff=options.prefs.lc.func.bphighcutoff;
hp_LowCutoff=options.prefs.lc.func.bplowcutoff;

disp('Bandpass-filtering...');
% sampleFreq   = 1/TR;
sampleLength = signallength;
paddedLength = rest_nextpow2_one35(sampleLength); %2^nextpow2(sampleLength);

mask=ones(size(interpol_tc,1),1);
mask(:)=1;
maskLowPass =	repmat(mask, [1, paddedLength]);
maskHighPass=	maskLowPass;

% GENERATE LOW PASS WINDOW	20070514, reference: fourior_filter.c in AFNI
% Revised by YAN Chao-Gan, 100420. Fixed the bug in calculating the frequency band.
% Low pass, such as freq < 0.08 Hz
idxCutoff	=round(lp_HighCutoff *paddedLength *TR + 1); % Revised by YAN Chao-Gan, 100420. Fixed the bug in calculating the frequency band. %idxCutoff	=round(ALowPass_HighCutoff *paddedLength *TR);
idxCutoff2	=paddedLength+2 -idxCutoff;				%Center Index =(paddedLength/2 +1)
maskLowPass(:,idxCutoff+1:idxCutoff2-1)=0; %High eliminate

% GENERATE HIGH PASS WINDOW
% high pass, such as freq > 0.01 Hz
idxCutoff	=round(hp_LowCutoff *paddedLength *TR + 1); % Revised by YAN Chao-Gan, 100420. Fixed the bug in calculating the frequency band. %idxCutoff	=round(AHighPass_LowCutoff *paddedLength *TR);
idxCutoff2	=paddedLength+2 -idxCutoff;				%Center Index =(paddedLength/2 +1)
maskHighPass(:,1:idxCutoff-1)=0;	%Low eliminate
maskHighPass(:,idxCutoff2+1:paddedLength)=0;	%Low eliminate

% 20070513	remove trend --> FFT --> filter --> inverse FFT --> retrend
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

disp('Done.');

if vizz
    subplot(3,2,pcnt)
    plot(interpol_tc(round(1:size(interpol_tc,1)/1000:size(interpol_tc,1)),:)');
    title('Bandpass filtered time series.');
end

%% average gmtc over ROI
aID = fopen([ea_space(options,'labeling'),options.lc.general.parcellation,'.txt']);
atlas_lgnd=textscan(aID,'%d %s');
dimensionality=length(atlas_lgnd{1}); % how many ROI.

gmtc=nan(size(interpol_tc,2),dimensionality);
cnt=1;
for c=double(atlas_lgnd{1}')
    gmtc(:,cnt)=nanmean(interpol_tc(voxelmask.vals==c,:),1);
    cnt=cnt+1;
end


%% add methods dump:
cits={
    'Horn, A., Ostwald, D., Reisert, M., & Blankenburg, F. (2014). The structural-functional connectome and the default mode network of the human brain. NeuroImage, 102 Pt 1, 142?151. http://doi.org/10.1016/j.neuroimage.2013.09.069'
    'Weissenbacher, A., Kasess, C., Gerstl, F., Lanzenberger, R., Moser, E., & Windischberger, C. (2009). Correlations and anticorrelations in resting-state functional connectivity MRI: a quantitative comparison of preprocessing strategies., 47(4), 1408?1416. http://doi.org/10.1016/j.neuroimage.2009.05.005'
    'Horn, A., & Kuehn, A. A. (2015). Lead-DBS: a toolbox for deep brain stimulation electrode localizations and visualizations. NeuroImage, 107, 127?135. http://doi.org/10.1016/j.neuroimage.2014.12.002'
    'Horn, A., Li, N., Dembek, T. A., Kappel, A., Boulay, C., Ewert, S., et al. (2019). Lead-DBS v2: Towards a comprehensive pipeline for deep brain stimulation imaging. NeuroImage, 184, 293?316. http://doi.org/10.1016/j.neuroimage.2018.08.068'
    };
ea_methods(options,['Resting-state fMRI data was preprocessed following the pipeline described in (Horn et al. 2014) as implemented in Lead-DBS software (Horn & K?hn 2015; Horn & Li et al. 2018; www.lead-dbs.org). Pre-processing steps broadly follow the suggestions made in ',...
    ' (Weissenbacher 2009). This involved realignment of data, regression of movement-parameters, a WM-, CSF- as well as global signal and band-pass filtering (',num2str(options.prefs.lc.func.bphighcutoff),'-',...
    num2str(options.prefs.lc.func.bplowcutoff),' Hz). No spatial smoothing was applied.'],...
    cits);


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


function sessvec=ea_gensessvec(sessvec)
X=zeros(size(sessvec,1),max(sessvec));
for sess=1:max(sessvec)
    X(:,sess)=sessvec==sess;
end
sessvec=X;
