function gmtc=ea_extract_timecourses(options)

% prepare voxelmask:

directory=[options.root,options.patientname,filesep];

ea_warp_parcellation(options.prefs.pprest,'rest',options);
vizz=0;

%% create voxelmask

Vatl=spm_vol([directory,'templates',filesep,'labeling',filesep,'rrestw',options.lc.general.parcellation,'.nii,1']);
Xatl=spm_read_vols(Vatl);

nonzeros=find(Xatl(:));
vv=Xatl(nonzeros);

[xx,yy,zz]=ind2sub(size(Xatl),nonzeros);

voxelmask.locsvx=[xx,yy,zz,ones(size(xx,1),1)]';
voxelmask.locsmm=[Vatl.mat*voxelmask.locsvx]';
voxelmask.locsvx=voxelmask.locsvx(1:3,:)';
voxelmask.locsmm=voxelmask.locsmm(:,1:3);
voxelmask.vals=round(vv);

%% set some initial parameters here:

TR=options.lc.func.prefs.TR;
directory=[options.root,options.patientname,filesep];
restfilename=options.prefs.pprest;
signallength=ea_detsiglength([directory,restfilename]);
stringnum=cell(signallength,1);

for i=1:signallength
    stringnum{i}=num2str(i);
end
single_s_files=cellfun(@(x) [directory,restfilename,',',x],stringnum,'Uniformoutput',false);
single_s_files=single_s_files';


%% Extract timecourses of specified ROI
V=spm_vol(single_s_files);

for i=1:signallength
    interpol_tc(i,:)=spm_sample_vol(V{i},double(voxelmask.locsvx(:,1)),double(voxelmask.locsvx(:,2)),double(voxelmask.locsvx(:,3)),1);

end

aID = fopen([options.earoot,'templates',filesep,'labeling',filesep,options.lc.general.parcellation,'.txt']);
atlas_lgnd=textscan(aID,'%d %s');
dimensionality=length(atlas_lgnd{1}); % how many ROI.


%% Extract timecourses of complete volume for signal regression..
rfile=[directory,restfilename];
alltc=spm_read_vols(spm_vol(rfile));

interpol_tc=interpol_tc';

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
     load([directory,restbase(3:end),'_sessvec.mat']);
 else
     sessvec=ones(size(interpol_tc,2),1);
 end
 
 sessvec=ea_gensessvec(sessvec);
 
if vizz
    figure
    subplot(4,2,1)
    plot(interpol_tc(round(1:size(interpol_tc,1)/1000:size(interpol_tc,1)),:)');
    title('Raw timeseries')
end
if vizz
    subplot(4,2,2)
    plot(sessvec);
    title('Session vector')
end
%% Data corrections steps








disp('Done. Regressing out nuisance variables...');
%% regress out sessions
X0(:,1)=ones(signallength,1);


X0=[X0,sessvec+1];

%% actual regression:
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
    subplot(4,2,3)
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


%% average Glob-, WM- and CSF-Timecourses
disp('Calculating Global, WM and CSF-signals for signal regression...');

[~,rf]=fileparts(options.prefs.rest);
% regression steps
c1=ea_load_nii([directory,'rr',rf,'c1',options.prefs.prenii_unnormalized]);
c2=ea_load_nii([directory,'rr',rf,'c2',options.prefs.prenii_unnormalized]);
c3=ea_load_nii([directory,'rr',rf,'c3',options.prefs.prenii_unnormalized]);

globmap=logical((c1.img>0.5)+(c2.img)>0.5+(c3.img>0.5));
ec2map=c2.img; ec2map(ec2map<0.6)=0; ec2map=logical(ec2map);
ec3map=c3.img; ec3map(ec3map<0.6)=0; ec3map=logical(ec3map);

WMTimecourse=zeros(signallength,1);
CSFTimecourse=zeros(signallength,1);
GlobTimecourse=zeros(signallength,1);
for tmpt = 1:signallength
    OneTimePoint=alltc(:,:,:,tmpt);
    GlobTimecourse(tmpt)=squeeze(nanmean(OneTimePoint(globmap(:))));
    WMTimecourse(tmpt)=squeeze(nanmean(OneTimePoint(ec2map(:))));
    CSFTimecourse(tmpt)=squeeze(nanmean(OneTimePoint(ec3map(:))));
end

if vizz
        subplot(4,2,4)
    plot([GlobTimecourse,WMTimecourse,CSFTimecourse]);
    title('Global/WM/CSF Timecourses (cleaned from session).');
end


%% regress out movement parameters

rp_rest = load([directory,'rp_',rf,'.txt']); % rigid body motion parameters.
X1(:,1)=ones(signallength,1);
X1(:,2)=rp_rest(1:signallength,1);
X1(:,3)=rp_rest(1:signallength,2);
X1(:,4)=rp_rest(1:signallength,3);
X1(:,5)=rp_rest(1:signallength,4);
X1(:,6)=rp_rest(1:signallength,5);
X1(:,7)=rp_rest(1:signallength,6);
X1(:,7)=WMTimecourse;
X1(:,7)=CSFTimecourse;
X1(:,7)=GlobTimecourse;

% regress sessions from movement parameters
for mov=2:size(X1,2)
    beta_hat  = X0reg*squeeze(X1(:,mov));
    if ~any(isnan(beta_hat))
    X1(:,mov)=squeeze(X1(:,mov))-X0*beta_hat;
    else
        warning('Regression of WM-/CSF-Signals could not be performed.');
    end 
end

if vizz
    subplot(4,2,5)
    plot(X1);
    title('Motion parameters (cleaned from session vector).');
end

%% actual regression:
X1reg=(X1'*X1)\X1';
for voxx=1:size(interpol_tc,1)
    beta_hat  = X1reg*squeeze(interpol_tc(voxx,:))';
    if ~any(isnan(beta_hat))
    interpol_tc(voxx,:)=squeeze(interpol_tc(voxx,:))'-X1*beta_hat;
    else
        warning('Regression of WM-/CSF-Signals could not be performed.');
    end 
end

if vizz
    subplot(4,2,6)
    plot(interpol_tc(round(1:size(interpol_tc,1)/1000:size(interpol_tc,1)),:)');
    title('Time series cleaned from motion parameters.');
end



X2(:,1)=ones(signallength,1);
if options.prefs.lc.func.regress_wmcsf
    X2(:,2)=WMTimecourse;
    X2(:,3)=CSFTimecourse;
end
if options.prefs.lc.func.regress_global
    X2(:,4)=GlobTimecourse;
end    
    
%% actual regression of cleaned X2 (WM/Global) from time courses:
X2reg=(X2'*X2)\X2';
for voxx=1:size(interpol_tc,1)

    beta_hat        = X2reg*squeeze(interpol_tc(voxx,:))';
    if ~any(isnan(beta_hat))
    interpol_tc(voxx,:)=squeeze(interpol_tc(voxx,:))'-X2*beta_hat;
    else
        warning('Regression of WM-/CSF-Signals could not be performed.');
    end
end

if vizz
    subplot(4,2,7)
    plot(interpol_tc(round(1:size(interpol_tc,1)/1000:size(interpol_tc,1)),:)');
    title('Time series (cleaned from global/csf/wm).');
end

clear X X2


%% begin rest bandpass

lp_HighCutoff=options.prefs.lc.func.bphighcutoff;
hp_LowCutoff=options.prefs.lc.func.bplowcutoff;


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

if vizz
    subplot(4,2,8)
    plot(interpol_tc(round(1:size(interpol_tc,1)/1000:size(interpol_tc,1)),:)');
    title('Bandpass filtered time series.');
end


%% average gmtc over ROI

gmtc=nan(size(interpol_tc,2),dimensionality);

for c=1:dimensionality
    gmtc(:,c)=mean(interpol_tc(voxelmask.vals==c,:));
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


function sessvec=ea_gensessvec(sessvec)

X=zeros(size(sessvec,1),max(sessvec));
for sess=1:max(sessvec)
    X(:,sess)=sessvec==sess;
end
sessvec=X;
