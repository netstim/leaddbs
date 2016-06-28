function gmtc=ea_extract_timecourses(options)

% prepare voxelmask:

directory=[options.root,options.patientname,filesep];

ea_warp_parcellation(options.prefs.pprest,'rest',options);


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


%% Data corrections steps


disp('Calculating Global, WM and CSF-signals for signal regression...');

% regression steps
[~,rf]=fileparts(options.prefs.rest);
c1=ea_load_nii([directory,'rr',rf,'c1',options.prefs.prenii_unnormalized]);
c2=ea_load_nii([directory,'rr',rf,'c2',options.prefs.prenii_unnormalized]);
c3=ea_load_nii([directory,'rr',rf,'c3',options.prefs.prenii_unnormalized]);

globmap=logical((c1.img>0.5)+(c2.img)>0.5+(c3.img>0.5));
ec2map=c2.img; ec2map(ec2map<0.6)=0; ec2map=logical(ec2map);
ec3map=c3.img; ec3map(ec3map<0.6)=0; ec3map=logical(ec3map);





disp('Done. Regressing out nuisance variables...');
%% regress out movement parameters

rp_rest = load([directory,'rp_',rf,'.txt']); % rigid body motion parameters.
X(:,1)=ones(signallength,1);
X(:,2)=rp_rest(1:signallength,1);
X(:,3)=rp_rest(1:signallength,2);
X(:,4)=rp_rest(1:signallength,3);
X(:,5)=rp_rest(1:signallength,4);
X(:,6)=rp_rest(1:signallength,5);
X(:,7)=rp_rest(1:signallength,6);
%% actual regression:
for voxx=1:size(interpol_tc,1)
    beta_hat  = (X'*X)\X'*squeeze(interpol_tc(voxx,:))';
    if ~isnan(beta_hat)
    interpol_tc(voxx,:)=squeeze(interpol_tc(voxx,:))'-X*beta_hat;
    else
        warning('Regression of WM-/CSF-Signals could not be performed.');
    end 
end

% do the same on whole brain tc to get WM/GM/Global signal:
for xx=1:size(alltc,1)
    for yy=1:size(alltc,2)
        for zz=1:size(alltc,3)
            beta_hat  = (X'*X)\X'*squeeze(alltc(xx,yy,zz,:));
    if ~isnan(beta_hat)
    alltc(xx,yy,zz,:)=squeeze(alltc(xx,yy,zz,:))-X*beta_hat;
    else
        warning('Regression of WM-/CSF-Signals could not be performed.');
    end 
        end
    end
end


%% average Glob-, WM- and CSF-Timecourses
WMTimecourse=zeros(signallength,1);
CSFTimecourse=zeros(signallength,1);
GlobTimecourse=zeros(signallength,1);
for tmpt = 1:signallength
    OneTimePoint=alltc(:,:,:,tmpt);
    GlobTimecourse(tmpt)=squeeze(nanmean(OneTimePoint(ec2map(:))));
    WMTimecourse(tmpt)=squeeze(nanmean(OneTimePoint(ec2map(:))));
    CSFTimecourse(tmpt)=squeeze(nanmean(OneTimePoint(ec3map(:))));
end

%clear X
    X2(:,1)=ones(signallength,1);
    X2(:,2)=GlobTimecourse;
    X2(:,3)=WMTimecourse;
    X2(:,4)=CSFTimecourse;

    % regress out X from X2
    
    for tc=1:size(X2,2)
        
           beta_hat        = (X'*X)\X'*squeeze(X2(:,tc));
 
        X2(:,tc)=X2(:,tc)-X*beta_hat;
    end
    
    
%% actual regression:
for voxx=1:size(interpol_tc,1)

    beta_hat        = (X2'*X2)\X2'*squeeze(interpol_tc(voxx,:))';
    if ~isnan(beta_hat)
    interpol_tc(voxx,:)=squeeze(interpol_tc(voxx,:))'-X2*beta_hat;
    else
        warning('Regression of WM-/CSF-Signals could not be performed.');
    end
end

clear X X2



% for voxx=1:size(interpol_tc,1)
% 
%     beta_hat        = (X'*X)\X'*squeeze(interpol_tc(voxx,:))';
%     if ~isnan(beta_hat)
%     interpol_tc(voxx,:)=squeeze(interpol_tc(voxx,:))'-X*beta_hat;
%     else
%         warning('Regression of Motion parameters could not be performed.');
%     end
% 
% end


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
