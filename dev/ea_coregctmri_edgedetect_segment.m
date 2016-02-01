function varargout=ea_coregctmri_edgedetect_segment(options)
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ischar(options) % return name of method.
    varargout{1}='Coregister postop-CT with preop-MRI (Edgedetection incl. Segment)';
    varargout{2}={'SPM8','SPM12'};
    return
end

Dopt=1;
maxiter=200;



%% segment MRI


matlabbatch{1}.spm.tools.preproc8.channel.vols = {[options.root,options.patientname,filesep,options.prefs.prenii_unnormalized,',1']};
matlabbatch{1}.spm.tools.preproc8.channel.biasreg = 0.0001;
matlabbatch{1}.spm.tools.preproc8.channel.biasfwhm = 60;
matlabbatch{1}.spm.tools.preproc8.channel.write = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(1).tpm = {[spm('dir'),filesep,'toolbox',filesep,'Seg',filesep,'TPM.nii,1']};
matlabbatch{1}.spm.tools.preproc8.tissue(1).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(1).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(1).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(2).tpm = {[spm('dir'),filesep,'toolbox',filesep,'Seg',filesep,'TPM.nii,2']};
matlabbatch{1}.spm.tools.preproc8.tissue(2).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(2).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(2).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(3).tpm = {[spm('dir'),filesep,'toolbox',filesep,'Seg',filesep,'TPM.nii,3']};
matlabbatch{1}.spm.tools.preproc8.tissue(3).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(3).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(3).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(4).tpm = {[spm('dir'),filesep,'toolbox',filesep,'Seg',filesep,'TPM.nii,4']};
matlabbatch{1}.spm.tools.preproc8.tissue(4).ngaus = 3;
matlabbatch{1}.spm.tools.preproc8.tissue(4).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(4).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(5).tpm = {[spm('dir'),filesep,'toolbox',filesep,'Seg',filesep,'TPM.nii,5']};
matlabbatch{1}.spm.tools.preproc8.tissue(5).ngaus = 4;
matlabbatch{1}.spm.tools.preproc8.tissue(5).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(5).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(6).tpm = {[spm('dir'),filesep,'toolbox',filesep,'Seg',filesep,'TPM.nii,6']};
matlabbatch{1}.spm.tools.preproc8.tissue(6).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(6).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(6).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.warp.mrf = 0;
matlabbatch{1}.spm.tools.preproc8.warp.reg = 4;
matlabbatch{1}.spm.tools.preproc8.warp.affreg = 'mni';
matlabbatch{1}.spm.tools.preproc8.warp.samp = 3;
matlabbatch{1}.spm.tools.preproc8.warp.write = [0 0];
jobs{1}=matlabbatch;
if 0
cfg_util('run',jobs);
end
clear jobs matlabbatch



% generate "skull-contrast" image
matlabbatch{1}.spm.util.imcalc.input = {[options.root,options.patientname,filesep,'c1',options.prefs.prenii_unnormalized];
    [options.root,options.patientname,filesep,'c2',options.prefs.prenii_unnormalized];
    [options.root,options.patientname,filesep,'c3',options.prefs.prenii_unnormalized];
[options.root,options.patientname,filesep,'c4',options.prefs.prenii_unnormalized];
[options.root,options.patientname,filesep,'c5',options.prefs.prenii_unnormalized];
[options.root,options.patientname,filesep,'c6',options.prefs.prenii_unnormalized]};
matlabbatch{1}.spm.util.imcalc.output = ['skullcon',options.prefs.prenii_unnormalized];
matlabbatch{1}.spm.util.imcalc.outdir = {[options.root,options.patientname,filesep]};
matlabbatch{1}.spm.util.imcalc.expression = 'i4-(i1+i2+i3+i5+i6)';
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    jobs{1}=matlabbatch;
    cfg_util('run',jobs);
    clear matlabbatch jobs;


% use skull volume here:
ea_reslice_nii([options.root,options.patientname,filesep,'skullcon',options.prefs.prenii_unnormalized],[options.root,options.patientname,filesep,'small_',options.prefs.prenii_unnormalized],[4 4 4]);
ea_reslice_nii([options.root,options.patientname,filesep,options.prefs.rawctnii_unnormalized],[options.root,options.patientname,filesep,'small_',options.prefs.rawctnii_unnormalized],[4 4 4]);


vizz=1; % visualization on
alphas=options.coregct.coregthreshs;
%% estimate
disp('Loading images...');
CT=ea_load_nii([options.root,options.patientname,filesep,'small_',options.prefs.rawctnii_unnormalized],'simple');
MR=ea_load_nii([options.root,options.patientname,filesep,'small_',options.prefs.prenii_unnormalized],'simple');
delete([options.root,options.patientname,filesep,'small_',options.prefs.rawctnii_unnormalized]);
delete([options.root,options.patientname,filesep,'small_',options.prefs.prenii_unnormalized]);

disp('Done. Smoothing...');
MR.img(MR.img<100)=0; % remove noise around brain.
MR.img=smooth3(MR.img,'gaussian',[11 11 11]);
CT.img(CT.img<0)=0; % remove negative hounsfield parts.
CT.img=smooth3(CT.img,'gaussian',[11 11 11]);

if vizz
    h=figure('color','w','name','Coregistering CT to MR...','NumberTitle','off');
    axis equal
    axis off
end

disp('Done.');

for alpha=1:length(alphas)

    disp(['*** Pass ',num2str(alpha),'/',num2str(length(alphas)),', alpha: ',num2str(alphas(alpha)),'.']);
    eMR=ea_detect_edges_3d(MR.img,alphas(alpha));
    eCT=ea_detect_edges_3d(CT.img,alphas(alpha));

    eMR(eMR<2)=0;
    eCT(eCT<2)=0;

    [xx,yy,zz]=ind2sub(size(eMR),find(eMR(:)));
    ptMR=[xx,yy,zz,ones(length(xx),1)]';
    ptMR=MR.hdr.mat*ptMR; % mm notation
    [xx,yy,zz]=ind2sub(size(eCT),find(eCT(:)));
    ptCT=[xx,yy,zz,ones(length(xx),1)]';
    ptCT=CT.hdr.mat*ptCT; % mm notation
    disp('Done.');


    % define initialization parameters

    MRTree=KDTreeSearcher(ptMR');
    if ~exist('M','var') % first grain run..
        % define M based on centroids
        M=eye(4);
        for dim=1:3
            M(dim,4)=mean(ptMR(dim,:))-mean(ptCT(dim,:));
        end
    end
    ptrCT=M*ptCT;
    [~,DCT]=knnsearch(MRTree,ptrCT');
    [~,DMR]=knnsearch(ptrCT',ptMR');
    D=mean([DCT(:);DMR(:)]);
    priorD=D;
    reuselastmodM=0;
    cnt=0;
    priorM=M; % reset params

     if Dopt>D
        Dopt=D/2;
    end
    while cnt<maxiter && D>Dopt % exit conditions (hard coded for now).

        % modify M
        if ~reuselastmodM
            modM=randn(3,4)*(randn(1)*(0.5/alpha));% *0.01;
        end
        M=M+[modM;zeros(1,4)]; % fuzzy improvement
        ptrCT=M*ptCT;


        [~,DCT]=knnsearch(MRTree,ptrCT');
        [~,DMR]=knnsearch(ptrCT',ptMR');
        D=mean([DCT(:);DMR(:)]);



        % check improvement.
        if D<priorD % improvement
            fprintf('\n %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f \n\n',[D;M(:)]');
            if vizz
                plot3(ptMR(1,:),ptMR(2,:),ptMR(3,:),'.','Color',[0.8,0.2,0.8]);
                hold on
                plot3(ptrCT(1,:),ptrCT(2,:),ptrCT(3,:),'.','Color',[0.2,0.8,0.8]);
                hold off
                drawnow
            end

            priorD=D; % update priors

            priorM=M;
            reuselastmodM=1;
            cnt=0;
        else % no improvement -> reset.
            %         fprintf('\n %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f \n',[D;M(:)]');
            %         if vizz
            %             plot3(ptMR(1,:),ptMR(2,:),ptMR(3,:),'g.');
            %             hold on
            %             plot3(ptCT(1,:),ptCT(2,:),ptCT(3,:),'r.');
            %             hold off
            %             drawnow
            %         end
            %disp(['D = ',num2str(D),'.']);
            %disp(['priorD = ',num2str(priorD),'.']);
            cnt=cnt+1;
            M=priorM; % reset M to last value.
            reuselastmodM=0;
        end


    end

    if D<Dopt % since this is a hard exit criterion, it wouldn't help to re-iterate.
        break
    end


end


% M has been estimated and maps from voxels in CT to voxels in MR.
%% export coregistered CT.

matlabbatch{1}.spm.util.reorient.srcfiles = {[options.root,options.patientname,filesep,options.prefs.rawctnii_unnormalized,',1']};
matlabbatch{1}.spm.util.reorient.transform.transM = M;
matlabbatch{1}.spm.util.reorient.prefix = 'r';
jobs{1}=matlabbatch;
cfg_util('run',jobs);
clear jobs matlabbatch


%% use spm coreg to finalize result.

% we might want to try normalizing the rCT using Albada 2009?
% disp('Done. Normalizing CT intensities...');
% rCT=ea_load_nii(
% CT.img=reshape(ea_normal(CT.img(:)),size(CT.img,1),size(CT.img,2),size(CT.img,3));

% 1. initial coreg
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[options.root,options.patientname,filesep,'',options.prefs.prenii_unnormalized]};
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {[options.root,options.patientname,filesep,'r',options.prefs.rawctnii_unnormalized,',1']};
    matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

    jobs{1}=matlabbatch;
    cfg_util('run',jobs);
    clear matlabbatch jobs;

    % 2.
% multiply both CT and MRI with c1-3 volume of images:
% skullstrip MR:
matlabbatch{1}.spm.util.imcalc.input = {[options.root,options.patientname,filesep,'c1',options.prefs.prenii_unnormalized];
    [options.root,options.patientname,filesep,'c2',options.prefs.prenii_unnormalized];
    [options.root,options.patientname,filesep,'c3',options.prefs.prenii_unnormalized];
    [options.root,options.patientname,filesep,options.prefs.prenii_unnormalized]};
matlabbatch{1}.spm.util.imcalc.output = ['b',options.prefs.prenii_unnormalized];
matlabbatch{1}.spm.util.imcalc.outdir = {[options.root,options.patientname,filesep]};
matlabbatch{1}.spm.util.imcalc.expression = '(i1+i2+i3).*i4';
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    jobs{1}=matlabbatch;
    cfg_util('run',jobs);
    clear matlabbatch jobs;

    % skullstrip CT:
matlabbatch{1}.spm.util.imcalc.input = {[options.root,options.patientname,filesep,'c1',options.prefs.prenii_unnormalized];
    [options.root,options.patientname,filesep,'c2',options.prefs.prenii_unnormalized];
    [options.root,options.patientname,filesep,'c3',options.prefs.prenii_unnormalized];
    [options.root,options.patientname,filesep,'r',options.prefs.rawctnii_unnormalized]};
matlabbatch{1}.spm.util.imcalc.output = ['br',options.prefs.rawctnii_unnormalized];
matlabbatch{1}.spm.util.imcalc.outdir = {[options.root,options.patientname,filesep]};
matlabbatch{1}.spm.util.imcalc.expression = '(i1+i2+i3).*i4';
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    jobs{1}=matlabbatch;
    cfg_util('run',jobs);
    clear matlabbatch jobs;


    % 3. Final coreg.
costfuns={'nmi','mi','ecc'};
for costfun=1:3
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[options.root,options.patientname,filesep,'b',options.prefs.prenii_unnormalized,',1']};
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {[options.root,options.patientname,filesep,'br',options.prefs.rawctnii_unnormalized,',1']};
    matlabbatch{1}.spm.spatial.coreg.estimate.other = {[options.root,options.patientname,filesep,'r',options.prefs.rawctnii_unnormalized,',1']};
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = costfuns{costfun};
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

    jobs{1}=matlabbatch;
    cfg_util('run',jobs);
    clear matlabbatch jobs;
end

disp('Cleaning up...');
delete([options.root,options.patientname,filesep,'br',options.prefs.rawctnii_unnormalized]);
delete([options.root,options.patientname,filesep,'b',options.prefs.prenii_unnormalized]);
delete([options.root,options.patientname,filesep,'skullcon',options.prefs.prenii_unnormalized]);

disp('Done.');

matlabbatch{1}.spm.util.checkreg.data = {[options.root,options.patientname,filesep,options.prefs.prenii_unnormalized];
    [options.root,options.patientname,filesep,'r',options.prefs.rawctnii_unnormalized,',1']};
jobs{1}=matlabbatch;
cfg_util('run',jobs);
clear matlabbatch jobs;
