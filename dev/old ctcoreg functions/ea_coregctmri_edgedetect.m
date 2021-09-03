function varargout=ea_coregctmri_edgedetect(options)
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ischar(options) % return name of method.
    varargout{1}='Edgedetection';
    varargout{2}={'SPM8','SPM12'};
    varargout{3}=['0.6,0.4']; % suggestion for alpha-parameter.
    return
end


Dopt=12;
maxiter=200;

    disp('Loading images...');

%ea_reslice_nii([options.subj.preopAnat.(options.subj.AnchorModality).coreg],[options.subj.preopAnat.(options.subj.AnchorModality).coreg],[0.5 0.5 0.5]);

% MR
if isfield(options,'usediffmr_coregct')
    ea_reslice_nii([options.root,options.patientname,filesep,options.usediffmr_coregct],[options.root,options.patientname,filesep,'small_',options.usediffmr_coregct],[2 2 2],0);
    MR=ea_load_nii([options.root,options.patientname,filesep,'small_',options.usediffmr_coregct]);
    delete([options.root,options.patientname,filesep,'small_',options.usediffmr_coregct]);

else
    ea_reslice_nii([options.subj.preopAnat.(options.subj.AnchorModality).coreg],[options.root,options.patientname,filesep,'small_',options.prefs.prenii_unnormalized],[2 2 2],0);
    MR=ea_load_nii([options.root,options.patientname,filesep,'small_',options.prefs.prenii_unnormalized]);
    delete([options.root,options.patientname,filesep,'small_',options.prefs.prenii_unnormalized]);

end

% CT
ea_reslice_nii([options.subj.postopAnat.(options.subj.postopModality).preproc],[options.root,options.patientname,filesep,'small_',options.prefs.rawctnii_unnormalized],[2 2 2],0);
CT=ea_load_nii([options.root,options.patientname,filesep,'small_',options.prefs.rawctnii_unnormalized]);

delete([options.root,options.patientname,filesep,'small_',options.prefs.rawctnii_unnormalized]);

    disp('Done. Smoothing...');
    MR.img(MR.img<100)=0; % remove noise around brain.
    MR.img=smooth3(MR.img,'gaussian',[11 11 11]);
    CT.img(CT.img<0)=0; % remove negative hounsfield parts.
    CT.img=smooth3(CT.img,'gaussian',[11 11 11]);


vizz=1; % visualization on
alphas=options.coregct.coregthreshs;
%% estimate




if vizz
    ctmr=figure('color','w','name',[options.patientname,': Coregistering CT to MR...'],'NumberTitle','off','Toolbar','none','MenuBar','none','DockControls','off','WindowButtonMotionFcn',@ea_mouseMove);
    axis equal
    axis off
end

disp('Done.');

for alpha=1:length(alphas)

    disp(['*** Pass ',num2str(alpha),'/',num2str(length(alphas)),', alpha: ',num2str(alphas(alpha)),'.']);
    eMR=ea_detect_edges_3d(MR.img,alphas(alpha));
    eCT=ea_detect_edges_3d(CT.img,alphas(alpha));

    eMR(eMR<8)=0;
    eCT(eCT<2)=0;
    eMR(isnan(eMR))=0;
    eCT(isnan(eCT))=0;

    [xx,yy,zz]=ind2sub(size(eMR),find(eMR(:)));
    ptMR=[xx,yy,zz,ones(length(xx),1)]';
    ptMR=MR.mat*ptMR; % mm notation
    [xx,yy,zz]=ind2sub(size(eCT),find(eCT(:)));
    ptCT=[xx,yy,zz,ones(length(xx),1)]';
    ptCT=CT.mat*ptCT; % mm notation
    disp('Done.');

    % define initialization parameters
    if ~exist('M','var') % first grain run..
        % define M based on centroids
        M=eye(4);
        for dim=1:3
            M(dim,4)=mean(ptMR(dim,:))-mean(ptCT(dim,:));
        end
    end
    ptrCT=M*ptCT;
    [~,DCT]=knnsearch(ptMR',ptrCT');
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


        [~,DMR]=knnsearch(ptrCT(1:3,:)',ptMR(1:3,:)');

        [~,DCT]=knnsearch(ptMR(1:3,:)',ptrCT(1:3,:)');

        D=mean([DCT;DMR]);



        % check improvement.
        if D<priorD % improvement
            fprintf('\n %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f \n\n',[D;M(:)]');
            if vizz
                set(0,'CurrentFigure',ctmr);
                plot3(ptMR(1,:),ptMR(2,:),ptMR(3,:),'.','Color',[0,188/255,226/255]);
                hold on
                plot3(ptrCT(1,:),ptrCT(2,:),ptrCT(3,:),'.','Color',[247/255,133/255,20/255]);
                hold off
                axis off
                axis equal
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

close(ctmr);

% M has been estimated and maps from voxels in CT to voxels in MR.
%% export coregistered CT.

matlabbatch{1}.spm.util.reorient.srcfiles = {[options.subj.postopAnat.(options.subj.postopModality).preproc,',1']};
matlabbatch{1}.spm.util.reorient.transform.transM = M;
matlabbatch{1}.spm.util.reorient.prefix = 'r';
jobs{1}=matlabbatch;
try
spm_jobman('run',jobs);
catch
    warning('Pre-coregistration did not work. Please choose a different threshold.');
end
clear jobs matlabbatch


%% use spm coreg to finalize result.

% we might want to try normalizing the rCT using Albada 2009?
% disp('Done. Normalizing CT intensities...');
% rCT=ea_load_nii(,)
% CT.img=reshape(ea_normal(CT.img(:)),size(CT.img,1),size(CT.img,2),size(CT.img,3));

costfuns={'nmi','mi','ecc'};
for costfun=1:3
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[options.subj.preopAnat.(options.subj.AnchorModality).coreg]};
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {[options.root,options.patientname,filesep,'r',options.prefs.rawctnii_unnormalized,',1']};
    matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = costfuns{costfun};
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

    jobs{1}=matlabbatch;
    spm_jobman('run',jobs);
    clear matlabbatch jobs;
end

matlabbatch{1}.spm.util.checkreg.data = {[options.subj.preopAnat.(options.subj.AnchorModality).coreg];
    [options.root,options.patientname,filesep,'r',options.prefs.rawctnii_unnormalized,',1']};
jobs{1}=matlabbatch;
spm_jobman('run',jobs);
clear matlabbatch jobs;

% keep users naming scheme:
try
movefile([options.root,options.patientname,filesep,'r',options.prefs.rawctnii_unnormalized],[options.subj.postopAnat.(options.subj.postopModality).coreg]);
end


function ea_mouseMove(object, eventdata)
[az,el]=view;
view(az+1,el);
