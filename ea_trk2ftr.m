function [fibers,idx] = ea_trk2ftr(trk_in)
% transforms .trk to fibers.
% CAVE: fiber information in the trk needs to conform to the MNI space!

sd=ea_getspacedef;
%% basis information
template_in = ea_niigz([ea_space,sd.templates{1}]);

%% read .trk file
[header,tracks]=ea_trk_read(trk_in);

%% create variables needed
ea_fibformat='1.0';
fourindex=1;
idx=[tracks(:).nPoints]';
idx2 =cumsum(idx);
fibers = NaN(sum(idx),4);

start = 1;
fibercount = numel(idx);
ea_dispercent(0,'Converting trk to fibers.');
for a=1:fibercount
    ea_dispercent(a/fibercount);
    stop = idx2(a);
    fibers(start:stop,1:3)=tracks(a).matrix;
    fibers(start:stop,4)=a;
    start = stop+1;
end
ea_dispercent(1,'end');
fibers=double(fibers);

fibers = fibers';

%% transform fibers to template origin
%nii=ea_load_nii('/PA/Neuro/_projects/lead/lead_dbs/templates/space/MNI_ICBM_2009b_NLIN_ASYM/atlases/Macroscale Human Connectome Atlas (Yeh 2018)/FMRIB58_FA_1mm.nii.gz')
answ=questdlg('Please select type of .trk data', 'Choose origin of .trk data', 'DSI Studio / QSDR', 'Normative Connectome / 2009b space', 'Select .nii file', 'DSI Studio / QSDR');

switch answ
    
    case 'DSI Studio / QSDR'
        % http://dsi-studio.labsolver.org/Manual/Reconstruction#TOC-Q-Space-Diffeomorphic-Reconstruction-QSDR-
        fibers(1,:)=78.0-fibers(1,:);
        fibers(2,:)=76.0-fibers(2,:);
        fibers(3,:)=-50.0+fibers(3,:);    
        
    otherwise
        
        switch answ
            case 'Normative Connectome / 2009b space'
                
                nii=ea_load_nii(template_in);
                
            case 'Select .nii file'
                [fname,pathname]=uigetfile({'.nii','.nii.gz'},'Choose Nifti file for space definition');
                nii=ea_load_nii(fullfile(pathname,fname));
        end
        % nii.mat(1,1)=1;
        % nii.mat(2,2)=1;
        % nii.mat(3,3)=1;
        tmat=header.vox_to_ras;
        tmat(1:3,4)=header.dim';
        tfib=[fibers(1:3,:);ones(1,size(fibers,2))];
        tfib=tmat*tfib;
        
        %tfib=nii.mat*tfib;
        fibers(1:3,:)=tfib(1:3,:);
end


%fibers = fibers + repmat([nii.mat(1,4) nii.mat(2,4) nii.mat(3,4) 0]',1,size(fibers,2));
