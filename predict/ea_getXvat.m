function [X,XR,nii]=ea_getXvat(M,options)
XR=nan;

selectedregressor=M.clinical.vars{M.ui.clinicallist};
if size(selectedregressor,2)==1
    bihemispheric=0;
elseif size(selectedregressor,2)==2
    bihemispheric=1;
else
    ea_error('Please select a regressor with entries for each hemisphere or each patient to perform this action.');
end

cnt=1;

for pt=1:length(M.patient.list)
    %for left side
    fname_l=[options.root,options.patientname,filesep,'statvat_results',filesep,'s',num2str(pt),'_lh.nii'];
    if exist(fname_l,'file')>0
        nii=ea_load_nii(fname_l);
        %init outputs X and XR if necessary
        if ~exist('X','var')
            X=nan(length(M.patient.list),numel(nii.img));
            if bihemispheric
                XR=X;
            end
        end
        X(cnt,:)=nii.img(:);
    end
    
    %for right side
    if bihemispheric
        fname_r=[options.root,options.patientname,filesep,'statvat_results',filesep,'s',num2str(pt),'_rh.nii'];
        if exist(fname_r,'file')>0
            nii=ea_load_nii(fname_r);
            %init outputs X and XR if necessary
            if ~exist('X','var')
                X=nan(length(M.patient.list),numel(nii.img));
                if bihemispheric
                    XR=X;
                end
            end
            XR(cnt,:)=nii.img(:);
            XR(cnt,:)=logical(XR(cnt,:));
        end
    else
        fname_rflip=[options.root,options.patientname,filesep,'statvat_results',filesep,'s',num2str(pt),'_rh_flipped.nii'];
        if exist(fname_rflip,'file')>0
            nii=ea_load_nii(fname_rflip);
            %init outputs X and XR if necessary
            if ~exist('X','var')
                X=nan(length(M.patient.list),numel(nii.img));
                if bihemispheric
                    XR=X;
                end
            end
            X(cnt,:)=X(cnt,:)+nii.img(:)';
        end
    end
    
    X(cnt,:)=logical(X(cnt,:));    
    cnt=cnt+1;
end
nii.dt(1) = 8;
