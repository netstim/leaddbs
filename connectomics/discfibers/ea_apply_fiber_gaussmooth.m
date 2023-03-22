function [smp1_gauss,idxi_1]=fib_gauss(smp1,sig,trgt_coor,vcnty_thr,dim)

% Min Jae Kim (mkim@bwh.harvard.edu)

%% Step 1: Imposing Vicinity Threshold

% Running Fibertract 1
rmv_1=[];
for hh=1:length(smp1) 
    vc_euc=sqrt(sum((smp1(hh,:) - trgt_coor) .^ 2));
    if vc_euc>vcnty_thr
        rmv_1=[rmv_1;hh];
    end
end 
smp1(rmv_1,:)=[];

%% MNI to Voxel
% Running Fiber Tract 1 
smp1_vox=round(ea_mm2vox(smp1(:,1:3),[ea_space,'t1.nii']));
smp1_MNI=zeros(dim);
for i=1:length(smp1_vox)
    smp1_MNI(smp1_vox(i,1),smp1_vox(i,2),smp1_vox(i,3))=1;
end 
smp1_gauss=imgaussfilt3(smp1_MNI,sig); %3D Gaussian
smp1_gauss=smp1_gauss(:);
idxi_1=find(smp1_gauss~=0);
smp1_gauss=smp1_gauss(idxi_1);
end 
