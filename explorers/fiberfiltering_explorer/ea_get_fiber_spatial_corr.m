function spatio_corr_mat = fiber_spatial_rev3(fibers,trgt_coor,sig,vcnty_thr)

% Min Jae Kim (mkim@bwh.harvard.edu)

%% Organizing Spatial Coordinates of Each Fibers 
num_fibers=length(unique(fibers(:,4))); % number of fibers
fiber_idxi=cell(num_fibers,2); % data structure
for i=1:num_fibers
    fiber_idxi{i,1}=i; %1st column: fiber number 
    idxi=find(fibers(:,4)==i);
    coors=fibers(idxi,:);
    fiber_idxi{i,2}=coors(:,1:3); %2nd column: fiber coordinates
end

%% Step 0: Initializing Parameters 
method='Spearman'; % Method for Correlation: (1) Spearman or (2) Pearson
MNI_temp=ea_load_nii([ea_space,'t1.nii']);
dim=MNI_temp.dim;
%% Method 1: Spatial Correlation Method
spatio_corr_mat=zeros(num_fibers,num_fibers);
%spatio_corr_mat=zeros(20,20);% sample matrix for first 20 fibers
%% Processing Fibers (3D Gaussianization)
fib_gaus_cell=cell(num_fibers,1);
idxi_gaus_cell=cell(num_fibers,1);
%tic
parfor i=1:num_fibers
  smp1=fiber_idxi{i,2};
  [smp1_gauss,idxi_1]=ea_apply_fiber_gaussmooth(smp1,sig,trgt_coor,vcnty_thr,dim);
  fib_gaus_cell{i,1}=smp1_gauss;
  idxi_gaus_cell{i,1}=idxi_1;
  disp("Processing Fiber Number #"+string(i))
end
%toc

%% Performing Spatial Correlation (3D Gaussianization)
for i=1:num_fibers
  smp1=fib_gaus_cell{i,1};
  idxi_1=idxi_gaus_cell{i,1};
  tic
  for j=1:num_fibers
    smp2=fib_gaus_cell{j,1};
    idxi_2=idxi_gaus_cell{j,1};
    [C,ia,ib] = union(idxi_1,idxi_2);
    [C_1,i_1,~]=intersect(C,idxi_1);
    [C_2,i_2,~]=intersect(C,idxi_2);
    temp1=zeros(length(C),1);
    temp1(i_1)=smp1;
    temp2=zeros(length(C),1);
    temp2(i_2)=smp2;
    rho=corr(temp1,temp2,'Type',method,'Rows','complete');
    spatio_corr_mat(i,j)=rho;
  end
  toc
end
%ea_spatial_corr
% Visualization
figure
heatmap(spatio_corr_mat);
grid off
colormap jet
xlabel("Fiber Number")
ylabel("Fiber Number");
title("Pairwise Spatial Correlation (r) Between Fibers");
end
