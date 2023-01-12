
function ea_mask_electrode(CT_path,reco_file)
lead path
load(char(reco_file));
CT=ea_load_nii(char(CT_path)); %should be normalized CT image (i.e glpostop_ct.nii') 
CT_img=CT.img;
disp(['Step 1: Running PACER Algorithm...']);
rcnt=PaCER(char(CT_path));
side={"right","left"};
%% 
for hh=1:2
full_trj=rcnt{1,hh}.skeleton;
scrf=reco.scrf.trajectory{1,hh};
ori=rcnt{1,hh}.skeleton;
if ori(1,3)>scrf(1,3)
    coor1=scrf(1,:);
    coor2=ori(1,:);
    dt2=ori(2,:)-ori(1,:);
    
    %int_1=[coor1(1):dt2(1):coor2(1)]';
    %int_2=[coor1(2):dt2(2):coor2(2)]';
    int_3=[coor1(3):dt2(3):coor2(3)]';
    
    dt1=(coor2(1)-coor1(1))/length(int_3);
    int_1=[coor1(1):dt1:coor2(1)]';
    int_1(end)=[];
    
    dt2=(coor2(2)-coor1(2))/length(int_3);
    int_2=[coor1(2):dt2:coor2(2)]';
    int_2(end)=[];
    
    int=[int_1,int_2,int_3];
    full_trj=[int;full_trj];
end 

XYZmm=[full_trj';ones(1,length(full_trj))];
pcavox = zeros(CT.dim);
transf=inv(CT.mat);
coors_col=[];
for i = 1:size(XYZmm,2)
    coors=transf*XYZmm(:,i);
    coors=floor(coors(1:3));
    coors_col=[coors_col;coors'];
end 
% 



%
temp= zeros(CT.dim);
z_cols=coors_col(:,3); %Voxel-Space in CT of the trajectory
thr=0.5;
for jj=1:length(z_cols)
img_slc=(squeeze(CT_img(:,:,z_cols(jj))));
thr_int=min(img_slc(:))+((max(img_slc(:))-min(img_slc(:)))*thr);
idxi=find(img_slc<thr_int);
idxi_2=find(img_slc>=thr_int);

img_slc_thr=img_slc;
img_slc_thr(idxi)=0;
img_slc_thr(idxi_2)=1;
%imshow(img_slc_thr)
%
[X,Y]=ind2sub(size(img_slc_thr),find(img_slc_thr==1));
candi=[X,Y];
gntrd=coors_col(jj,:);
gntrd=[gntrd(1),gntrd(2)];

    euc_col=[];
    for ii=1:length(candi)
        x_dist=abs(gntrd(1)-candi(ii,1)).*CT.voxsize(1);
        y_dist=abs(gntrd(2)-candi(ii,2)).*CT.voxsize(2);
        euc=sqrt(x_dist.^2+y_dist.^2);
        if euc<=2
            euc_col=[euc_col;ii];
        end 
    end 
    
    for gg=1:length(euc_col)
        idxi=euc_col(gg,:);
        coor=[X(idxi),Y(idxi),z_cols(jj)];
        temp(coor(1),coor(2),coor(3))=1;
    end    
end 
%% Writing Nifiti Files

CT_rev=CT;
CT_rev.img=temp;
out_tit=side{hh}+"_electrode_mask.nii";
CT_rev.fname=char(out_tit);
ea_write_nii(CT_rev);

end 
