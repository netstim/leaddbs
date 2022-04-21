function ea_pca_templates
vizz=0;
usemask=1;




% load data
t1=ea_load_nii([ea_space,'t1.nii']);
t2=ea_load_nii([ea_space,'t2.nii']);
pd=ea_load_nii([ea_space,'pd.nii']);

if usemask
    msk=ea_load_nii([ea_space,'mask.nii']);
    maskidx=find(msk.img(:)>0);
    fullidx=[1:length(t1.img(:))]';
else
    maskidx=[1:length(t1.img(:))]';
    fullidx=[1:length(t1.img(:))]';
end
% feature scale data
t1.img(fullidx)=zscore(t1.img(fullidx));
t2.img(fullidx)=zscore(t2.img(fullidx));
pd.img(fullidx)=zscore(pd.img(fullidx));

% generate joint pointlist
idx=[t2.img(maskidx),t1.img(maskidx),pd.img(maskidx)];
allidx=[t2.img(fullidx),t1.img(fullidx),pd.img(fullidx)];
if vizz
    figure
    plot3(idx(1:1000:end,1),idx(1:1000:end,2),idx(1:1000:end,3),'.');
end

% train PCA only on brain part of the image:
[wcoeff,~,~,~,exp]=pca(idx);
% orthogonalize PC coeffs:
coefforth = inv(diag(std(idx)))*wcoeff;
% apply orthogonalized coefficients to full images:
cscores = zscore(allidx)*coefforth;
disp([num2str(exp(1)),' % of variance explained by first PC.']);


if vizz
    hold on
    for v=1:3
        plot3([0,V(1,v)],[0,V(2,v)],[0,V(3,v)],'r-');
    end
end

joint=t1;
joint.img(:)=0;
joint.img(fullidx)=cscores(:,1);

joint.img=joint.img-min(joint.img(:));
joint.img=joint.img/max(joint.img(:));
joint.img=joint.img*(100); % 
joint.fname=[ea_space,'pca.nii'];
joint.dt(1) = 4;
ea_write_nii(joint);

%% uncomment to also write out 2nd and 3rd PCs
% joint=t1;
% joint.img(:)=0;
% joint.img(fullidx)=cscores(:,2);
% 
% joint.img=joint.img-min(joint.img(:));
% joint.img=joint.img/max(joint.img(:));
% joint.img=joint.img*(100); % 
% joint.fname=[ea_space,'pca2.nii'];
% joint.dt(1) = 4;
% ea_write_nii(joint);
% 
% joint=t1;
% joint.img(:)=0;
% joint.img(fullidx)=cscores(:,3);
% 
% joint.img=joint.img-min(joint.img(:));
% joint.img=joint.img/max(joint.img(:));
% joint.img=joint.img*(100); % 
% joint.fname=[ea_space,'pca3.nii'];
% joint.dt(1) = 4;
% ea_write_nii(joint);

