function ea_pca_templates
vizz=0;
usemask=0;

% load data
t1=ea_load_nii([ea_space,'t1.nii']);
t2=ea_load_nii([ea_space,'t2.nii']);
pd=ea_load_nii([ea_space,'pd.nii']);

    maskidx=[1:length(t1.img(:))]';

% feature scale data
t1.img(maskidx)=zscore(t1.img(maskidx));
t2.img(maskidx)=zscore(t2.img(maskidx));
pd.img(maskidx)=zscore(pd.img(maskidx));

% generate joint pointlist
idx=[t1.img(maskidx),t2.img(maskidx),pd.img(maskidx)];
if vizz
    figure
    plot3(idx(1:1000:end,1),idx(1:1000:end,2),idx(1:1000:end,3),'.');
end

[coeff,~,~,~,exp]=pca(idx');
disp([num2str(exp),' % of variance explained by first PC.']);
if vizz
    hold on
    for v=1:3
        plot3([0,V(1,v)],[0,V(2,v)],[0,V(3,v)],'r-');
    end
end

% % get largest eigenvector:
% primeval=S(logical(eye(3)));
% [~,primeval]=max(primeval);
% primevector=V(:,primeval);

if vizz
    figure
    plot3(U(1:1000:end,1),U(1:1000:end,2),U(1:1000:end,3),'.');
end



joint=t1;
joint.img(:)=0;
joint.img(maskidx)=coeff(:,1);
joint.fname=[ea_space,'pca.nii'];
joint.dt=[16,0];
ea_write_nii(joint);

