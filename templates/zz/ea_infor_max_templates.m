vizz=0;
usemask=1;

% load data
t1=ea_load_nii('mni_hires_t1.nii');
t2=ea_load_nii('mni_hires_t2.nii');
pd=ea_load_nii('mni_hires_pd.nii');
mask=ea_load_nii('bmni_hires.nii.gz');

if usemask
    maskidx=mask.img(:);
    maskidx=maskidx>0;
else
    maskidx=[1:length(mask.img(:))]';
end

% ea_normal data
t1.img(maskidx)=ea_normal(t1.img(maskidx));
t2.img(maskidx)=ea_normal(t2.img(maskidx));
pd.img(maskidx)=ea_normal(pd.img(maskidx));

% generate joint pointlist
idx=[t1.img(maskidx),t2.img(maskidx),pd.img(maskidx)];
if vizz
    figure
    plot3(idx(1:1000:end,1),idx(1:1000:end,2),idx(1:1000:end,3),'.');
end

[U,S,V]=svd(idx,0);

if vizz
    hold on
    for v=1:3
        plot3([0,V(1,v)],[0,V(2,v)],[0,V(3,v)],'r-');
    end
end

% get largest eigenvector:
primeval=S(logical(eye(3)));
[~,primeval]=max(primeval);
primevector=V(:,primeval);

if vizz
    figure
    plot3(U(1:1000:end,1),U(1:1000:end,2),U(1:1000:end,3),'.');
end



joint=t1;
joint.img(maskidx)=U(:,1);
joint.img(maskidx)=ea_normal(joint.img(maskidx));
% fill non-brain regions with mean information.
if usemask
joint.img(~maskidx)=ea_normal(mean([t1.img(~maskidx),t2.img(~maskidx),pd.img(~maskidx)],2));
end
%joint.img(~maskidx)=0;
joint.fname='mni_hires_svd.nii';
ea_write_nii(joint);
keyboard
