

function hr = nifti_to_HARDI_phasecorr(nifti_mag,nifti_phase);

op = cd(nifti_mag);
x = dir('*.nii'); hr = nifti_to_hardi({x.name});
cd(op);
cd(nifti_phase);
x = dir('*.nii'); hr_ph = nifti_to_mrstruct('volume',{x.name});
cd(op);

b0idx = find(hr.user.bfactor==0);

hr.dataAy = hr.dataAy.*exp(i*pi*double(hr_ph.dataAy)/4096);

im = hr.dataAy;


%%
[X Y] = ndgrid(-12:12);
R2 = X.^2 + Y.^2;
sig = 2;
gauss = exp(-R2/(2*sig^2)); %.*R2.^2;
gauss = gauss/sum(gauss(:));
gauss(13,13) = 0;
b0 = mean(im(:,:,:,b0idx),4);
%b0 = b0*0+1;

ims = imfilter(im.*repmat(conj(b0),[1 1 1 size(im,4)]),gauss,'circular');
ims = im .* conj(ims) ./ (eps+abs(ims)) ;
ims = ims .*repmat((conj(b0)./(eps+abs(b0))),[1 1 1 size(im,4)]);


hr.dataAy = single(real(ims));


mask = not(imfilter(im(:,:,:,1),fspecial('gaussian',[15 15],5))>10);
%mask = logical(1+0*mask);



b =80;
%q = ims(2:20,2:20,:,b);
%q1 = im(2:20,2:20,:,b);
%q = ims(2:end-1,2:end-1,2:end-1,b);
%q1 = im(2:end-1,2:end-1,2:end-1,b);
q = ims(:,:,:,b);
q1 = im(:,:,:,b);
b = -60:2:60; 
figure(100);
subplot(2,2,1);
h = log(1+hist3([real(q(mask(:))) imag(q(mask(:)))],{b b}));
h(b==0,b==0) =0;
imagesc(b,b,h);
subplot(2,2,2);
imagesc(gauss);
subplot(2,2,3);
h = log(1+hist3([real(q1(mask(:))) imag(q1(mask(:)))],{b b}));
h(b==0,b==0) =0;
imagesc(b,b,h);


figure(200);
r = [0 300];
s = 5;
subplot(2,2,1);
imagesc(real(q(:,:,s)),r);
subplot(2,2,2);
imagesc(imag(q(:,:,s)),r);
subplot(2,2,3);
imagesc(abs(q(:,:,s)),r);
subplot(2,2,4);
imagesc(angle(q1(:,:,s)),[-pi pi]);

