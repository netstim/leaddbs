function ea_compare_Rmaps(rmap1,rmap2,n1,n2,alpha)
if ~exist('alpha','var')
    alpha=0.05;
end
imr1=ea_load_nii(rmap1);
imr2=ea_load_nii(rmap2);

t_r1 = 0.5*log((1+imr1.img)./(1-imr1.img));
t_r2 = 0.5*log((1+imr2.img)./(1-imr2.img));
z = (t_r1-t_r2)./sqrt(1./(n1-3)+1./(n2-3));
p = (1-normcdf(abs(z),0,1))*2;


imr1.img(p>alpha)=nan;
imr2.img(p>alpha)=nan;

[pth1,r1fname]=fileparts(rmap1);
[pth2,r2fname]=fileparts(rmap2);

R1largerR2=imr1;
R1largerR2.img(imr1.img<imr2.img)=nan;
R1largerR2.img(isnan(R1largerR2.img))=0;

R1smallerR2=imr1;
R1smallerR2.img(imr1.img>imr2.img)=nan;
R1smallerR2.img(isnan(R1smallerR2.img))=0;

R1sigR2=R1largerR2;
R1sigR2.img=R1sigR2.img+R1smallerR2.img;
R1sigR2.fname=fullfile(pth1,[r1fname,'_sig_',r2fname,'.nii']);



R2largerR1=imr2;
R2largerR1.img(imr2.img<imr1.img)=nan;
R2largerR1.img(isnan(R2largerR1.img))=0;

R2smallerR1=imr2;
R2smallerR1.img(imr2.img>imr1.img)=nan;
R2smallerR1.img(isnan(R2smallerR1.img))=0;

R2sigR1=R2largerR1;
R2sigR1.img=R2sigR1.img+R2smallerR1.img;
R2sigR1.fname=fullfile(pth2,[r2fname,'_sig_',r1fname,'.nii']);

ea_write_nii(R1sigR2);
ea_write_nii(R2sigR1);


% This function compare if two correlation coefficients are significantly
% different.
% The correlation coefficients were tansfered to z scores using fisher's r
% to z transformation. 
% ref: http://core.ecu.edu/psyc/wuenschk/docs30/CompareCorrCoeff.pdf
%--------------------------------------------------------------------------
% Inputs: (1) r1: correlation coefficient of the first correlation (2) r2:
% correlation coefficient of the second correlation (3) n1: number of
% samples used to compute the first correlation (4) n2: number of samples
% used to compute the second correlation
%--------------------------------------------------------------------------
% Output: (1) p: p value, the probability that H0 (the correlation
% coefficiets are not different) is correct
%--------------------------------------------------------------------------
% Example :
% x = rand(20,1); 
% y1= x+rand(20,1)*0.05;
% y2= x+rand(20,1)*0.5;
% r1=corr(x,y1);
% r1=corr(x,y2);
% p = compare_correlation_coefficients(r1,r2,length(x),length(x));
%--------------------------------------------------------------------------