function ea_3dgauss_from_pointlist(ptlist,fname)  

[gauss3d,mat]=d3gauss(ptlist(:,1),ptlist(:,2),ptlist(:,3));

V=ea_synth_nii(ea_niigz(fname),mat,[16,1],gauss3d);
gauss3d=gauss3d/max(gauss3d(:));
gauss3d(gauss3d<0.001)=0;
spm_write_vol(V,gauss3d);



function [P_Gaussian,mat]=d3gauss(var1,var2,var3)


Data=[var1,var2,var3];
Data(isnan(Data(:,1)),:)=[];
Data(isnan(Data(:,2)),:)=[];
Data(isnan(Data(:,3)),:)=[];
offset=10;
Mu = mean(Data);
X = linspace(Mu(1)-offset,Mu(1)+offset,250);
Y = linspace(Mu(2)-offset,Mu(2)+offset,250);
Z = linspace(Mu(3)-offset,Mu(3)+offset,250);

D = length(Data(1,:));

Sigma = cov(Data);

P_Gaussian = zeros(length(X),length(Y),length(Z));
cnt=1;
mapfrom=zeros(10000,4);
mapto=zeros(10000,4);
incnt=1;
for i=1:length(X)
    for j=1:length(Y)
        for k=1:length(Z)
            x = [X(i),Y(j),Z(k)];
            P_Gaussian(i,j,k) = 1/((2*pi)^(D/2)*sqrt(det(Sigma)))...
                *exp(-1/2*(x-Mu)*Sigma^-1*(x-Mu)');
            
            incnt=incnt+1;
            if (incnt/100)==round(incnt/100); % sample every hundred points
                if cnt<10001
                    mapfrom(cnt,:)=[i,j,k,1];
                    mapto(cnt,:)=[x,1];
                    cnt=cnt+1;
                end
            end
        end
    end
end

%mapfrom(:,1:3)=mapfrom(:,1:3)+offset;

mat=linsolve(mapfrom,mapto)';




