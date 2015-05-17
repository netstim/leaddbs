function createModelModelLUT_classical(ten,alpha,tissue_params,Pstruc,fid)


alpha =0;
dirs = load('dwidirections.mat');
dirs = [dirs.dirs128 -dirs.dirs128];

[~, W] = createWeightingScheme(ten,Pstruc.b_weighting);

Da = tissue_params(1);
Dp = tissue_params(2);
vf = tissue_params(3);


P = squeeze(ten(1,1,:))*dirs(1,:).^2 + squeeze(ten(2,2,:))*dirs(2,:).^2 + squeeze(ten(3,3,:))*dirs(3,:).^2 + ...
   + 2*(squeeze(ten(1,3,:))*(dirs(1,:).*dirs(3,:))+ squeeze(ten(1,2,:))*(dirs(1,:).*dirs(2,:)) +squeeze(ten(3,2,:))*(dirs(3,:).*dirs(2,:)));
Q = repmat(squeeze(ten(1,1,:)) + squeeze(ten(2,2,:)) + squeeze(ten(3,3,:)),[1 size(dirs,2)]);

M = exp(-Da*P)*vf + exp(-Da*P - Dp*(Q-P))*(1-vf);
M = M - alpha*repmat(sum(W*M),[size(M,1) 1]) / trace(W);


MM = M'*W*M;


C = abs(dirs'*dirs);
C(C>1) = 1;
[sd idx] = sort(C(:));
beta = MM(idx);

T = [sd.^0 sd.^2 sd.^4 sd.^6];
coeff = pinv(T)*beta*tissue_params(5)^2;

x = (0:0.01:1)';
approx = [x.^0  x.^2 x.^4 x.^6]*coeff;


%figure(11); plot(acos(sd),beta,acos(x),[x.^0  x.^2 x.^4 x.^6]*coeff);

createCcode(fid,coeff);





function createCcode(fid,coeff)



fprintf(fid,'REAL EnergyComputerBase::modelmodelCorrFixTis(REAL dot)\n {\n');

fprintf(fid,'\tREAL x2 = dot*dot; \n');
fprintf(fid,'\treturn (%.8f + x2*%.8f + x2*x2*%.8f + x2*x2*x2*%.8f)*prior_strength; \n }\n',coeff(1),coeff(2),coeff(3),coeff(4));
return;









function h = labelSymMonoms(idx)

P = [2 1 3 ; 3 2 1 ; 1 3 2 ; 3 1 2 ; 2 3 1];

h = zeros(size(idx,1),1);
cnt = 1;
for j = 1:size(idx,1),
    if h(j) == 0,
        e = zeros(size(idx,1),1);
        e(j) = 1;
        for k = 1:5,
            Pidx =  reshape(P(k,idx),size(idx));
            e = e + all(Pidx==repmat(idx(j,:),[size(idx,1) 1]),2);
        end;
        h(e>0) = cnt;
        cnt = cnt + 1;
    end;
end;    
    
    





