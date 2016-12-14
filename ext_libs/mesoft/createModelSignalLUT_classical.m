function Am = createModelSignalLUT_classical(ten,alpha,tissue_parameters,sinterp,Pstruc,fid)

alpha = 0;
dirs = sinterp.bDir;

Da = tissue_parameters(1);
Dp = tissue_parameters(2);
vf = tissue_parameters(3);


[~, W] = createWeightingScheme(ten,Pstruc.b_weighting);

E = squeeze(sum(cmult(dirs,ten,1,1).*repmat(dirs',[1 1 size(ten,3)]),2));
Erad = repmat(permute(ten(1,1,:)+ten(2,2,:)+ten(3,3,:),[1 3 2]),[size(dirs,2) 1]);
Am = vf*exp(-(E)*Da) + (1-vf)*exp(-(E)*Da - (Erad-E)*Dp);
Am = Am - alpha *repmat(sum(cmult(Am,W,2,1),2),[1 size(Am,2)])/trace(W);
Am = cmult(Am,W,2,1)*tissue_parameters(5);

%%
%Am = (exp(-15*E)-0.8*exp(-10*E))*7;
%  figure(9898); 
%  clf; 
%  di = [bten(:,bv>1.5) -bten(:,bv>1.5)];
%  r = 6;
%  q = [E(r,bv>1.5) E(r,bv>1.5)]; 
%  q = exp(-15*q)-0.8*exp(-q*10);
%  showGlyph(di,q);
 
 %plot(abs(bten(:,bv>1.5)'*dirs(:,r)),q(1:end/2),'*')
 