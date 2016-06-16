
function maps = computeParameterMaps(ftr,M);

      

  
      idxmask = ftr.user.S2(:,:,:,2)>0;
      
      % kernel for gathering values withtin original DWI-voxel
      osamp =  (ftr.trackParam.params.p_wid); 
      osamp = osamp/M;
      st = max(ceil(osamp/2),1);
      ker = ones(osamp,osamp,osamp); %ker = ker /sum(ker(:));
      
      % compute averages
      %  vf(:,:,:,1) - raw vfi 
      %  vf(:,:,:,2) - projection of raw model onto constant
      %  vf(:,:,:,3) - segment's contributions within voxel  \sum_i w_i
      
      
%%      
      vf = ftr.user.vf;
      mask = double(not(isnan(vf(:,:,:,1))));
      vf(isnan(vf)) = 0;


      if 0,
          sigma = 1*osamp;
          ng = ceil(sigma*3); ng = ng + mod(ng,2) + 1;
          gau = fspecial('gaussian',[ng 1],sigma); [gx gy gz] = ndgrid(gau); gau = gx.*gy.*gz;
          sz = size(mask);
          gau = padarray(gau,floor((sz-ng)/2),'post');
          gau = padarray(gau,floor((sz-ng)/2)-mod(sz,2) +1,'pre');
          gau = fftn(ifftshift(gau));

          masksm = (real(ifftn(fftn(mask).*gau)));
          err = (real(ifftn(fftn(ftr.user.stderr).*gau)));
          vfsm = vf;
          for k = 1:size(vfsm,4),
              for j = 1:size(vfsm,5),
                vfsm(:,:,:,k,j) = real(ifftn(fftn(vf(:,:,:,k,j)).*gau));   
              end;
          end
      else
      
          gau = fspecial('gaussian',[5 1],1*osamp); [gx gy gz] = ndgrid(gau); gau = gx.*gy.*gz;
          masksm = imfilter(double(mask),gau);
          vfsm = imfilter(vf,gau);
      end;
      
      vfsm = vfsm ./ repmat(eps+masksm,[1 1 1 size(vfsm,4) size(vfsm,5)]);
      
      
      vfsm = vfsm(st:osamp:end,st:osamp:end,st:osamp:end,:,:);
      
      mask = mask(st:osamp:end,st:osamp:end,st:osamp:end);
      vfsm = vfsm(:,:,:,:,2:end) .* repmat(1./(eps+vfsm(:,:,:,:,1)),[1 1 1 1 size(vfsm,5)-1]);
      vfsm = vfsm.*repmat(mask,[1 1 1 size(vfsm,4) size(vfsm,5)]);
      

      %stackview(vfsm(:,:,:,1),@(x)imagesc(x,[0 0.6]))
      
  %%    
      lmax = size(vfsm,4); 
      B = ftr.trackParam.bvals;      
      
      
      
      
%%      
      if false,

%%          
%           lmax_cut = min(5,lmax);
% 
%           B = ftr.trackParam.bvals;     
%          %B = B(2:end);
%           LM = 0:2:(lmax_cut*2-1); LM = (2*LM+1);%.*exp(-0.4*LM);
%           
%           LM(2) = LM(2);
%           t = sqrt(0:0.05:1);
%           Mm = myleg(lmax_cut*2-1,t'); Mm = Mm(:,1:2:end); Mm = Mm.*repmat(LM,[size(Mm,1) 1]);
%           S = cmult(vfsm(:,:,:,1:lmax_cut,:),Mm,4,2);
%           
%           
%           S = abs(cat(4,ones(size(S,1),size(S,2),size(S,3),1,size(S,5)),S));
%           S = abs(S);
%           S = log(S);
% 
%           [B T] = ndgrid(B,t.^2); 
%           B=B(:); T=T(:); W = 1+0*sqrt(1-T(:));
%           Ra = B-B.*T;
%           Ax = B.*T;
% 
%           P = [Ax Ra Ax.^2 Ra.^2 Ra.*Ax  Ra*0+1 ];
%           Re = diag([0 0 0 0 0 0]);
%           Pinv = inv(P'*diag(W)*P + Re*1)*P'*diag(W);
%           A = cmult(S(:,:,:,:),Pinv,4,2);
%           Spred2 = cmult(S(:,:,:,:),P*Pinv,4,2);
%           
%           
%           
          
%%          

M = reshape(vfsm,[size(vfsm,1)*size(vfsm,2)*size(vfsm,3)  size(vfsm,4) size(vfsm,5)]);
M = M(mask(:)>0,:,:);

modelstruc = evalin('base','modelstruc');


%SNR = ftr.user.b0avg(mask(:)>0)./(modelstruc.maxnz*ftr.user.b0avg(mask(:)>0) +err(mask(:)>0));
%SNR = ftr.user.b0avg(mask(:)>0)./(modelstruc.maxnz*ftr.user.b0avg(mask(:)>0) +ftr.user.stderr(mask(:)>0));
snr = ftr.user.b0avg / mean(err(mask(:)>0));

SNR = snr(mask(:)>0);

est_tmp = modelstruc.apply(SNR,M);
est = nan([prod(size(mask)) size(est_tmp,2)]);
est(mask(:)>0,:) = est_tmp;
est = reshape(est,[(size(mask)) size(est_tmp,2)]);


% 
% lmax = 4;
% M = reshape(vfsm(:,:,:,1:(lmax/2+1),:),[size(vfsm,1)*size(vfsm,2)*size(vfsm,3) (lmax/2+1) size(vfsm,5)]);
% 
% tmp = squeeze(M(:,2:end,1:end-1)./repmat(sum(M(:,2:end,:),3),[1 1 size(M,3)-1]));
% M = [squeeze(M(:,1,:)) tmp(:,:) ];
% 
% M = M(:,:);
% maxN = size(M,2);
% P = [M(:,1)*0+1 ];
% M = (M);
% for k = 1:maxN,
%     P = cat(2,P,M(:,k));
%     for j = k:maxN,
%         P = cat(2,P,M(:,k).*M(:,j));
%         for r = j:maxN,
%             P = cat(2,P,M(:,k).*M(:,j).*M(:,r));
%         end
%     end;        
% end;
% 
% 
% alpha = evalin('base','alpha');
% est = reshape(P*alpha,size(vfsm,1),size(vfsm,2),size(vfsm,3),4);


  %%        
          
          
          
          
          
          
          
          
          
          
          
          

%           
%           P = [Ax Ra B.^2 B.^2.*T B.^2.*T.^2  Ra*0+1 ];
%           %P = [Ax Ra B.^2 B.^2.*T.^2  Ra*0+1 ];
%           Pinv = inv(P'*diag(W)*P)*P'*diag(W);
%           A = cmult(S(:,:,:,:),Pinv,4,2);
%           Spred2 = cmult(S(:,:,:,:),P*Pinv,4,2);
%           
%           
%           Dax = -A(:,:,:,1);
%           Drad = -A(:,:,:,2);
%           Kax = 2*A(:,:,:,3);
%           Krad = 2*A(:,:,:,4);
%           K3 = 2*A(:,:,:,5);
%           S0 = 0*A(:,:,:,4);
%           FW = 1-exp(S0);
%          S0 = A(:,:,:,6);
%          FW = 1-exp(S0);

          
%           wm = evalin('base','WMmask(:)');
%           
%           
%           nD1 = volume2mosaic(Dax); nD1 = nD1(wm);
%           nD2 = volume2mosaic(Drad); nD2 = nD2(wm);
%           nK1 = volume2mosaic(Kax); nK1 = nK1(wm);
%           nK2 = volume2mosaic(Krad); nK2 = nK2(wm);
%                     
%           sols2 = evalin('base','sols2');
%           
% 
%           s1 = numevalsol(sols2,'D1',Dax(:),'D2',Drad(:),'K1',Kax(:),'K2',Krad(:));         
%           ns = 2;
%           
%          
%           Dax_int = real(reshape(s1.Di(ns,:),size(Dax)));
%           Drad_int = Dax_int*0;
%           Dax_ext = real(reshape(s1.De(ns,:),size(Dax)));
%           Drad_ext = real(reshape(s1.De(ns,:),size(Dax)));
%           volfrac = real(reshape(s1.f(ns,:),size(Dax)));
%           volfrac_csf = FW; %real(reshape(s1.fw(ns,:),size(Dax)));
%           disper = real(reshape(s1.a(ns,:),size(Dax)));
          %volfrac = disper+0.5;
          
%           D1 = Dax; D2 = Drad; K1 = Kax; K2 = Krad; sg= -1;
%           Dax_int =  D1 + sg*D2.*abs(K1./K2).^(1/2) ;%Dax.*Krad - Drad.*(sg*abs(Kax.*Krad).^(1/2) ))./Krad;
%           Drad_int =  Dax_int*0;
%           Dax_ext =  D1 - sg*(K2.*abs(K1./K2).^(1/2))./D2; %(sg*abs(Kax.*Krad).^(1/2) + Dax.*Drad)./Drad;
%           Drad_ext = (Drad.^2 + Krad)./Drad;
%           volfrac   = Krad./(Drad.^2 + Krad);
%           disper = volfrac*0;
          
          
%           
%          vl = [7 14 9 ; 5 20 10; 5 15 12; 5 19 11];
%          % vl = [20 13 2 ; 13 20 2; 12 4 2; 12 20 2];
%           
%           
%           figure(3849);
%           clf;
%           
%           for k = 1:4,
%             subplot(2,2,k);
%             v = vl(k,:);
%             %bar([vfsm(v(1),v(2),v(3),1,1),vfsm(v(1),v(2),v(3),1,2) ... 
%             %    vfsm(v(1),v(2),v(3),1,1)./vfsm(v(1),v(2),v(3),1,1) ]);
%             bar((squeeze([vfsm(v(1),v(2),v(3),:,:)])));
%             %axis([0 4 -10 0]);
%           end;
%           
%           
%           
%           figure(3847);
%           clf;
%           
%           for k = 1:4,
%             subplot(2,2,k);
%             v = vl(k,:);
%           
%           col = 'rgb';
%          s0 = squeeze(S0(v(1),v(2),v(3)));
%          Di = Dax_int(v(1),v(2),v(3));
%          Dip = Drad_int(v(1),v(2),v(3));
%          De = Dax_ext(v(1),v(2),v(3));
%          Dep = Drad_ext(v(1),v(2),v(3));
%          f =  volfrac(v(1),v(2),v(3));
%          fw =  volfrac_csf(v(1),v(2),v(3));
%          a =  disper(v(1),v(2),v(3));
%       %   [Di De Dep f fw a]
%          a = [Dax(v(1),v(2),v(3)) Drad(v(1),v(2),v(3)) Kax(v(1),v(2),v(3)) Krad(v(1),v(2),v(3)) K3(v(1),v(2),v(3))];
%          sn = @(x) log(x);
%          for b = 1:3,
%               idx = B==b;
%               plot(T(idx),squeeze(sn(exp(S(v(1),v(2),v(3),idx)))),[col(b) '*']); hold on;
%               
% %              Spred = f.*exp(-Di*Ax-Dip*Ra) + (1-f).*exp(-De*Ax-Dep*Ra);
% %              Spred = exp(-Di*Ax)*f + (1-f).*exp(-De*B) ;
% 
%   %           Spred = (dispstick(Di,a*2,[1 0 0]',[Ax(:)' ; Ra(:)'])*f + (1-f-fw).*exp(-De*Ax-Dep*Ra)) + fw*exp(-3*(Ax+Ra)) ;
% %              plot(T(idx),sn(Spred(idx)),[col(b) '-']); hold on;
%               
%               Sp2 = squeeze(Spred2(v(1),v(2),v(3),idx));
%               
%               plot(T(idx),sn(exp(Sp2)),[col(b) '--']); hold on;
%               axis([0 1 -3 0]);
%               %axis([0 1 0 1]);
%               title(a');
% 
%           end;
%           end;
          %%
          %volfrac = disper+0.5;
          
%           
%           
%           assignin('base','nD1',nD1);
%           assignin('base','nD2',nD2);
%           assignin('base','nK1',nK1);
%           assignin('base','nK2',nK2);
%           
%           % (D1*D2 - (D1*K2/D2 + (K1*K2)^(1/2)) + D1*K2/D2)/D2
% 
%           sg = -1;
%           
%           D1 = Dax;
%           D2 = Drad;
%           K1 = Kax ;
%           K2 = Krad ;
          
%           Dax_int =  D1 + D2.*abs(K1./K2).^(1/2) ;%Dax.*Krad - Drad.*(sg*abs(Kax.*Krad).^(1/2) ))./Krad;
%           Dax_ext =  D1 - (K2.*abs(K1./K2).^(1/2))./D2; %(sg*abs(Kax.*Krad).^(1/2) + Dax.*Drad)./Drad;
%           Drad_ext = (Drad.^2 + Krad)./Drad;
%           volfrac   = Krad./(Drad.^2 + Krad);
      elseif 0
          LM = 0:2:(lmax*2-1); LM = (2*LM+1);%.*exp(-0.05*LM);
          LM(2) = LM(2);
          t = 0:0.05:1;
          Mm = myleg(lmax*2-1,t'); Mm = Mm(:,1:2:end); Mm = Mm.*repmat(LM,[size(Mm,1) 1]);
          S = permute(cmult(vfsm,Mm,4,2),[1 2 3 5 4]);
          S = log(abs(cat(4,ones(size(S,1),size(S,2),size(S,3),1,size(S,5)),S)));

          [B T] = ndgrid(B,t.^2); 
          B=B(:); T=T(:); W = sqrt(1-T(:));
          Ra = B-B.*T;
          Ax = B.*T;

          P = [Ra Ax ];
          Pinv = inv(P'*diag(W)*P)*P'*diag(W);
          A = cmult(S(:,:,:,:),Pinv,4,2);

          Dax = -A(:,:,:,3-1);
          Drad = -A(:,:,:,2-1);

          sg = -1;
          Dax_int = A(:,:,:,1)*0;
          Dax_ext = -A(:,:,:,3-1);
          Drad_ext = -A(:,:,:,2-1);
          volfrac   = A(:,:,:,1)*0;
      
      end;
      
      
      
           
      % get lowres map with #voxels in highres map
      acpervox = imfilter(double(vf(:,:,:,1)>0),ker); 
      acpervox = acpervox(st:osamp:end,st:osamp:end,st:osamp:end,:);
      
%      vfsm = vfsm ./ (eps+repmat(acpervox,[1 1 1 size(vfsm,4)]));
      
      % mean of signal 
      vfmsig = zeros(size(ftr.user.b0avg)); 
      vfmsig(idxmask) = ftr.user.meansignal;
      sz = size(vfsm);
      if any(size(vfmsig) ~= sz(1:3)),
          fac = sz(1)/size(vfmsig,1);
          tmp = vfmsig;
          vfmsig = zeros(sz(1:3));
          for k = 1:fac,
              for j = 1:fac,
                  for i = 1:fac,
                    vfmsig(k:fac:end,j:fac:end,i:fac:end) = tmp;
                  end;
              end;
          end;          
      end;
      
   
          
      
      % mask out with nans
%       maski = vfsm(:,:,:,1)>0;
%       for k = 1:4,
%           tmp = vfsm(:,:,:,k);
%           tmp(not(maski)) = nan;
%           vfsm(:,:,:,k) = tmp;
%       end;
      

      % compute average segment count in voxel
      sz = size(vfsm);
      pos = ftr.user.P(1:3,:)';
      pos = pos ./ repmat(ftr.vox(:)',[size(pos,1) 1]);
      fdsm = hist4(pos,sz(1:3));
      fdsm = fdsm ./ (acpervox+eps);
            
      
      % compute average endpoint count in voxel
      idx = ftr.user.P(9,:) == -1 | ftr.user.P(10,:) == -1;
      ep = hist4(pos(idx,:),sz(1:3));      
      ep = ep ./ (acpervox+eps);

      %V = double(mask); V(mask>0) = Dax;
      maps.vfi = fdsm*0; %est(:,:,:,1);%vfsm(:,:,:,1);
      %V = double(mask); V(mask>0) = Drad;
      maps.vfsw = fdsm*0; %est(:,:,:,3)+est(:,:,:,2);%$vfsm(:,:,:,2);
      %V = double(mask); V(mask>0) = Kax;
      maps.vfe = fdsm*0; %est(:,:,:,3);%vfsm(:,:,:,3);
      %V = double(mask); V(mask>0) = Krad;
      maps.S0 = fdsm*0; % est(:,:,:,4);%vfsm(:,:,:,4);
      maps.segcount = fdsm;
      maps.termcount = ep;
      
      
      maps.Din = fdsm*0; %est(:,:,:,1);%vfsm(:,:,:,1);
      %V = double(mask); V(mask>0) = Drad;
      maps.Dexax = fdsm*0; %est(:,:,:,3)+est(:,:,:,2);%$vfsm(:,:,:,2);
      %V = double(mask); V(mask>0) = Kax;
      maps.Dexrad = fdsm*0; %est(:,:,:,3);%vfsm(:,:,:,3);
      %V = double(mask); V(mask>0) = Krad;
      maps.vf = fdsm*0; %est(:,:,:,4);%vfsm(:,:,:,4);
      maps.vf_csf =fdsm*0; % est(:,:,:,5);%vfsm(:,:,:,4);
      maps.snr = fdsm*0; %snr;%vfsm(:,:,:,4);
      maps.segcount = fdsm;
      maps.termcount = ep;
      
      
%       est(est<0) = 0;
%      assignin('base','est',est);
   
      
      
function h = hist4(d,sz)
d = double(floor(d)+1);

d = d(d(:,1)>=1 & d(:,1)<=sz(1) & d(:,2)>=1 & d(:,2)<=sz(2) & d(:,3)>=1 & d(:,3)<=sz(3) ,:);

didx = sub2ind(sz,d(:,1),d(:,2),d(:,3));
h = sparse([didx(:);sz(1)*sz(2)*sz(3)],ones(length(didx(:))+1,1),ones(length(didx(:))+1,1));
h = reshape(full(h),sz);

      


      

function p = myleg(n,x);
if n == 0,
    p = x*0+1;
    return;
end
if n == 1
    p = [x*0+1 x];
    return;
end;
p = zeros(size(x,1),n+1);
p(:,1:2) = [x*0+1 x];
for k = 2:n,
    p(:,k+1) = ((2*k-1)*x.*p(:,k) - (k-1)*p(:,k-1))/k;    
end;



    
    
    
            
            
            
      
      
      
      