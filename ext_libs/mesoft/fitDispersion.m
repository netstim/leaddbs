

 
 
function res = fitDispersion(Dpara_in,Dpara_ex,Dorth,vf,vf_csf,snr,noisedeg, proto,Smeas,dirs256,disptypes)


    b = proto.b;
    buni = proto.buni;
    dirs = proto.dirs;
    ten = proto.ten;
 
    
    
    
    
       %% recompute signal
        
         lag = @(x) exp(x/2).*((1-x).*besseli(0,-x/2) - x.* besseli(1,-x/2));
         riceexp = @(x,sigma) sigma*sqrt(pi/2) *lag(-x.^2./(2*sigma^2));        

      
         S = simSig(Dpara_in,Dpara_ex,Dorth,vf,vf_csf,ten,[0 0 1]',1,'nodisp');
        
         highb = round(b/100) == round(buni(end)/100);
         
         lmax = 8;
         [ n weights fod] =  csd(Smeas(highb),dirs(:,highb),dirs256, S(highb),lmax);
  
        weights(weights<0.25) = 0;
        if all(weights==0),
            weights = weights+1; % just to avoid nans
            n = randn(size(n));
        end;
        weights = weights/sum(weights);
        
        Mten = 0;
        for k = 1:size(n,2), 
            w = weights(k);
            Mten = Mten + n(:,k)*n(:,k)' *w;
        end;
        if any(isnan(Mten(:))) | any(isinf(Mten(:))),
            Mten = randn(3,3);
        end;
 
        [U D] = eigs(Mten);
        [~,ix]=sort(D(logical(eye(length(D)))));
        pd = U(:,ix(3));

      %%       
      
      
        lmax = 2;
        kappa = 1;
        S = simSig(Dpara_in,Dpara_ex,Dorth,vf,vf_csf,ten,n,weights,'nodisp');
        S = riceexp(S,1/snr);
        for k = 2:length(buni),
            bidx = round(b/100) == round(buni(k)/100);                  
            [M shidx] = SH(dirs(:,bidx),lmax);
            pM = M*S(bidx); 
            pp = M*Smeas(bidx); 
            for l = 2:length(shidx)
                pMP(l-1,k-1) = pp(shidx{l})'*pM(shidx{l});
                pPred2(l-1,k-1) = sum(pM(shidx{l}).^2);
                pMeas2(l-1,k-1) = sum(pp(shidx{l}).^2);
            end;
        end;
                         
        aa = sum(pMP,2)./sum(pPred2,2);
        kappa = aa(1);
        kappa(kappa>1) = 1;
        kappa(kappa<0) = 0;
        
        

        
        
    
    clear S
    for k = 1:length(disptypes);
        S{k} = simSig(Dpara_in,Dpara_ex,Dorth,vf,vf_csf,ten,n,weights,disptypes{k},kappa);  
    
        errRic(k) = chiLogLik(double(Smeas), double(S{k}),noisedeg*2,double( 1/snr(1)));
        errGS(k) = sum((Smeas-S{k}).^2)*snr(1)^2;

        S{k} = riceexp(S{k},1/snr(1));
    end;
 
    res.S = S;
    res.errRic = errRic;
    res.errGS = errGS;
    res.fod = fod;
    res.weights = weights;
    res.n = n;
    res.pd = pd;
    res.kappa = kappa;
    
 
 
 
      

function [M idx] = SH(n,lmax)
n = n ./ repmat(sqrt(sum(n.^2)),[3 1]);

M = [];
for l = 0:2:lmax,
    m = legendre(l,n(3,:),'sch');
    n2 = n(1,:)+i*n(2,:);
    n2 = n2 ./ abs(n2);
    m = m.* (repmat(n2,[size(m,1) 1]).^repmat((0:l)',[1 size(n2,2)]))*sqrt(2*l+1);
    idx1 = size(M,1);
    M = [M ; m(1,:) ; real(m(2:end,:)) ; imag(m(2:end,:))];
    idx2 = size(M,1);
    idx{l/2+1} = idx1+1:idx2;
end;

M = M/sqrt(size(n,2));


 function S = simSig(Dpara_in,Dpara_ex,Dorth,vf,vf_csf,ten,n,weights,disptype,kappa)
 S = 0;
 tr = squeeze(ten(1,1,:)+ten(2,2,:)+ten(3,3,:))/1000;
 for k = 1:size(n,2),
     %q2 = squeeze(cmult(n(:,k),cmult(n(:,k),ten,1,1),1,2))/1000;
     q2 = squeeze((n(1,k).^2 *ten(1,1,:) + n(2,k).^2 *ten(2,2,:) + n(3,k).^2 *ten(3,3,:) + 2*(n(1,k)*n(2,k)*ten(1,2,:) + n(3,k)*n(2,k)*ten(3,2,:) + n(1,k)*n(3,k)*ten(1,3,:)))/1000); 
     if strcmp(disptype,'nodisp') ~= 1,
         Eintra = dispstick(Dpara_in,kappa,[],[q2(:) (tr(:)-q2(:))]',disptype);
         Eextra = dispstick(Dpara_ex-Dorth,kappa,[],[q2(:) (tr(:)-q2(:))]',disptype) .* exp(-(tr-q2)*Dorth); 
     else
        Eintra = exp(-q2*Dpara_in);
        Eextra = exp(-q2*Dpara_ex-(tr-q2)*Dorth); 
     end;
     Ecsf = exp(-tr*3);
     S = S + weights(k)*(vf*Eintra+ (1-vf-vf_csf)*Eextra + vf_csf*Ecsf );
 end;

 
 
 
 
 