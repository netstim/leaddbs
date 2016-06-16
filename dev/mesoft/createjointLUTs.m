%test
function Errs = createjointLUTs(ten);


dirs = load('dwidirections.mat');
dirs = dirs.dirs128;

W = createWeightingScheme(ten);
W = diag(W);


orderMM = 8;
orderMS = 4;

alpha = 0;

%% model model

n = 25;
Dmin = -0.3;
Dmax = 6;
D1 = [Dmin Dmax];
D2 = [Dmin Dmax];
D3 = [Dmin Dmax];
l1 = (0:(n-1))/(n-1)*(D1(2)-D1(1)) + D1(1);
l2 = (0:(n-1))/(n-1)*(D2(2)-D2(1)) + D2(1);
l3 = (0:(n-1))/(n-1)*(D3(2)-D3(1)) + D3(1);
[L1 L2 L3] = ndgrid(l1,l2,l3);
idx = L1>=L2 & L2>=L3;
L1 = L1(idx(:));
L2 = L2(idx(:));
L3 = L3(idx(:));


%% model signal

Dpa = (-0.2:0.1:3.5)';
Dor = Dpa';

d1 = repmat(Dpa,[1 size(Dor,2)]);
d2 = repmat(Dor,[size(Dpa,1) 1]);




%% corr ms - mm
[Dpara1 Dorth1 Dpara2 Dorth2] = ndgrid(Dpa(:));
lam1 =  Dpara1+Dpara2;
lam23 = Dorth1+Dorth2;


%kappaMM = 1:0.5:8;
%kappaMS = 1:0.5:8;

kappaMM = 1:0.2:2;
kappaMS = 1:0.2:2;

%kappaMM = 3.5:0.1:4.5;
%kappaMS = 3.5:0.1:4.5;

Errcmat = zeros(length(kappaMM),length(kappaMS));
ErrMM = zeros(length(kappaMM),length(kappaMS));
ErrMS = zeros(length(kappaMM),length(kappaMS));

for kmm = 1:length(kappaMM),
    ppMM = kappaMM(kmm);
    [MS] = createTotSymPolMat(1./([lam1(:) lam23(:) lam23(:)]'+ppMM),1:orderMM);    
    
    for kms = 1:length(kappaMS),
        ppMS = kappaMS(kms);


        %% model model

        [M basis] = createTotSymPolMat(1./([L1(:) L2(:) L3(:)]'+ppMM),1:orderMM);
        S = sum(diag(W)*exp(-(squeeze(ten(1,1,:))*L1(:)'+squeeze(ten(2,2,:))*L2(:)'+squeeze(ten(3,3,:))*L3(:)') ),1);
        a =  pinv(M')*S(:);
        Sest = M'*a;



        %% model signal
        [P fun basis] = createSymPolMat([1./(ppMS+([d1(:) d2(:)]))]',[1:1:orderMS]);
        E = repmat(permute(squeeze(sum(cmult(dirs,ten,1,1).*repmat(dirs',[1 1 size(ten,3)]),2)),[3 4 1 2]),[size(Dpa,1) size(Dor,2)]);
        Erad = repmat(permute(ten(1,1,:)+ten(2,2,:)+ten(3,3,:),[1 2 4 3]),[size(Dpa,1), size(Dor,2) size(dirs,2) 1]);
        Am = exp(-(E).*repmat(Dpa,[1 size(Dor,2) size(dirs,2) size(ten,3)])-(Erad-E).*repmat(Dor,[size(Dpa,1) 1 size(dirs,2) size(ten,3)]));
        Am = Am - alpha *repmat(sum(cmult(Am,W,4,1),4),[1 1 1 size(Am,4)]);
        pinvP = pinv(P');

        Bm = reshape(pinvP*reshape(Am,[size(Am,1)*size(Am,2) size(Am,3)*size(Am,4)]),[size(pinvP,1) size(Am,3) size(Am,4)]); %/size(ten,3);
        Bm = cmult(Bm,diag(W),3,1);


        %% eval 
        Cmat = reshape(a'*MS,[1 1]*length(Dpa)^2);

        % [MRW1] = createTotSymPolMat(1./([Dpara1(:) Dorth1(:) Dorth1(:)]'+ppMM),1:orderMM);
        % [MRW2] = createTotSymPolMat(1./([Dpara2(:) Dorth2(:) Dorth2(:)]'+ppMM),1:orderMM);
        % Cmat = reshape(a'*MS - (a'*MRW1).*(a'*MRW2),[1 1]*length(Dpa)^2);


                %Cmat2 = reshape(cmult(Am(:,:,k,:),cmult(Am(:,:,k,:),diag(W),4,1),4,4),[1444 1444]);
         
         k = 1;
         Cmat2 = cmult(Bm(:,k,:),Am(:,:,k,:),3,4);
         Cmat2 = reshape(cmult(P,Cmat2,1,1),[1444 1444]);
         Cmat3 = reshape(cmult(Am(:,:,k,:),cmult(Am(:,:,k,:),diag(W),4,1),4,4),[1444 1444]);

         errcmat = sum(abs(Cmat(:)-Cmat2(:)))/1444^2*100;
         errMM = sum(abs(S(:)-Sest(:))./S(:))*100;
         errMS =  sum(abs(Cmat3(:)-Cmat2(:)))/1444^2*100;
         
         Errcmat(kmm,kms) = errcmat;
         ErrMM(kmm,kms) = errMM;
         ErrMS(kmm,kms) = errMS;
         fprintf('ppMM %.1f ppMS %.1f  errcmat %f errMM %f errSM %f\n',ppMM, ppMS,errcmat,errMM,errMS);
         
         sfigure(3);
         subplot(1,3,1);
         imagesc(kappaMM,kappaMS,log(Errcmat));
         subplot(1,3,2);
         imagesc(kappaMM,kappaMS,log(ErrMM));
         subplot(1,3,3);
         imagesc(kappaMM,kappaMS,log(ErrMS));
         
    end;
end; 
Errs.cmat = Errcmat;
Errs.MM = ErrMM;
Errs.MS = ErrMS;
Errs.kappaMM = kappaMM;
Errs.kappaMS = kappaMS;

 [mini idx] = min(Errs.cmat(:)); [i j] = ind2sub(size(Errs.cmat),idx); [Errs.kappaMM(i) Errs.kappaMS(j)]
 mini

function [M basis] = createTotSymPolMat(scheme,orders)

syms x y z;

Q = [x; y; z];
basis = [];


for order = orders,
        idx = [1 2 3]';
        clear w;
        for k = 1:order-1,
            idx = [idx ones(size(idx,1),1)*1 ; ...
                   idx ones(size(idx,1),1)*2 ; ...
                   idx ones(size(idx,1),1)*3];
            idx = idx(find(idx(:,end)>=idx(:,end-1)),:);
        end;
       
        B = prod(Q(idx),2);
        
        l = labelSymMonoms(idx);
        clear M;
        for k = 1:max(l),
            M(k,1) = sum(B(find(l==k)));
        end;
        
        basis = [basis ; factor(M)];
end;    
basis = arrayfun(@(x)  strrep(strrep(char(x),'^','.^'),'*','.*'),basis,'uniformoutput',false);
%(basis)
f = @(x,y,z) cellfun(@(c) eval(c),basis, 'UniformOutput',false);
M = (f(scheme(1,:),scheme(2,:),scheme(3,:)));
M = cat(1,M{:});
M = [M ; ones(1,size(scheme,2))];

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
    



   
function [M f basis] = createSymPolMat(scheme,orders)

syms x y;

Q = [x; y];
basis = [];


for order = orders,
        idx = [1 2]';
        clear w;
        for k = 1:order-1,
            idx = [idx ones(size(idx,1),1)*1 ; ...
                   idx ones(size(idx,1),1)*2 ];
                   
            idx = idx(find(idx(:,end)>=idx(:,end-1)),:);
        end;
      
        basis = [basis ; prod(Q(idx),2)];
end;    
basis = arrayfun(@(x)  strrep(strrep(char(x),'^','.^'),'*','.*'),basis,'uniformoutput',false);
basis;
f = @(x,y) cellfun(@(c) eval(c),basis, 'UniformOutput',false);
M = (f(scheme(1,:),scheme(2,:)));
M = cat(1,M{:});
M = [M ; ones(1,size(scheme,2))];

return;






