function trainmodel
%%
%scheme = evalin('base','scheme');
d1 = 0.5+(0:0.05:1)*2.5;
d2 = 0+(0:0.05:1)*3;
d3 = 0.15+(0:0.05:1)*2.5;
v = 0:0.05:1;
[D1 D2 D3 V] = ndgrid(d1,d2,d3,v);
nt = 2*10^5;
D1 = 0.5+rand(nt,1)*2.5;
D2 = 0.5+rand(nt,1)*2.5;
D3 = 0+sqrt(rand(nt,1))*2.5;
V = rand(nt,1)*1;
Vw = (1-V).*rand(nt,1)*0;

%idx = find(abs(D1-(D2+3*D3))<1.5);
idx = find(abs(D1-(D2+D3))<2.5 & abs(D1-(D2+3*D3))<2.5);
D1 = D1(idx);
D2 = D2(idx);
D3 = D3(idx);
V = V(idx);
Vw = Vw(idx);
length(idx)

tries = length(D1(:));

n = randn(3,tries);
n = n./repmat(sqrt(sum(n.^2)),[3 1]);


%ds = get(findobj('tag','fiberGT_main'),'userdata');
%ten = ds.original_bTensor;
ten = evalin('base','hr.user.bTensor')/1000;
for k = 1:size(ten,3),

    
                       [U D] = eigs(ten(:,:,k));
               [~,ix]=sort(D(logical(eye(length(D)))));
                   scheme(:,k) = sqrt(D(1,1))*U(:,ix(3));   

end;

pb = (n'*scheme).^2;
b = sum(scheme.^2);

nz = 0.01;
N = 1;
B0 = (0.25:0.05:1).^4;


S1 =repmat(V(:),[1 size(pb,2)]).*exp(-repmat(D1(:),[1 size(pb,2)]).*pb) + ...
    repmat((1-V(:)-Vw(:)),[1 size(pb,2)]).*exp(-repmat(D2(:),[1 size(pb,2)]).*pb-repmat(D3(:),[1 size(pb,2)]).*repmat(b,[size(pb,1) 1]) ) +...
    repmat(Vw(:),[1 size(pb,2)]).*exp(-repmat(D1(:)*0+3,[1 size(pb,2)]).*repmat(b,[size(pb,1) 1])) ;



for j = 1:length(B0),
    

    S0 = B0(j)*S1;
    S = 0;
    for k = 1:N,
       S = S + abs(S0 + nz*sqrt(N)*(randn(size(S0))+1i*randn(size(S0)))).^2;
    end;
    S = sqrt(S/N);
    S = S/B0(j);








    M = [];

    buni = unique(round(b*10))/10;
    lmax = 6;



    % 
    % 
    % n = n+0.2*randn(3,tries);
    % n = n./repmat(sqrt(sum(n.^2)),[3 1]);
    % pb = (n'*scheme).^2;
    % b = sum(scheme.^2);


    for k = 2:length(buni),
        bval = buni(k);
        idx = find(round(b*10)==bval*10);
        c2 = pb(:,idx)/bval;
        m = myleg(lmax,sqrt(c2(:)))/length(idx);
        m = reshape(m(:,1:2:(lmax+1)),[size(c2,1) size(c2,2) lmax/2+1]);
        legproj = squeeze(sum(repmat(S(:,idx),[1 1 lmax/2+1]).*m,2));
        M = cat(3,M,legproj);
    end;

    P = prepFeats(M,lmax);

    p = 0.8;
    nl = @(x) x; %power(abs(x),p);
    inl = @(x) x; %power(abs(x),1/p);

    target = nl([D1(:) D2(:) D3(:) V(:) Vw(:)]);
    alpha = pinv(P)*target;
    pred = inl(P*alpha);
    ran = [3.5 3 3 1 1];
    sfigure(1000);
    for a = 1:5,
        subplot(2,3,a);
     imagesc(hist3([pred(:,a),target(:,a)],{(0:0.03:1)*ran(a) (0:0.03:1)*ran(a)}))
     colormap hot
    end;
    drawnow;
    ALPHA{j} = alpha;
end;
%assignin('base','alpha',alpha);


modelstruc.lmax = lmax;
modelstruc.maxnz = nz;
modelstruc.prepFeats = @prepFeats;
modelstruc.apply = @(x,y) apply(B0,x,y,ALPHA,lmax,nz);


assignin('base','modelstruc',modelstruc);


%%
return;



function est = apply(B0,SNR,x,ALPHA,lmax,nz)

S0 = SNR*nz;
P = prepFeats(x,lmax);
B0edge = [-inf (B0(2:end)+B0(1:end-1))*0.5 inf];
est = zeros(size(P,1),size(ALPHA{1},2));
for k = 1:length(B0),
    idx = S0>B0edge(k) & S0<=B0edge(k+1);
    est(idx,:) = P(idx,:)*ALPHA{k};
end;










function P = prepFeats(M,lmax)


%    kappa = 0.5;
%    tmp = squeeze(M(:,2:(lmax/2+1),1:end)./repmat(sum(kappa*M(:,(1:lmax/2)*0+1,:) + M(:,2:(lmax/2+1),:),3),[1 1 size(M,3)]));
    tmp = squeeze(M(:,2:(lmax/2+1),1:end-1)./repmat(sum(abs(M(:,2:(lmax/2+1),:)),3),[1 1 size(M,3)-1]));
    M = [squeeze(M(:,1,:)) tmp(:,:) ];

    M = abs(M);
    size(M);
    maxN = size(M,2);
    P = [M(:,1)*0+1 ];
    M = (M);
    for k = 1:maxN,
        P = cat(2,P,M(:,k));
        for j = k:maxN,
            P = cat(2,P,M(:,k).*M(:,j));
            for r = j:maxN,
                P = cat(2,P,M(:,k).*M(:,j).*M(:,r));
            end
        end;        
    end;












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



    
    
    
            
            
            
            
            
            
            