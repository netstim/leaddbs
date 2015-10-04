function bestpp = createModelModelLUT(ten,Pstruc,fixedTissuePrior,fid,verbose)

ordermax = Pstruc.ordermax(1);
kappa = Pstruc.ordermax(2);
alpha = Pstruc.alpha;


% get the likelihood norm we are optimiznig for
W = createWeightingScheme(ten,Pstruc.b_weighting);




% the grid of diffusion values we are optimizing on
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
dirs = load('dwidirections.mat');
dirs = dirs.dirs128;



%pp = 0.4:0.2:8;
pp = kappa;

order = 1:ordermax;


if length(pp) > 1, % serach for best kappa value

    for k = 1:length(pp),
        M = createTotSymPolMat(1./([L1(:) L2(:) L3(:)]'+pp(k)),order);
        invM = pinv(M');
        S = sum(W*exp(-(squeeze(ten(1,1,:))*L1(:)'+squeeze(ten(2,2,:))*L2(:)'+squeeze(ten(3,3,:))*L3(:)') ),1);

        Sest = (M'*invM*(S(:)));
        err(k) = max(abs(S(:)-Sest(:))./(S(:)));

        if k > 1,
            if err(k) > err(k-1),
              %  break;
            end;
        end;
        sfigure(999); plot(pp(1:length(err)),err); drawnow;
        fprintf('.');
    end;

    idx = (find(err<imerode(err,[1 0 1]))); 
    idx = idx(end);
    
    bestpp = pp(idx);
else  % or take the suggested one
    bestpp = pp;
end;


[M basis basis_eqyz] = createTotSymPolMat(1./([L1(:) L2(:) L3(:)]'+bestpp),order);
S = sum(W*exp(-(squeeze(ten(1,1,:))*L1(:)'+squeeze(ten(2,2,:))*L2(:)'+squeeze(ten(3,3,:))*L3(:)') ),1);
a =  pinv(M')*S(:);
Sest = M'*a;
error = max(abs(S(:)-Sest(:))./(S(:)));


dataPrior = [];
createCcode(fid,basis,basis_eqyz,a,bestpp,order,dataPrior);
createModelModelLUT_classical(ten,1,fixedTissuePrior,Pstruc,fid);


fprintf('#monomials: %i maximal error:%f kappa:%f\n',length(a),error,bestpp);

if verbose,
    figure(10);
    subplot(2,1,1);
    plot(pp(1:length(err)),err);
    subplot(2,1,2);
    plot(a,'*-');

    S = reshape(S,size(L1));
    Sest = reshape(Sest,size(L1));

    figure(3);
    for k = 1:min(n,12),
        subplot(3,4,k);        
        plot(1:n,log(S(:,:,k))-log(Sest(:,:,k)))
    end;
end;


return



function createMfun(fid,basis,alpha,pp,order)

fprintf(fid,'function r = modelmodelCorrTest(Dpara1,Dorth1,Dpara2,Dorth2,dot)\n');

fprintf(fid,'del1 = Dpara1 - Dorth1;\ndel2 = Dpara2 - Dorth2;\n');
fprintf(fid,'sumdel = del1 + del2;\n');
fprintf(fid,'sq = sqrt( sumdel*sumdel - 4*del1*del2*(1-dot*dot) );\n');
fprintf(fid,'x = 1.0/(%.2f+0.5*(sumdel+sq) + Dorth1+Dorth2); \n',pp);
fprintf(fid,'y = 1.0/(%.2f+0.5*(sumdel-sq) + Dorth1+Dorth2); \n',pp);
fprintf(fid,'z = 1.0/(%.2f+Dorth1+Dorth2); \n',pp);

W = {'x','y','z'};

for k = 2:max(order),
    for j = 1:3,
        if k == 2,
            last{j} = W{j};
        end;
        str = [ W{j} num2str(k) ' = ' last{j} ' .* ' W{j} '; '];
        fprintf(fid,'%s \n',str);
        last{j} = [W{j} num2str(k)];
    end;
end;

fprintf(fid, 'r =  %.5f + ', alpha(end));
for k = 1:length(basis),
    fprintf(fid,'... \n       + %.5f * ( %s )  ',alpha(k),strrep(basis{k},'.^',''));
end;
fprintf(fid,';\n');





function createCcode(fid,basis,basis_eqyz,alpha,pp,order,dataPrior)



fprintf(fid,'REAL EnergyComputerBase::modelmodelCorr(REAL Dpara1,REAL Dorth1,REAL Dpara2,REAL Dorth2,REAL dot)\n {\n');
fprintf(fid,'\tif (dot == 1) {\n');
fprintf(fid,'\tREAL x = 1.0/(%.5f+Dpara1+Dpara2); \n',pp);
fprintf(fid,'\tREAL y = 1.0/(%.5f+Dorth1+Dorth2); \n',pp);


W = {'x','y'};

for k = 2:max(order),
    for j = 1:2,
        if k == 2,
            last{j} = W{j};
        end;
        str = ['REAL ' W{j} num2str(k) ' = ' last{j} ' * ' W{j} '; '];
        fprintf(fid,'\t%s \n',str);
        last{j} = [W{j} num2str(k)];
    end;
end;


fprintf(fid, '\tREAL res = %.8f + ', alpha(end));
uniquebasis_eqyz = unique(basis_eqyz);

for k = 1:length(uniquebasis_eqyz),
    alpha_accu = sum(alpha(cellfun(@(x) strcmp(x,uniquebasis_eqyz{k}),basis_eqyz)));
    fprintf(fid,'\n \t     + %.8f * ( %s ) ',alpha_accu,strrep(strrep(uniquebasis_eqyz{k},'.^',''),'.',''));
end;
fprintf(fid,';\n ');


fprintf(fid, '\treturn res;\n');
fprintf(fid,'\t} else {\n');
fprintf(fid,'\tREAL del1 = Dpara1 - Dorth1;\n\tREAL del2 = Dpara2 - Dorth2;\n');
fprintf(fid,'\tREAL sumdel = del1 + del2;\n');
fprintf(fid,'\tREAL sq = sqrt( sumdel*sumdel - 4*del1*del2*(1-dot*dot) );\n');
fprintf(fid,'\tREAL x = 1.0/(%.5f+0.5*(sumdel+sq) + Dorth1+Dorth2); \n',pp);
fprintf(fid,'\tREAL y = 1.0/(%.5f+0.5*(sumdel-sq) + Dorth1+Dorth2); \n',pp);
fprintf(fid,'\tREAL z = 1.0/(%.5f+Dorth1+Dorth2); \n',pp);




W = {'x','y','z'};

for k = 2:max(order),
    for j = 1:3,
        if k == 2,
            last{j} = W{j};
        end;
        str = ['REAL ' W{j} num2str(k) ' = ' last{j} ' * ' W{j} '; '];
        fprintf(fid,'\t%s \n',str);
        last{j} = [W{j} num2str(k)];
    end;
end;

fprintf(fid, '\tREAL res = %.8f + ', alpha(end));
for k = 1:length(basis),
    fprintf(fid,'\n \t     + %.8f * ( %s ) ',alpha(k),strrep(strrep(basis{k},'.^',''),'.',''));
end;
fprintf(fid,';\n ');


fprintf(fid, '\treturn res; } \n }\n\n ');

if not(isempty(dataPrior)),
    fprintf(fid,'REAL EnergyComputerBase::modelmodelCorrFixTis(REAL dot)\n {\n');
    fprintf(fid, '\tREAL res = 0; \n');
    fprintf(fid,'%s\n',dataPrior);
    fprintf(fid,' \treturn res;\n } \n');
end;




function [M strbasis strbasis_eqyz] = createTotSymPolMat(scheme,orders)

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
strbasis = arrayfun(@(x)  strrep(strrep(char(x),'^','.^'),'*','.*'),basis,'uniformoutput',false);
strbasis_eqyz = arrayfun(@(x)  strrep(strrep(char(simplify(subs(x,'z','y'))),'^','.^'),'*','.*'),basis,'uniformoutput',false);
%(basis)
f = @(x,y,z) cellfun(@(c) eval(c),strbasis, 'UniformOutput',false);
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
    
    





