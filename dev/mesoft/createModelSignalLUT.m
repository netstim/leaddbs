
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% creates the c-code of the modelsignal correlation
% and returns the matrix needed to map the signal onto approaximation
% coeffcients

function Bm = createModelSignalLUT(ten,Pstruc,fixedTissuePrior,sinterp,fid,verbose)


maxorder = Pstruc.ordermax(3);
pp = Pstruc.ordermax(4);
alpha = Pstruc.alpha;


% the discrete direction on the sphere we along which we do the approximation
dirs = sinterp.bDir;

% the grid we are optimiznig on
D1 = (-0.2:0.1:3.5)';
D2 = -0.2:0.1:3.5;
d1 = repmat(D1,[1 size(D2,2)]);
d2 = repmat(D2,[size(D1,1) 1]);
[P fun basis] = createSymPolMat([1./(pp+([d1(:) d2(:)]))]',[1:1:maxorder]);

W = createWeightingScheme(ten,Pstruc.b_weighting);

% create all possible exp-model functions on the grid 
E = repmat(permute(squeeze(sum(cmult(dirs,ten,1,1).*repmat(dirs',[1 1 size(ten,3)]),2)),[3 4 1 2]),[size(D1,1) size(D2,2)]);
Erad = repmat(permute(ten(1,1,:)+ten(2,2,:)+ten(3,3,:),[1 2 4 3]),[size(D1,1), size(D2,2) size(dirs,2) 1]);
Am = exp(-(E).*repmat(D1,[1 size(D2,2) size(dirs,2) size(ten,3)])-(Erad-E).*repmat(D2,[size(D1,1) 1 size(dirs,2) size(ten,3)]));
Am = Am - alpha *repmat(sum(cmult(Am,W,4,1),4),[1 1 1 size(Am,4)]);

% the inverse
pinvP = pinv(P');

% and finally the matrix
Bm = reshape(pinvP*reshape(Am,[size(Am,1)*size(Am,2) size(Am,3)*size(Am,4)]),[size(pinvP,1) size(Am,3) size(Am,4)]); %/size(ten,3);
Bm = cmult(Bm,W,3,1);
if not(isempty(fixedTissuePrior)),
    Am = createModelSignalLUT_classical(ten,1,fixedTissuePrior,sinterp,Pstruc,fid);
    Bm = cat(1,Bm,permute(Am,[3 1 2]));
end;


createCcode(fid,basis,pp,maxorder,fixedTissuePrior)



if verbose,

     n1 = [1 1 3]'; n1 = n1 / norm(n1);
     n2 = [0.5 1 1]'; n2 = n2 /norm(n2);
     for k = 1:size(ten,3),
        S(k) = exp(-1*n1'*ten(:,:,k)*n1 - 1*(trace(ten(:,:,k))-n1'*ten(:,:,k)*n1)) + 0*exp(-2*n2'*ten(:,:,k)*n2 - 0.1*trace(ten(:,:,k)));
     end;

     nz = single(randn(size(S))*0.0);
     S = abs(S + nz)';

     alpha = cmult(Bm,S,3,1);

     cost = @(x) modelmodelCorrTest(x(1),x(2),x(1),x(2),1) - 2*modelsignalCorr(x(1),x(2),fun,alpha,pp,sinterp,x(3:5));
      D = 0:0.2:3; 

      for k = 1:length(D), 
          for j = 1:length(D),
          c(k,j) = double(cost([D(k) D(j) n1'])); 
          end;
      end; 
      figure(99);
      imagesc(D,D,exp(-c(:,:,1)*1000));

end;    
    


function r = modelsignalCorr(Dpara,Dorth,fun,alpha,pp,sinterp,n)
   n = n/norm(n);
   M = fun(1/(pp+Dpara),1/(pp+Dorth));
   M = cat(1,M{:});   
   a = evalSinterp(permute(single(repmat(n(:),[1 size(alpha,1)])),[1 3 2]),single(alpha'),sinterp);
   r = a*[M;1];

   
%%%%%%%% creates C-code with estimate approximation coeffcients
function createCcode(fid,basis,pp,order,fixedTissuePrior)


fprintf(fid,'REAL EnergyComputerBase::modelSignalCorr(REAL Dpara,REAL Dorth,REAL *alpha)\n {\n');

fprintf(fid,'\tREAL x = 1.0/(%.5f+ Dpara); \n',pp);
fprintf(fid,'\tREAL y = 1.0/(%.5f+ Dorth); \n',pp);

W = {'x','y'};

for k = 2:order,
    for j = 1:2,
        if k == 2,
            last{j} = W{j};
        end;
        str = ['REAL ' W{j} num2str(k) ' = ' last{j} ' * ' W{j} '; '];
        fprintf(fid,'\t%s \n',str);
        last{j} = [W{j} num2str(k)];
    end;
end;

fprintf(fid, '\treturn  ');
for k = 1:length(basis),
    fprintf(fid,'\n \t     + alpha[%i] * ( %s ) ',k-1,strrep(strrep(basis{k},'.^',''),'.',''));
end;
fprintf(fid,'\n \t     + alpha[%i] ',length(basis));
fprintf(fid,';\n }\n');


fprintf(fid,'REAL EnergyComputerBase::modelSignalCorrGuide(REAL *alpha)\n {\n');

if not(isempty(fixedTissuePrior)),
    fprintf(fid,'\n \t return alpha[%i]*prior_strength ',length(basis)+1);
end;
fprintf(fid,';\n }\n');
   
   
   

%%%%%%%% creates second order mononomials
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





    