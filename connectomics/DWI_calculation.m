function [EigenVal, eigVect, M_error, sqErr, kurtParams, kurtTensor] = do_calculation(DEscheme,  singleSliceDWI_set, matrix,threshold,kurto)

if exist('kurto') ~= 1
   kurto = false;
end;


if length(size(DEscheme)) == 2,
    for k = 1:size(DEscheme,1),
        B_tensor(:,:,k) = DEscheme(k,:)'*DEscheme(k,:);
    end;
else
    B_tensor = DEscheme;
end;

% nsteps =length(B_tensor);
% X1 = [reshape(B_tensor(1,1,:),1,nsteps,1)]; %1st diag element
% X2 = [reshape(B_tensor(2,2,:),1,nsteps,1)]; %2nd diag element
% X3 = [reshape(B_tensor(3,3,:),1,nsteps,1)]; %3rd diag element
% X4 = [reshape(B_tensor(1,2,:),1,nsteps,1)];
% X5 = [reshape(B_tensor(1,3,:),1,nsteps,1)];
% X6 = [reshape(B_tensor(2,3,:),1,nsteps,1)];
% anisoX=[ones(size(X1')), X1', X2', X3',2*X4', 2*X5', 2*X6'];
% iAnisoX=pinv(anisoX);





%%
bval = squeeze(B_tensor(1,1,:)+ B_tensor(2,2,:)+ B_tensor(3,3,:))/1000;


for k = 1:size(B_tensor,3),
    [U D] = eigs(B_tensor(:,:,k),1);
    b_Dir(:,k) = U(:,1)*sqrt(D);;
end;


if length(unique(round(bval*10))) < 3 | not(kurto) ,        
    kurto = false;
 %   monfun = createSymPolMat(2);
 %   M = monfun(b_Dir(1,:),b_Dir(2,:),b_Dir(3,:));
  
    x = b_Dir(1,:); y = b_Dir(2,:); z = b_Dir(3,:);
    M = [ x*0+1; x.^2 ; x.*y*2 ;  y.^2 ;  x.*z*2 ;  y.*z*2 ;z.^2];
    
else
    kurto = true;
    monfun = createSymPolMat(4);
    M = monfun(b_Dir(1,:),b_Dir(2,:),b_Dir(3,:));
    mC = [1 1/3 0 1/3 0 0 1/3 1/5 0 2/5 0 1/5 0 0 0 0 2/5 0 2/5 0 0 1/5]';
end;

if kurto
    display('DKI estimation is performed! Please be patient....')
end;

M = M';
%%
Minv = pinv(M);
A = -M; A(:,1:7)=0*A(:,1:7); 
%W = diag(0.1+1./(1+bval).^4);
%Minv = inv(M'*A*M)*M'*A;

%%



%%   do dti calculation for the particular slice
mask= mean(singleSliceDWI_set,3)>threshold;


%% initiate variables
eigVect=zeros(matrix(1), matrix(2), 3,3);
kurtParams=zeros(matrix(1), matrix(2), 8);
difComp=zeros(matrix(1), matrix(2), 3,3);
M_error= zeros(matrix(1:2));
sqErr= zeros(matrix(1:2));

if nargout == 6,
    kurtTensor = zeros(matrix(1), matrix(2),22);
end;


maxD = 5*10^-3;
maxK = 5;

%% tensor calculation
Y1 =singleSliceDWI_set;
for x_coor =1: matrix(2)
    for y_coor =1: matrix(1)
        Y  = [];	% clear
        if mask(y_coor, x_coor) ~= 0;
            %%
            tempY = Y1(y_coor,x_coor,:);
            tempY(isnan(tempY)) = 0;
                   
            y=squeeze(log(abs(tempY)+eps )); 
            Y = [Y y'];

            estS = exp(M*pinv(M)*y);
                        
            W =diag(double(estS));
            if kurto,
                [aniso_a renorm resid exitflag]=lsqlin(W*M,W*double(y),A,zeros(size(M,1),1),[],[],[],[],[],optimset('largescale','off','display','off'));
            else
                aniso_a= Minv*y;
            end;
            YY = M * aniso_a;
            M_error(y_coor,x_coor) = max(abs(YY - Y')); % MaxErr
            d_tens = [aniso_a(2),aniso_a(3),aniso_a(5);...
                      aniso_a(3),aniso_a(4),aniso_a(6);...
                      aniso_a(5),aniso_a(6),aniso_a(7)]; 
            if any(isinf(d_tens(:))) | any(isnan(d_tens(:))),
                d_tens;
            end;
            [vv, dd] = eig(-d_tens);
            eigVect(y_coor,x_coor,:,:) = vv(:,:);
            difComp(y_coor,x_coor,:,:) = dd(:,:);
            sqErr(y_coor,x_coor) = mean((YY - Y').^2);
            
            if kurto,

                [d idx] = sort(diag(dd));
                d = abs(d);
                pvec = vv(:,idx(3));
                Pvec = monfun(pvec(1),pvec(2),pvec(3));
                D_ax = d(3);
                K_ax =  (6*Pvec(8:end)'*aniso_a(8:end) / D_ax^2);
                
                K_mean = 6*aniso_a(8:end)'*mC(8:end) / (aniso_a(2:7)'*mC(2:7))^2;
                
                
%                 q = 0.25; phi= (q:q:1)*pi;
%                 rvecs = vv(:,idx(1))*cos(phi) + vv(:,idx(2))*sin(phi);
%                 AVorth = monfun(rvecs(1,:),rvecs(2,:),rvecs(3,:));
%                 K_AVvec = -6*mean(AVorth(8:end,:)'*aniso_a(8:end)) / (d(1)+ d(2))^2
                
                q = 0.05; phi= (q:q:1)*pi;
                rvecs = vv(:,idx(1))*cos(phi) + vv(:,idx(2))*sin(phi);
                AVorth = monfun(rvecs(1,:),rvecs(2,:),rvecs(3,:));
                Dorth = -(AVorth(2:7,:)'*aniso_a(2:7));
                [K_orth idx] = max( (6*(AVorth(8:end,:)'*aniso_a(8:end)) ./ Dorth.^2));
                D_orth = Dorth(idx);
                
                sg = 1;
                C_ax = 2*K_ax/6*D_ax^2;
                C_orth = 2*K_orth/6*D_orth^2;
                Dpara_int = (D_ax*C_orth + D_ax*D_orth^2 - D_orth*(sg*abs(C_ax*C_orth)^(1/2) + D_ax*D_orth))/C_orth;
                Dpara_ext = (sg*abs(C_ax*C_orth)^(1/2) + D_ax*D_orth)/D_orth;
                Dorth_ext = (D_orth^2 + C_orth)/D_orth;
                volfrac   = C_orth/(D_orth^2 + C_orth);
                
                
                
                D_orth(D_orth<0) = 0; D_orth(D_orth>maxD) = maxD;
                Dpara_int(Dpara_int<0) = 0; Dpara_int(Dpara_int>maxD) = maxD;
                Dpara_ext(Dpara_ext<0) = 0; Dpara_ext(Dpara_ext>maxD) = maxD;
                Dorth_ext(Dorth_ext<0) = 0; Dorth_ext(Dorth_ext>maxD) = maxD;
                                
                kurtParams(y_coor,x_coor,:) = [K_ax K_orth K_mean D_orth Dpara_int Dpara_ext Dorth_ext  volfrac];
                kurtTensor(y_coor,x_coor,:) = aniso_a;
                
            end;
            
            
            
        end
    end
end

 kurtParams = real(kurtParams);
 kurtParams(kurtParams<0) = 0;
 kurtParams(kurtParams>maxK) = maxK;

EigenVal= sum(difComp, 3);
EigenVal(EigenVal<0) = 0;
EigenVal(EigenVal>maxD) = maxD;

if not(kurto)
    kurtParams = [];
end;


%end of function  do_calculation







% function [f fac] = createSymPolMat(maxorder)
% syms x y z;
% Q = [x y z];
% basis = [];
% fac = [];
% for order = 2:2:maxorder,
%         idx = [1 2 3]';
%         clear w;
%         for k = 1:order-1,
%             idx = [idx ones(size(idx,1),1)*1 ; ...
%                    idx ones(size(idx,1),1)*2 ; ...
%                    idx ones(size(idx,1),1)*3];
%             idx = idx(find(idx(:,end)>=idx(:,end-1)),:);
%         end;
%         for k = 1:size(idx,1),
%             q = hist(idx(k,:),[1 2 3]);
%             w(k) = multinom(order,q); 
%         end;
%         
%         fac = [fac w];
%         basis = [basis ; w'.*prod(Q(idx),2)];
%         %basis = [basis ; prod(Q(idx),2)];
% end;    
% basis = arrayfun(@(x)  strrep(strrep(char(x),'^','.^'),'*','.*'),basis,'uniformoutput',false);
% 
% dec = @(x) cat(1,x{:});
% app = @(x) cat(1,ones(1,size(x,2)),x);
% f = @(x,y,z) cellfun(@(c) eval(c),basis, 'UniformOutput',false);
% f = @(x,y,z) app(dec(f(x,y,z)));
% 
% fac = [1 fac];
% 
% return;


function r = multinom(N,q)
r = 1;
n = N;
for k = 1:length(q),
    if q(k)>0,
        r = r* factorial(n)/factorial(n-q(k))/factorial(q(k));
        n = n-q(k);
    end;
end;




