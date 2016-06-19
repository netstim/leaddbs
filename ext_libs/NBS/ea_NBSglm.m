function test_stat=NBSglm(varargin)
%NBSglm  Permutation-based non-parametric inference using the general
%linear model (GLM). Vectorial implementation, in the sense that a for loop 
%is not used for each seperate GLM. 
%
%   Test_Stat=NBSglm(GLM) operates on each GLM defined by the structure
%   GLM. 
%
%   Test_Stat=NBSglm(GLM,H) attempts to write out progress to uiwaitbar 
%   with handle H.
%
%   A GLM structure contains the following fields:
%       GLM.perms:        Number of permutations to generate
%       GLM.X:            n x p design matrix, including a column of ones
%                         if necessary. p is the number of independent
%                         variables, n is the number of observations. 
%       GLM.y:            n x M array, where each column stores the 
%                         dependent variables for a seperate GLM
%       GLM.contrast:     1 x p contrast vector, only one contrast can be
%                         specified
%       GLM.test:         Type of test: {'onesample','ttest','ftest'}
%       GLM.exchange:     n x 1 vector specifying valid exchange blocks for 
%                         a repeated measures design. An exchange block of
%                         [1 2 3 1 2 3] confines permutation to
%                         observations 1 & 4, 2 & 5 and 3 & 6 [optional]
%
%   See also GLMexamples
%
%   Test_Stat is a GLM.perms + 1 x M array of test-statistics. The first 
%   row is the oberved test statistic. The remaining rows are samples of 
%   the test-statistic under the null hypothesis. Each column corresponds 
%   to a seperate GLM. 
%
%   Remarks: 
%       A column of ones is always included in the F-statistic reduced 
%       model, unless doing so results in a rank deficient design matrix.
%       However, a column of ones is not included in the full model unless
%       specified by the user in GLM.X. 
%
%       Any independent variable with a zero contrast is treated as a 
%       nuisance regressor. The method of Freedman & Lane is used to deal 
%       with nuisance regressors. This method is described in Anderson & 
%       Robinson (2001) Permutation tests for linear models. 43(1):75-88
%
%       The same permutation sequence is applied to every GLM. This is
%       important if spatial dependencies between each GLM are used in 
%       subsequent cluster-based inference. 
%       
%   azalesky@unimelb.edu.au

GLM=varargin{1};
if nargin==2
    H=varargin{2};
end

%Number of predictors (including intercept)
p=length(GLM.contrast);

%Number of independent GLM's to fit
M=size(GLM.y,2);

%Number of observations
n=size(GLM.y,1);

%Determine nuisance predictors not in contrast
ind_nuisance=find(~GLM.contrast);

if isfield(GLM,'exchange')
    %Set up exchange blocks
    blks=unique(GLM.exchange); 
    %Number of blocks
    n_blks=length(blks);
    %Number of observations per block
    sz_blk=n/n_blks; 
    blk_ind=zeros(n_blks,sz_blk);
    for i=1:length(blks)
        blk_ind(i,:)=find(blks(i)==GLM.exchange);
    end
end

if isempty(ind_nuisance)
    %No nuisance predictors
    
else
    %Regress out nuisance predictors and compute residual
    b=zeros(length(ind_nuisance),M);
    resid_y=zeros(n,M); 
    b=GLM.X(:,ind_nuisance)\GLM.y;
    resid_y=GLM.y-GLM.X(:,ind_nuisance)*b; 
end


test_stat=zeros(GLM.perms+1,M);
for i=1:GLM.perms+1
    y_perm=zeros(n,M);
    %Permute
    if i==1
        %Don't permute the first run
        y_perm=GLM.y; 
    else
        if isempty(ind_nuisance) 
            %Permute signal 
            if exist('blk_ind','var')
                %Use the same permutation for every GLM
                for j=1:n_blks
                    y_perm(blk_ind(j,:),:)=...
                    GLM.y(blk_ind(j,randperm(sz_blk)),:);
                end
            else                
                %Use the same permutation for every GLM
                y_perm=GLM.y(randperm(n)',:);
                %Use a different permutation for every GLM
                %[tmp,perm_vec]=sort(rand(n,M)); clear tmp;
                %GLM.y=GLM.y(perm_vec); 
            end
        else
            %Permute residuals
            if exist('blk_ind','var')
                %Use the same permutation for every GLM
                for j=1:n_blks
                    resid_y(blk_ind(j,:),:)=...
                    resid_y(blk_ind(j,randperm(sz_blk)),:);
                end
            else
                %Use the same permutation for every GLM
                resid_y=resid_y(randperm(n)',:);
                %Use a different permutation for every GLM
                %[tmp,perm_vec]=sort(rand(n,M)); clear tmp;
                %resid_y=resid_y(perm_vec); 
            end
        end
    end
    %Freedman & Lane
    %Add permuted residual back to nuisance signal, giving a realisation 
    %of the data under the null hypothesis  
    if ~isempty(ind_nuisance)
        y_perm=resid_y+[GLM.X(:,ind_nuisance)]*b;
    end
       
    b_perm=zeros(p,M);
    b_perm=GLM.X\y_perm;
    
    %Compute statistic of interest
    if strcmp(GLM.test,'onesample')
        %Flip signs
        if i==1
            %Don't permute first run
            test_stat(i,:)=mean(y_perm); 
        else
            test_stat(i,:)=mean(y_perm.*repmat(sign(rand(n,1)-0.5),1,M)); 
        end
    elseif strcmp(GLM.test,'ttest')
        resid=zeros(n,M);
        mse=zeros(n,M);
        resid=y_perm-GLM.X*b_perm;
        mse=sum(resid.^2)/(n-p);
        se=sqrt(mse*(GLM.contrast*inv(GLM.X'*GLM.X)*GLM.contrast'));
        test_stat(i,:)=(GLM.contrast*b_perm)./se;
    elseif strcmp(GLM.test,'ftest')
        sse=zeros(1,M);
        ssr=zeros(1,M);
        %Sum of squares due to error
        sse=sum((y_perm-GLM.X*b_perm).^2);
        %Sum of square due to regression
        ssr=sum((GLM.X*b_perm-repmat(mean(y_perm),n,1)).^2);
        if isempty(ind_nuisance)
            test_stat(i,:)=(ssr/(p-1))./(sse/(n-p));
        else
            %Get reduced model
            %Column of ones will be added to the reduced model unless the
            %resulting matrix is rank deficient
            X_new=[ones(n,1),GLM.X(:,ind_nuisance)];
            %+1 because a column of 1's will be added to the reduced model
            b_red=zeros(length(ind_nuisance)+1,M);
            %Number of remaining variables
            v=length(find(GLM.contrast))-1;
            [n,ncolx]=size(X_new);
            [Q,R,perm]=qr(X_new,0);
            rankx = sum(abs(diag(R)) > abs(R(1))*max(n,ncolx)*eps(class(R)));
            if rankx < ncolx
                %Rank deficient, remove column of ones
                X_new=GLM.X(:,ind_nuisance);
                b_red=zeros(length(ind_nuisance),M);
                v=length(find(GLM.contrast));
            end
   
%OLD            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
%             if ~any(sum(GLM.X(:,ind_nuisance))==n)
%                 %Add a column of 1's only if one doesn't exist
%                 X_new=[ones(n,1),GLM.X(:,ind_nuisance)];
%                 %+1 because a column of 1's will be added to the reduced model
%                 b_red=zeros(length(ind_nuisance)+1,M);
%                 %Number of remaining variables
%                 v=length(find(GLM.contrast))-1;
%             else
%                 X_new=GLM.X(:,ind_nuisance);
%                 b_red=zeros(length(ind_nuisance),M);
%                 v=length(find(GLM.contrast));
%             end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            sse_red=zeros(1,M);
            ssr_red=zeros(1,M);
            b_red=X_new\y_perm;
            sse_red=sum((y_perm-X_new*b_red).^2);
            ssr_red=sum((X_new*b_red-repmat(mean(y_perm),n,1)).^2);
            test_stat(i,:)=((ssr-ssr_red)/v)./(sse/(n-p));
        end
    end
    try uiwaitbar(H,i/(GLM.perms+1)); catch; end
end

%Added to v1.1.2
%Covers the case where the dependent variable is identically zero for all
%observations. The test statistic in this case in NaN. Therefore, force any
%NaN elements to zero. This case is typical of connectivity matrices
%populated using streamline counts, in which case some regional pairs are
%not interconnected by any streamlines for all subjects. 
test_stat(isnan(test_stat))=0; 
    