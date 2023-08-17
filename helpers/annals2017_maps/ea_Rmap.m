function [R,Rperm,Rnaned,Rpermnaned,h]=ea_Rmap(varargin)
% creates a correlation nifti file given a set of images and a
% regressor ("R map" in Horn 2017 AoN)

% ea_Rmap(fis,regressor,outputname,mask,sk,corrtype) % sk can be 'k','s','sk'
% for smoothing and normalization options. Mask only necessary if
% choosing 'k' option.

pthresh=0.05;
itercount=1000;
regressor=varargin{2};
output=varargin{3};
h=nan;

if nargin<6
    corrtype='Spearman';
else
    corrtype=varargin{6};
end

[X,n]=ea_genX(varargin{1:5}); % accumulates all images into image matrix X.
X=X';
nnanix=~isnan(nansum(X,1)).*(abs(nansum(X,1)))>0;

if nargin>6
    if ismember(varargin{7},{'permute','permuteplot'})
        Rperm=nan(1000,size(X,2));
        regressorperm=repmat(regressor,1,itercount);
        for i=1:itercount
            regressorperm(:,i)=regressorperm(randperm(numel(regressor)),i);
        end
        Rperm(:,nnanix)=corr(regressorperm,X(:,nnanix),'type',corrtype,'rows','pairwise');
    end
end

mask=varargin{4};
switch lower(corrtype)
    case {'skippedpearson','skipped_pearson','robustpearson','robust_pearson','skipped'}
        R=nan(1,size(X,2));
        R(mask)=ea_skipped_correlation(regressor,X(:,mask),'pearson');
    case {'skippedspearman','skipped_spearman','robust','robust_spearman','robustspearman'}
        R=nan(1,size(X,2));
        R(mask)=ea_skipped_correlation(regressor,X(:,mask),'spearman');
    case {'bend'}
        R=nan(1,size(X,2));
        R(mask)=ea_bendcorr(regressor,X(:,mask));
    case {'circular'}
        R=nan(1,size(X,2));
        R(mask)=ea_circcorr(regressor,X(:,mask));
    case {'stable'}
        R=nan(1,size(X,2));
        rob_r=R;
        T=R;
        sig_T=R;
        
        % run stable R
        [R(mask),rob_r(mask),T(mask),sig_T(mask)]=ea_stablecorr(regressor,X(:,mask)); 
        
%         % already export additional maps since these are special for stable:
%         params=varargin(1:5);
%         [pth,fn,ext]=fileparts(params{3});
%         params{3}=fullfile(pth,[fn,'_robavg',ext]);
%         ea_exportmap(n,rob_r,params{:});
%         
%         params{3}=fullfile(pth,[fn,'_T',ext]);
%         ea_exportmap(n,T,params{:});
%         
%         params{3}=fullfile(pth,[fn,'_sig_T',ext]);
%         ea_exportmap(n,sig_T,params{:});
    otherwise
        R=corr(regressor,X,'type',corrtype,'rows','pairwise');
end
ea_exportmap(n,R,varargin{1:5});

if exist('Rperm','var') % permutation test
    Rpermnaned=[R;Rperm]; % for now, first entry is the unpermuted one.
    sRd=sort([R;Rperm],1,'descend');
    % delete values from Rpermnaned that are not significant (uncorrected):
    delp=Rpermnaned<...
        repmat(sRd(round((pthresh/2)*itercount),:),itercount+1,1);
    deln=Rpermnaned>...
        repmat(sRd(round((1-(pthresh/2))*itercount),:),itercount+1,1);
    del=logical(delp.*deln);
    Rpermnaned(del)=nan;

    Rnaned=Rpermnaned(1,:);
    Rpermnaned=Rpermnaned(2:end,:);

    [pth,fn,ext]=fileparts(varargin{3});
    varargin{3}=fullfile(pth,[fn,'_sig',ext]);
    ea_exportmap(n,Rnaned,varargin{1:5});

    if strcmp(varargin{7},'permuteplot')
        ea_dispercent(0,'Iterating patients');
        for pt=1:length(varargin{1})
            Xthispt=X(pt,:)';
            Ihat(pt)=atanh(corr(Xthispt(varargin{4}),R(varargin{4})','rows','pairwise','type',corrtype)); % real predictions
            Ihat_Rperm(pt,:)=atanh(corr(Xthispt(varargin{4}),...
                Rperm(:,varargin{4})','rows','pairwise','type',corrtype)); % permuted predictions
            ea_dispercent(pt/length(varargin{1}));
        end
        ea_dispercent(1,'end');

        [R]=corr(Ihat',regressor); % Predictive R of unpermuted values

        [Rperm]=corr(Ihat_Rperm,regressor); % Predictive R of permuted values

        Rperm=sort(Rperm,'descend'); % build distribution
        RlargerRperm=R>Rperm;
        p_predict_perm=sum(~RlargerRperm)./size(Rperm,1); % calculate final permutation based p value
        disp(['Permutation based p for overall prediction = ',num2str(p_predict_perm),'.']);

        h=ea_corrplot(regressor,Ihat',p_predict_perm,{'Empirical vs. Predicted','Empirical','Predicted'});
    end
end
