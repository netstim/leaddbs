function varargout=ea_corr(X,Y,corrtype)
% wrapper for different correlation types
if ~exist('corrtype','var')
    corrtype='pearson';
end

% Set inf to nan
X(isinf(X)) = nan;
Y(isinf(Y)) = nan;

switch lower(corrtype)
    case 'pearson'
        [varargout{1},varargout{2}]=corr(X,Y,'rows','pairwise','type','Pearson');
    case 'spearman'
        [varargout{1},varargout{2}]=corr(X,Y,'rows','pairwise','type','Spearman');
    case 'kendall'
        [varargout{1},varargout{2}]=corr(X,Y,'rows','pairwise','type','Kendall');
    case 'bend'
        if (any(isnan(X(:))) && size(X,2)>1) || (any(isnan(Y(:))) && size(Y,2)>1) % feed in column wise
            for xx=1:size(X,2)
                for yy=1:size(Y,2)
                    [R(xx,yy),p(xx,yy)]=ea_bendcorr(X(:,xx),Y(:,yy));
                end
            end
            varargout{1}=R;
            varargout{2}=p;
        else
        [varargout{1},varargout{2}]=ea_bendcorr(X,Y);
        end
    case {'skipped spearman','skipped'}
        [varargout{1}]=ea_skipped_correlation(X,Y,'Spearman');
    case 'skipped pearson'
        [varargout{1}]=ea_skipped_correlation(X,Y,'Pearson');
    case {'dist','distance'}
        [varargout{1}]=ea_distcorr(X,Y);
end

