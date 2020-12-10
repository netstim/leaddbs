function varargout=ea_corr(X,Y,corrtype)
% wrapper for different correlation types
switch lower(corrtype)
    case 'pearson'
        [varargout{1},varargout{2}]=corr(X,Y,'rows','pairwise','type','Pearson');
    case 'spearman'
        [varargout{1},varargout{2}]=corr(X,Y,'rows','pairwise','type','Spearman');
    case 'kendall'
        [varargout{1},varargout{2}]=corr(X,Y,'rows','pairwise','type','Kendall');
    case 'bend'
        [varargout{1},varargout{2}]=ea_bendcorr(X,Y);
    case 'skipped spearman'
        [varargout{1}]=ea_skipped_correlation(X,Y,'Spearman');
    case 'skipped pearson'
        [varargout{1}]=ea_skipped_correlation(X,Y,'Pearson');
end