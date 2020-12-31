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
        [varargout{1},varargout{2}]=corr(X,Y,'rows','complete','type','Pearson');
    case 'spearman'
        [varargout{1},varargout{2}]=corr(X,Y,'rows','complete','type','Spearman');
    case 'kendall'
        [varargout{1},varargout{2}]=corr(X,Y,'rows','complete','type','Kendall');
    case 'bend'
        [varargout{1},varargout{2}]=ea_bendcorr(X,Y);
    case {'skipped spearman','skipped'}
        [varargout{1}]=ea_skipped_correlation(X,Y,'Spearman');
    case 'skipped pearson'
        [varargout{1}]=ea_skipped_correlation(X,Y,'Pearson');
end

