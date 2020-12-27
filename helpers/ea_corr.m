function varargout=ea_corr(X,Y,corrtype)
% wrapper for different correlation types
if ~exist('corrtype','var')
    corrtype='pearson';
end

if any(isnan(X)) || any(isnan(Y))
    notnan=(~isnan(X)).*(~isnan(Y));
else
    notnan=ones(length(X),1);
end
if any(isinf(X)) || any(isinf(Y))
    notinf=(~isinf(X)).*(~isinf(Y));
else
    notinf=ones(length(X),1);
end

notnan=logical(notnan.*notinf);

switch lower(corrtype)
    case 'pearson'
        [varargout{1},varargout{2}]=corr(X(notnan),Y(notnan),'type','Pearson');
    case 'spearman'
        [varargout{1},varargout{2}]=corr(X(notnan),Y(notnan),'type','Spearman');
    case 'kendall'
        [varargout{1},varargout{2}]=corr(X(notnan),Y(notnan),'type','Kendall');
    case 'bend'
        [varargout{1},varargout{2}]=ea_bendcorr(X(notnan),Y(notnan));
    case {'skipped spearman','skipped'}
        [varargout{1}]=ea_skipped_correlation(X(notnan),Y(notnan),'Spearman');
    case 'skipped pearson'
        [varargout{1}]=ea_skipped_correlation(X(notnan),Y(notnan),'Pearson');
end

