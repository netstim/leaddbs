function [h,R,p,g]=ea_corrbox(varargin)
% (X,Y,labels,corrtype,group1,group2,pperm,colors,markers)
% decider function deciding to do a glmplot / corrplot or raincloud plot for results
% visualization

if size(varargin{1},2)>1
    [h,R,p,g]=ea_glmplot(varargin{:});
    return
end

if ea_isbinary(varargin{1}) || ea_isbinary(varargin{2})
    if ea_isbinary(varargin{1})
        [h,R,p,g]=ea_raincloud(varargin{:});
    elseif ea_isbinary(varargin{2})
        [h,R,p,g]=ea_raincloud(varargin{2},varargin{1},varargin{3:end});
    end
else
    [h,R,p,g]=ea_corrplot(varargin{:});
end



