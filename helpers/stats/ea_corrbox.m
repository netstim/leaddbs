function [h,R,p,g] = ea_corrbox(varargin)
% Entry function deciding to do a glmplot, corrplot or raincloud plot
% ea_glmplot input: (X,Y,permutation,labels,group1,group2,colors,markers)
% ea_raincloud input: (X,Y,labels)
% ea_corrplot input: (X,Y,permutation,labels,group1,group2,colors,markers)

if size(varargin{1},2) > 1
    [h,R,p,g] = ea_glmplot(varargin{:});
    return
end

I = varargin{1}(~isnan(varargin{1}));
Ihat = varargin{2}(~isnan(varargin{2}));
if ea_isbinary(I) || ea_isbinary(Ihat)
    if ea_isbinary(I)
        [h,R,p,g] = ea_raincloud(I, Ihat, varargin{4});
    elseif ea_isbinary(Ihat)
        [h,R,p,g] = ea_raincloud(Ihat, I, varargin{4});
    end
else
    [h,R,p,g] = ea_corrplot(varargin{:});
end
