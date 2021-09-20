function varargout=ea_normalize_apply_normalization(options)
% Wrapper function to apply normalization to pre- an post-op images
%
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn


if ischar(options) % return name of method.
    varargout{1}='(Re-)apply (priorly) estimated normalization';
    varargout{2}=1; % dummy output
    varargout{3}=0; % hassettings.
    varargout{4}=0; % dummy output
    return
end

ea_apply_normalization(options);
