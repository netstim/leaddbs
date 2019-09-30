function ea_javacomponent(varargin)
% Warpper for javacomponent, warning disabled

warnState = warning('off', 'MATLAB:ui:javacomponent:FunctionToBeRemoved');
javacomponent(varargin{:});
warning(warnState);
