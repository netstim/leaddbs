function ea_javacomponent(varargin)
% Warpper for javacomponent, warning disabled

warnState = warning('off', 'MATLAB:ui:ea_javacomponent:FunctionToBeRemoved');
javacomponent(varargin{:})
warning(warnState);
