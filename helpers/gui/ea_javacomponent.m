function varargout = ea_javacomponent(varargin)
% Warpper for javacomponent, warning disabled

warnState = warning('off', 'MATLAB:ui:javacomponent:FunctionToBeRemoved');
[varargout{1:nargout}] = javacomponent(varargin{:});
warning(warnState);
