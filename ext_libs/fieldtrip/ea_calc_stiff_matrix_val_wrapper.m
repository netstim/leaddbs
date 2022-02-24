function varargout = ea_calc_stiff_matrix_val_wrapper(varargin)
% Wrapper for 'ea_calc_stiff_matrix_val'. Fix the runtime path before
% calling it. No need to install addtional dependencies.

ea_fix_runtimepath;
[varargout{1:nargout}] = ea_calc_stiff_matrix_val(varargin{:});
ea_delete([pwd, filesep, 'fort.6']);
