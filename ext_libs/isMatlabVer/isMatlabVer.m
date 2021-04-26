function Dummy = isMatlabVer(varargin)
% isMatlabVer - Compare Matlab version to specified number
% Match = isMatlabVer(Relop, N)
% INPUT:
%   Relop: Comparison operator as string: '<', '<=', '>', '>=', '=='.
%   N:     Number to compare with as DOUBLE vector with 1 to 4 elements.
%
% OUTPUT:
%   Match: Locical scalar, TRUE for matching comparison, FALSE otherwise.
%
% EXAMPLES:
%   version ==> '7.8.0.342 (R2009a)'  (different results for other version!)
%   isMatlabVer('<=', 7)               % ==> TRUE
%   isMatlabVer('>',  6)               % ==> TRUE
%   isMatlabVer('<',  [7, 8])          % ==> FALSE
%   isMatlabVer('<=', [7, 8])          % ==> TRUE
%   isMatlabVer('>',  [7, 8, 0, 342])  % ==> FALSE
%   isMatlabVer('==', 7)               % ==> TRUE
%   isMatlabVer('==', [7, 10])         % ==> FALSE
%   isMatlabVer('>',  [7, 8, 0])       % ==> FALSE (the 342 is not considered!)
%
% NOTES: The C-Mex function takes about 0.6% of the processing time needed
%   by Matlab's VERLESSTHAN, which can check other toolboxes also.
%
% Compile with: mex isMatlabVer.c
%
% Compiler: GCC 9.3.0, Clang 12.0, MSVC 2019 (M/T)
% Author:   Jan Simon, Heidelberg, (C) 2010 matlab.THISYEAR(a)nMINUSsimon.de
% License:  BSD - use, copy, modify on own risk, mention the author.
%
% See also: VER, VERLESSTHAN.

% Rev 14-April-2021 (Ningfei Li, ningfei.li(a)gmail.com):
%   Clean up.
%   Use "u" instead of "L" for string literal (char16_t). 
%   Update compiled binaries.

% This is a dummy M-file to feed the HELP command.
error(['JSimon:', mfilename, ':MexNotFound'], ...
   'Cannot find Mex file. Please compile it!');
