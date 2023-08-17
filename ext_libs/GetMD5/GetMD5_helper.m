function S = GetMD5_helper(V)
% GetMD5_helper: Convert non-elementary array types for GetMD5
% The C-Mex function GetMD5 calls this function to obtain meaningful unique data
% for function handles, java or user-defined objects and sparse arrays. The
% applied processing can depend on the needs of the users, therefore it is
% implemented as an M-function, which is easier to modify than the C-code.
%
% INPUT:
%   V: Array of any type, which is not handled in the C-Mex.
% OUTPUT:
%   S: Array or struct containing elementary types only.
%      The implementation migth be changed by the user!
%      Default:
%      - Sparse arrays:   Struct containing the indices and values.
%      - Function handle: The reply of FUNCTIONS and the size and date of the
%        file.
%      - User defined and java objects: V.hashCode if existing, else: struct(V).
%
% NOTE:
%   For objects the function getByteStreamFromArray() might be exhaustive and
%   efficient, but unfortunately it is not documented.
%
% Tested: Matlab/64 7.8, 7.13, 8.6, 9.1, Win7/64
% Author: Jan Simon, Heidelberg, (C) 2016-2019 matlab.2010(a)n(MINUS)simon.de

% $JRev: R5n V:013 Sum:I9OCbLsTN/AG Date:17-Feb-2019 19:01:00 $
% $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
% $File: Tools\GLFile\GetMD5_helper.m $
% History:
% 001: 28-Jun-2015 19:19, Helper for GetMD5.
% 012: 27-Aug-2018 01:38, String class considered.
% 013: 09-Feb-2019 18:08, ismethod(class(V),) to support R2018b.

% Initialize: ==================================================================
% Do the work: =================================================================
% isa(V, 'class') is surprisingly expensive. Get the class as string once:
classV = class(V);

% The dimensions, number of dimensions and the name of the class is considered
% already in the MEX function!

if strcmp(classV, 'function_handle')  %#ok<*STISA>
   % Please adjust the subfunction ConvertFuncHandles to your needs.
   
   % The Matlab version influences the conversion by FUNCTIONS:
   % 1. The format of the struct replied FUNCTIONS is not fixed,
   % 2. The full path of toolbox function e.g. for @mean differ.
   S = functions(V);
   
   % Include modification file time and file size. Suggested by Aslak Grinsted:
   if ~isempty(S.file)
      d = dir(S.file);
      if ~isempty(d)
         S.filebytes = d.bytes;
         S.filedate  = d.datenum;
      end
   end
   
   % ALTERNATIVE: Use name and path. The <matlabroot> part of the toolbox
   % functions is replaced such that the hash for @mean does not depend on the
   % Matlab version.
   % Drawbacks: Anonymous functions, nested functions...
   % funcStruct = functions(FuncH);
   % funcfile   = strrep(funcStruct.file, matlabroot, '<MATLAB>');
   % S          = uint8([funcStruct.function, ' ', funcfile]);
   
   % Finally I'm afraid there is no unique method to get a hash for a function
   % handle. Please adjust this conversion to your needs.
   
elseif (isobject(V) || isjava(V)) && ismethod(classV, 'hashCode')  % ===========
   % Java or user-defined objects might have a hash function:
   S = char(V.hashCode);
   
elseif issparse(V)  % ==========================================================
   % Create struct with indices and non-zero values:
   [S.Index1, S.Index2, S.Value] = find(V);
   
elseif strcmp(classV, 'string')  % Since R2016b: ===============================
   % Convert all strings to UINT8 byte stream including the dimensions:
   classUint8 = uint8([117, 105, 110, 116, 49, 54]);
   C          = cell(1, numel(V));
   for iC = 1:numel(V)
      % Engine = CoreHash(uint16(Data{iS}), Engine);
      aString = uint16(V{iC});
      C{iC}   = [classUint8, ...
         typecast(uint64([ndims(aString), size(aString)]), 'uint8'), ...
         typecast(uint16(aString), 'uint8')];
   end
   S = [uint8([115, 116, 114, 105, 110, 103]), ...  % 'string'
      typecast(uint64([ndims(V), size(V)]), 'uint8'), ...
      cat(2, C{:})];
   
else  % Most likely this is a user-defined object: =============================
   try    % Perhaps a direct conversion is implemented:
      S = uint8(V);
      
      % Matt Raum had this excellent idea - unfortunately this function is
      % undocumented and might not be supported in te future:
      % S = getByteStreamFromArray(DataObj);

   catch ME  % Or perhaps this is better:
      fprintf(2, ['### %s: Convert object to struct as fallback.', char(10), ...
         '    %s\n'], ME.message);
      WarnS = warning('off', 'MATLAB:structOnObject');
      S     = struct(V);
      warning(WarnS);
   end
end

% return;
