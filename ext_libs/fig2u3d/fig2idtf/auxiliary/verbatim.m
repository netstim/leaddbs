function str = verbatim
%VERBATIM  Get the text that appears in the next comment block.
%  Returns the text of the first comment block following the call.  The
%  block comment delimiters, %{ and %}, must appear on lines by themselves
%  (optionally preceded by white space).
%
%  If you want a final end-of-line character then leave a blank line before
%  the %}.
%
%  If both comment delimiters are preceded by the same white space (same
%  combination of spaces and tabs) then that white space will be deleted
%  (if possible) from the beginning of each line of the commented text.
%  This is so the whole block can be indented.
%
%  Example,
%
%      str = verbatim;
%          %{
%          This is the text
%          that will be returned by verbatim.
%          %}
%
%  VERBATIM can only be used in an m-file.


% Version: 1.1, 22 March 2009
% Author:  Douglas M. Schwarz
% Email:   dmschwarz=ieee*org, dmschwarz=urgrad*rochester*edu
% Real_email = regexprep(Email,{'=','*'},{'@','.'})


% Get the function call stack.
[dbs,thisWorkspace] = dbstack('-completenames');
assert(length(dbs) > 1,'VERBATIM must be called from an M-file.')
dbs = dbs(thisWorkspace + 1);
lines = repmat({''},1,100);

% Open the file.
fid = fopen(dbs.file);

try
	% Skip lines up to the current line in the calling function.
	textscan(fid,'%*s',0,'HeaderLines',dbs.line,'Delimiter','');

	% Read lines until one begins with '%{'.
	line = '';
	while ~isequal(line,-1) && isempty(regexp(line,'^\s*%{','once'))
		line = fgetl(fid);
	end

	% Read and save lines until one begins with '%}'.  Leave first cell
	% empty.
	k = 1;
	while true
		k = k + 1;
		lines{k} = fgetl(fid);
		if isequal(lines{k},-1) || ...
				~isempty(regexp(lines{k},'^\s*%}','once'))
			break
		end
	end

	% Close the file.
	fclose(fid);

catch err
	% Close the file and rethrow the error.
	fclose(fid);
	rethrow(err)
end

% If white space preceeding '%{' and '%}' is the same then delete it from
% the beginning of each line of the commented text so you can have indented
% code.
white_space1 = regexp(line,'^\s*','match','once');
white_space2 = regexp(lines{k},'^\s*','match','once');
if strcmp(white_space1,white_space2)
	lines = regexprep(lines,['^',white_space1],'');
end

% Construct the output string.
str = [sprintf('%s\n',lines{2:k-2}),lines{k-1}];
