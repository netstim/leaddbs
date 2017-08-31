function ea_error(msg, title, backtrace)
% The general error message function.
% 
% set backtrace to 'dbstack' if needed.
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ~exist('title', 'var') || isempty(title)
    title = '';
end

errordlg(msg, title);

if ~exist('backtrace', 'var') % simple mode
    error(msg);
else % suppress complete backtrace, only keep the first one
    fprintf('\n');
    errstruct.message = msg;
    errstruct.identifier = '';
    errstruct.stack.file = backtrace(1).file;
    errstruct.stack.name = backtrace(1).name;
    errstruct.stack.line = backtrace(1).line;
    error(errstruct);
end
