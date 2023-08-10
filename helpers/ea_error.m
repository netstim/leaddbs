function ea_error(msg, options)
% Wrapper function for error which can also show errordlg and support
% simplied stack information.

arguments
    msg {mustBeTextScalar}
    options.title {mustBeTextScalar} = 'Error'
    options.showdlg {mustBeNumericOrLogical} = true
    options.simpleStack {mustBeNumericOrLogical} = false
end

msg = sprintf(msg);

if options.showdlg
    errordlg(msg, options.title);
end

% Construct errorStruct
fprintf('\n');
errorStruct.message = msg;
errorStruct.identifier = '';
errorStruct.stack = dbstack('-completenames');

if length(errorStruct.stack) > 1 % Call from other function
    % Override the first stack (which is ea_error itself)
    errorStruct.stack(1) = errorStruct.stack(2);
    
    % Only keep the top stacks if simpleStack is set to true
    if options.simpleStack
        errorStruct.stack = errorStruct.stack(1:2);
    end

    error(errorStruct);
else % Call directly
    ea_cprintf('CmdWinErrors', '%s\n\n', msg);
end
