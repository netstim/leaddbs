function mrCloseDlg(h)
%
% function mrCloseDlg(h)
%
% Closes a dialog box/message box opened with mrMsgBox or mrWaitBar.
% h is a handle to the window.
% Checks to make sure that h exists before trying to close it.
%
% djh, 7/2005

if ishandle(h)
	close(h);
elseif isfield(h,'disppercent')
  disppercent(inf);
end

