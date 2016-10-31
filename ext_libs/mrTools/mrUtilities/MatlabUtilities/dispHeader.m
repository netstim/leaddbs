% dispHeader.m
%
%        $Id:$ 
%      usage: dispHeader(header,<len=40>,<c='='>)
%         by: justin gardner
%       date: 04/22/11
%    purpose: displays a line of 60 characters as a header line, for example
%
%             >> dispHeader('header text')
%             ============= header text ==============
%
%             You can change the text length and or the separator character:
%             >> dispHeader('header text',20,'+')
%             +++ header text ++++
%
function retval = dispHeader(header,len,c)

% check arguments
if ~any(nargin == [0 1 2 3])
  help dispHeader
  return
end

% default header is just a full line
if nargin < 1, header = '';end

% default length
if (nargin < 2) || isempty(len),len = 60;end

% default separator character
if (nargin < 3) || isempty(c),c = '=';end

% get length of texgt
headerLen = length(header);

% if it is longer than the desired header length, then
% display two lines of separators one above and below the header
if (headerLen+2) >= len
  disp(repmat(c,1,len));
  disp(header)
  disp(repmat(c,1,len));
elseif headerLen == 0
  % if the header is empty, just display a full line
  disp(repmat(c,1,len));
else
  % otherwise put header inside separator characters
  fillerLen = ((len-(headerLen+2))/2);
  
  % display the first part
  disp(sprintf('%s %s %s',repmat(c,1,floor(fillerLen)),header,repmat(c,1,ceil(fillerLen))));
end
