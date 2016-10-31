% gethdr.m
%
%        $Id$	
%      usage: gethdr(d,x,y,s)
%         by: justin gardner
%       date: 08/21/03
%    purpose: retrieve the hdrs in a matrix
%             form from the data structure
%
function [hdr,time,hdrste] = gethdr(d,x,y,s)

if ~any(nargin == [2 4])
  help gethdr;
  return
end

% if only two arguments, it must be a 3 tuple
if nargin == 2
  if length(x) == 3
    y = x(2);
    s = x(3);
    x = x(1);
  else
    help gethdr;
    return
  end
end

hdr = shiftdim(d.ehdr(x,y,s,:,:), 3);
hdrlen = d.hdrlen;

% check if this is a glm analysis
% then we scale the hrd 
if isfield(d, 'hrf')
%    hdr = hdr*(d.hrf'-mean(d.hrf));
    hdr = hdr*d.hrf';
    hdrlen = size(d.hrf,1);
end

% return the time as well if asked for
if nargout >= 2
  time = d.tr/2:d.tr:(hdrlen*d.tr);
end

% return the standard error
if nargout >=3
  hdrste = shiftdim(d.ehdrste(x,y,s,:,:), 3);
  % not sure about this part
  if isfield(d, 'hrf')
     %hdrste = sqrt(hdrste.^2*abs(d.hrf)');
     %JB: this is the right way of doing it:
     hdrste = sqrt(hdrste.^2*(d.hrf.^2)'); %the variance (squared stddev) of a random variable multiplied by a constant 
     %is equal to the variance multiplied by the square of the constant. 
  end

end
% if we only have one response, then we will have
% to take the transpose to make sure that we have
% row array
% if size(d.ehdr,4) == 1 & nargout >=3
%  hdr = hdr';
%  hdrste = hdrste';
% end

