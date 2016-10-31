% hrfDiffGamma.m
%
%      usage: hrfDiffGamma('params'), hrfDiffGamma(tr, params)
%         by: farshad moradi
%       date: 14/06/07
%    purpose: returns a canonical hrf
%
function [hrf] = hrfDiffGamma(tr, params)

if ~any(nargin == [1 2])
  help hrfDiffGamma
  return
end

if tr=='params'
    hrf = {...
        {'description', 'hrfDiffGamma', 'comment describing the hdr model'},...
        {'x', 6, 'hrf = gampdf(...,x,1)-gampdf(...,y,1)/z'},...
        {'y', 16, 'hrf = gampdf(...,x,1)-gampdf(...,y,1)/z'},...
        {'z', 6, 'hrf = gampdf(...,x,1)-gampdf(...,y,1)/z'},...
        {'stimDur', 0.01, 'duration of stimulation/event (seconds, min=0.01s). a boxcar function that is convolved with hrf'},...
        {'incDeriv',0,'type=checkbox','include derivative of the hrf in the model?'},...
    };
    return
end

tmax = max(params.y*3, 20);

if isfield(params, 'tmax')
    tmax = params.tmax;
end

shift = 0;
if isfield(params, 'shift')
    shift = params.shift;
end

dt = 0.01;

t = 0:dt:tmax;
HRF = gampdf(t, params.x, 1) - gampdf(t, params.y, 1)/params.z;
HRF = convn(HRF, ones(1, max(1, ceil(params.stimDur/dt))) );

if shift<0
    HRF = [zeros(1, ceil(-shift/dt)), HRF];
elseif shift>0
    HRF = HRF( ceil(shift/dt):end );
end
    
% subsample hrf
t = [0:length(HRF)-1]*dt;
h_intp = interp1(t, HRF, tr/2:tr:max(t));

% remove mean
h_intp = h_intp - mean(h_intp);
% normalize
h_intp = h_intp / norm(h_intp'); 


if params.incDeriv
    
    % take the derivative
    HRFD = [diff(HRF), 0];
    
    % subsample hrf derivative
    hd_intp = interp1(t, HRFD, tr/2:tr:max(t));
    
    % remove mean
    hd_intp = hd_intp - mean(hd_intp);
    % orthogonalize
    hd_intp = hd_intp - h_intp*(hd_intp/h_intp);
    % normalize
    hd_intp = hd_intp / norm(hd_intp');

    % return as column vectors, zero at time zero
    hrf = [h_intp; hd_intp]';
else
    % return as column vector, zero at time zero
    hrf = h_intp';
end
