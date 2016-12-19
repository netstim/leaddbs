% rotate DE_scheme in patient coordinate system

% created by: Susanne Schnell (11.09.2008) on a PC
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [varargout] = rotate_de_scheme(DE_scheme, edges)

% check input
if nargin == 0
    errStr = 'Not enough input arguments given');
    DE_scheme = [];
elseif nargin == 1
    errStr = 'Not enough input arguments given');
    DE_scheme = [];
elseif nargin == 2
    DE_scheme = squeeze(DE_scheme);
    edges = squeeze(edges);
    errStr = [];
end

% rotate DE gradients
for m = 1 : size(DE_scheme,1)
    if isempty(edges)
        errStr = 'Error in rotate_de_scheme: parameter "edges" missing.';
        return;
    end
    [p_org,errStr]= hm_analyse(edges);
    p_rot= [0 0 0 p_org(4:6) sign(p_org(7:9)) 0 0 0];
    rotHM= hm_create(p_rot);
    rotM= rotHM(1:3, 1:3)';
    DE_scheme(m,:) = rotM*DE_scheme(m,:)';
    DE_scheme = DE_scheme';
    errStr = [];
end


% check output
if nargout == 1
    varargout{1} = DE_scheme; 
elseif nargout == 2
    varargout{1} = DE_scheme;
    varargout{2} = errStr;
elseif nargout > 2
    varargout{1} = DE_scheme;
    varargout{2} = errStr;
    for m = 3 : nargout
        varargout{m} = [];
    end
else
    mes = 'Error in rotate_de_scheme: not enough output arguments given';
    error(mes);
end
