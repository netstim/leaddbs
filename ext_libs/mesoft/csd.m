

function [n weights fod] = csd(data,bDirs,HRdirs, Strue,lmax)


SHscheme = genSH(double(bDirs'),lmax);
SHschemeOUT = genSH(HRdirs',lmax);
SHscheme300 = genSH(HRdirs',lmax);


shhres = SH2RH(amp2SH(Strue,SHscheme));

%shhres/shhres(1)

[n weights fod] = csdeconv (shhres, data, SHscheme,SHscheme300);



function [ locmax weights FOD] = csdeconv (R_RH, S, DW_scheme, HR_scheme, lambda, tau)

% function [ F_SH, num_it ] = csdeconv (R_RH, S, DW_scheme, HR_scheme, lambda, tau)
%
% Constrained spherical deconvolution:
%
% deconvolves the axially symmetric response function 'R_RH' (in RH
% coefficients) from the function defined by the measurements 'S' taken along
% the directions of 'DW_Scheme' to give the surface 'F_SH' (in SH
% coefficients), by constraining the corresponding amplitudes of 'F_SH'
% (calculated along the directions given by 'HR_scheme') to be non-negative.
% The optional parameters 'lambda' and 'tau' correspond to the regularisation
% weight (1 by default) and the threshold on the FOD amplitude used to identify
% negative lobes (0.1 by default).



if ~exist('lambda'), lambda = 1; end
if ~exist('tau'), tau = 0.1; end


% guess appropriate value of lmax (assuming no super-resolution):
lmax = lmax_for_SH (S);
if 1 | HR_scheme.lmax < lmax
  lmax = HR_scheme.lmax; 
end

% build forward spherical convolution matrix up to lmax:
m = [];
for l = 0:2:lmax
  m = [ m R_RH(l/2+1)*ones(1,2*l+1) ];
end
fconv = DW_scheme.sh.*m(ones(size(DW_scheme.sh,1),1),:);

% generate initial FOD estimate, truncated at lmax = 4:
F_SH = fconv(:,1:nSH_for_lmax(4))\S;
F_SH (end+1:size(HR_scheme.sh,2),1) = 0;

% set threshold on FOD amplitude used to identify 'negative' values:
threshold = tau*mean (HR_scheme.sh*F_SH);


% scale lambda to account for differences in the number of 
% DW directions and number of mapped directions;
lambda = lambda * size(fconv,1)*R_RH(1)/size(HR_scheme.sh,1);

% main iteration loop:
fconv(:,end+1:size(HR_scheme.sh,2)) = 0;
k = [];
for num_it = 1:5,
  A = HR_scheme.sh*F_SH;
  k2 = find (A < threshold);
  if size(k2,1) + size (fconv,1) < size(HR_scheme.sh,2)
    disp ('too few negative directions identified - failed to converge');
    break;
  end
  if size(k) == size(k2), if k == k2, break; end; end
  k = k2;

  M = [ fconv; lambda*HR_scheme.sh(k,:) ];
  F_SH = M\[ S; zeros(size(k,1),1) ];
end


maxdir = 4;
FOD = SH2amp(F_SH,HR_scheme);   
locmax= reshape(getLocalMax(FOD,HR_scheme.vert',maxdir),[3 maxdir]);
weights = sqrt(sum(locmax.^2));
locmax = locmax ./ repmat(weights,[3 1]); 


function lmax = lmax_for_SH (SH)

lmax = 2*(floor(sqrt(1+8*size(SH,1))-3)/4);

function S = c2s(X)
S(:,3) = sqrt(sum(X.^2, 2));
S(:,1) = acos(X(:,3)./S(:,3));
S(:,2) = atan2(X(:,2), X(:,1));


function scheme = genSH(N, lmax)

P = c2s(N);

scheme.el = P(:,1);
scheme.az = P(:,2);

scheme.sh = [];
scheme.lmax = lmax;

for l = 0:2:lmax
  scheme.sh = [ scheme.sh eval_SH(l, scheme.el, scheme.az)' ];
end

scheme.vert= s2c([ scheme.el scheme.az 1+0*scheme.az]);
scheme.mesh = convhulln(scheme.vert);

function s = eval_ALP(l, el)

  s = legendre(l, cos(el'));
  for m = 0:l
    s(m+1,:) = s(m+1,:).*sqrt((2*l+1)*factorial(l-m) / ((4*pi)*factorial(l+m)));
  end

  if l
    s = [ s(end:-1:2,:); s ];
  end

function RH = SH2RH (SH)

lmax = lmax_for_SH (SH);
D_SH = gen_delta (0, 0, lmax);
k = find (D_SH);
RH = SH(k)./D_SH(k);

function SH = gen_delta(el, az, lmax)
SH = [];

for l = 0:2:lmax
  s = eval_SH(l, el, az);
  SH = [ SH; s ];
end

function SH = amp2SH (S, scheme);

if ~isfield (scheme, 'shinv')
  scheme.shinv = pinv(scheme.sh);
end

sz = size(S);
SH = reshape(scheme.shinv * S(:,:),[size(scheme.shinv,1) sz(2:end)]);


function s = eval_SH(l, el, az)

s = ones(size(az,1),1);

if l > 0  
  s = [ sqrt(2)*sin(az*(l:-1:1)) s sqrt(2)*cos(az*(1:l)) ];
end

s = eval_ALP(l, el).*s';

function n = nSH_for_lmax (lmax)

% function n = nSH_for_lmax (lmax)
%
% returns the number of even SH coefficients for a 
% given harmonic order 'lmax'

n = (lmax+1)*(lmax+2)/2;



function dir = getLocalMax(sigsym,bDir,nummax)

if evalin('base','exist(''getlocmaxstruct'')')
    gls = evalin('base','getlocmaxstruct');
    if size(bDir,2) ~= size(bDir,2),
        clear gls;
    elseif any(abs(bDir(:)-gls.bDir(:))>eps),
        clear gls;
    end;
end;
    


if not(exist('gls')),
    ncliques = computeNeighbors(double(bDir'));
    locFitPS = zeros(6,size(ncliques.neighborsarray,1)+1,size(bDir,2));
    bDirt = bDir';
    for k = 1:length(ncliques.neighbors),
         Ns = [ncliques.neighborsarray(:,k); k];
         Ms = [bDirt(Ns,1).^2 bDirt(Ns,2).^2 bDirt(Ns,3).^2 2*bDirt(Ns,1).*bDirt(Ns,2) 2*bDirt(Ns,1).*bDirt(Ns,3) 2*bDirt(Ns,2).*bDirt(Ns,3)];
         locFitPS(:,:,k) = pinv(Ms);
    end;
    gls.ncliques = ncliques;
    gls.locFitPS = locFitPS;
    gls.bDir = bDir;
    assignin('base','getlocmaxstruct',gls);
end;
    
usesym = 1;
dir = getLocalMaxC(double(sigsym),double(bDir),double(gls.ncliques.neighborsarray),double(gls.locFitPS),double([nummax usesym]));




function S = SH2amp (SH, scheme)

% function S = SH2amp (SH, scheme)
%
% Maps SH coefficients 'SH' to amplitudes along directions
% in 'scheme'

S = scheme.sh(:, 1:size(SH,1))*SH;


