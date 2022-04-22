function varargout=ea_genvat_kuncel(varargin)
% This function generates a volume of activated tissue around for each
% electrode.
% Usage: VAT=ea_genvat(coords_mm,stimparams,options).
% ? stimparams is a struct variable with fields U (8*1 with voltage
% entries) and Im (8*1 with Impedance measurements).
%
% This function only touches the .VAT entry of stimparams struct of the
% given side.

if nargin>=5
    coords=varargin{1};
    S=varargin{2};
    side=varargin{3};
    options=varargin{4};
    stimname=varargin{5};
elseif nargin==1
    if ischar(varargin{1}) % return name of method.
        varargout{1}='Kuncel 2008';
        varargout{2} = false; % Doesn't support directed lead
        return
    end
end

switch side
    case 1
        sidec='R';
        cnts={'k0','k1','k2','k3','k4','k5','k6','k7'};
    case 2
        sidec='L';
        cnts={'k8','k9','k10','k11','k12','k13','k14','k15'};
end

[xx,yy,zz]=psphere(100);

try
    if isnan(S)
        clear stimparams
    end
end

radius=repmat(1.5,options.elspec.numel,1); % some default setting.

%try % if stimparams are set.
if ~isfield(S, 'sources')
    S.sources=1:4;
end

for source=S.sources

    stimsource=S.([sidec,'s',num2str(source)]);

    for cnt=1:length(cnts)
        U(cnt)=stimsource.(cnts{cnt}).perc;
    end
    Acnt=find(U>0);
    if length(Acnt)>1
        ea_error('In the Maedler and Kuncel models, only one active contact can be selected in each source');
    end
    U=stimsource.amp;

    radius(source)=kuncel08_eq1(U);
    volume(source)=(4/3)*pi*radius(source)^3;

    VAT{source}=[xx*radius(source)+coords{side}(Acnt,1);...
        yy*radius(source)+coords{side}(Acnt,2);...
        zz*radius(source)+coords{side}(Acnt,3)]';
    K{source}=convhulln(VAT{source}+randn(size(VAT{source}))*0.000001); % create triangulation.

    for dim=1:3
        ivx(source,dim,:)=[min(VAT{source}(:,dim)),max(VAT{source}(:,dim))];
    end
end

aivx=zeros(3,2);
aivx(:,1)=min(ivx(:,:,1));
aivx(:,2)=max(ivx(:,:,2));

voxspace=zeros(abs(floor(aivx(:,1))-ceil(aivx(:,2)))'*5);

[xxv,yyv,zzv]=ind2sub(size(voxspace),1:numel(voxspace));
XYZv=[xxv;yyv;zzv;ones(1,length(xxv))];
gvmm{1}=linspace(floor(aivx(1,1)),ceil(aivx(1,2)),size(voxspace,1));
gvmm{2}=linspace(floor(aivx(2,1)),ceil(aivx(2,2)),size(voxspace,2));
gvmm{3}=linspace(floor(aivx(3,1)),ceil(aivx(3,2)),size(voxspace,3));
XYZmm=[gvmm{1}(XYZv(1,:));...
    gvmm{2}(XYZv(2,:));...
    gvmm{3}(XYZv(3,:));...
    ones(1,length(XYZv))];
mat=mldivide(XYZv',XYZmm')';

for source=S.sources
    in=ea_intriangulation(VAT{source},K{source},XYZmm(1:3,:)');
    voxspace(sub2ind(size(voxspace),XYZv(1,in),XYZv(2,in),XYZv(3,in)))=1;
end

% write nifti of VAT
[~, ~, endian] = computer;
switch endian
    case 'L'
        endian = 0;
    case 'B'
        endian = 1;
end

Vvat.mat=mat;
%voxspace=permute(voxspace,[2,1,3]);
Vvat.dim=size(voxspace);
Vvat.dt = [4, endian];
Vvat.n=[1 1];
Vvat.descrip='lead dbs - vat';

stimDir = fullfile(options.subj.stimDir, ea_nt(options), stimname);
ea_mkdir(stimDir);
filePrefix = ['sub-', options.subj.subjId, '_sim-'];

if side == 1
    Vvat.fname = [stimDir, filesep, filePrefix, 'binary_model-kuncel_hemi-R.nii'];
elseif side == 2
    Vvat.fname = [stimDir, filesep, filePrefix, 'binary_model-kuncel_hemi-L.nii'];
end

ea_savestimulation(S,options);
spm_write_vol(Vvat,voxspace);

varargout{1}=VAT;
varargout{2}=volume;
varargout{3}=radius;


function r=kuncel08_eq1(U)
% This function calculates the  radius of Volume of Activated Tissue for
% stimulation settings U (Kuncel 2008).

r=0; %
if U %(U>0)
    k=0.22;
    Uo=0.1;
    r=sqrt((U-Uo)/k);
end


function [x,y,z,avgr] = psphere(n,rad)
% [x,y,z,r] = psphere(N)
%
% Distributes N points "equally" about a unit sphere.
%
% N The number of points to distribute
% x,y,z Each is 1 x N vector
% r The smallest linear distance between two neighboring
% points. If the function is run several times for the
% same N, r should not change by more than the convergence
% criteria, which is +-0.01 on a unit sphere.
%
%
% Distributes N points about a unit sphere so that the straight line
% distance between neighboring points is roughly the same. The
% actual criteria for stopping is slightly different. The difference
% between a given point and every other point is stored in a matrix.
% The iterations stop once the maximum difference between any element
% in successive distance matrices is less than 0.01. An absolute
% criteria was chosen due to self-distances being 0, and programming
% around this for relative convergence seemed too much work for too
% little reward.
%
% The algorithm first generates N random points. Then a repulsive
% force vector, based on 1/r^2, is calculated for each point due to
% all of the other points. The resultant force vector for each point
% is normalized, and then each point is displaced a distance S = 1
% in the direction of the force. Any value of S between 0.0 and 1.
% 0 seems to work with most values between 0.2 and 1 taking an average
% of 20 to 30 iterations to converge. If S is too high, too much
% "energy" is being added to the system and it won't converge. If S is
% too low, the points won't be evenly distributed even though the
% convergence criteria is met. Generally speaking, the larger N is
% the larger S can be without causing instabilities. After
% displacement, the point is projected back down onto the unit sphere.
% When the system nears convergence, the displacement vector for a
% given point is nearly in the same direction as the radius vector for
% that point due to the points being equally distributed. A good check
% to make sure the code is working is to specify N = 4 and then check
% that the resulting points form a regular tetrahedron (or whatever
% it's called). How you would do this I don't know (check angles
% maybe). That's why I supplied the demo option so you could look at
% it in progress.
%
% Jason Bowman
% jbowman90@hotmail.com
% Last Revised June 2000
%
% You may freely distribute this code. I only ask that you give credit
% where credit is due.

%Since rand produces number from 0 to 1, subtract off -0.5 so that
%the points are centered about (0,0,0).
x = rand(1,n) - 0.5;
y = rand(1,n) - 0.5;
z = rand(1,n) - 0.5;

%Make the matrix R matrices for comparison.
rm_new = ones(n);
rm_old = zeros(n);

%Scale the coordinates so that their distance from the origin is 1.
r = sqrt(x.^2 + y.^2 + z.^2);

x = x./r;
y = y./r;
z = z./r;

not_done = 1;

s = 1;

%Turns off the divide by 0 warning
warning off

while not_done
    for i = 1:n
        %Calculate the i,j,k vectors for the direction of the repulsive forces.
        ii = x(i) - x;
        jj = y(i) - y;
        kk = z(i) - z;

        rm_new(i,:) = sqrt(ii.^2 + jj.^2 + kk.^2);

        ii = ii./rm_new(i,:);
        jj = jj./rm_new(i,:);
        kk = kk./rm_new(i,:);

        %Take care of the self terms.
        ii(i) = 0;
        jj(i) = 0;
        kk(i) = 0;

        %Use a 1/r^2 repulsive force, but add 0.01 to the denominator to
        %avoid a 0 * Inf below. The self term automatically disappears since
        %the ii,jj,kk vectors were set to zero for self terms.
        f = 1./(0.01 + rm_new(i,:).^2);

        %Sum the forces.
        fi = sum(f.*ii);
        fj = sum(f.*jj);
        fk = sum(f.*kk);

        %Find magnitude
        fn = sqrt(fi.^2 + fj.^2 + fk.^2);

        %Find the unit direction of repulsion.
        fi = fi/fn;
        fj = fj/fn;
        fk = fk/fn;

        %Step a distance s in the direciton of repulsion
        x(i) = x(i) + s.*fi;
        y(i) = y(i) + s.*fj;
        z(i) = z(i) + s.*fk;

        %Scale the coordinates back down to the unit sphere.
        r = sqrt(x(i).^2 + y(i).^2 + z(i).^2);

        x(i) = x(i)/r;
        y(i) = y(i)/r;
        z(i) = z(i)/r;
    end

    %Check convergence
    diff = abs(rm_new - rm_old);

    not_done = any(diff(:) > 0.01);

    rm_old = rm_new;
end %while

%Find the smallest distance between neighboring points. To do this
%exclude the self terms which are 0.
tmp = rm_new(:);

indices = find(tmp~=0);

avgr = min(tmp(indices));

%Turn back on the default warning state.
warning backtrace
