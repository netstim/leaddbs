function ea_computeGM(options,modes,finas,threshs,fs)
expfolder=[options.root,options.patientname,filesep,'connectomics',filesep,options.lc.general.parcellation,filesep];

for mode=1:length(modes)
    load([expfolder,finas{mode},'_CM.mat']);
    X=eval([modes{mode},'_CM']);
    if ~isnan(threshs{mode}) % apply threshold.
        X(X<threshs(mode))=0;
        X=logical(X);
    end
    clear([modes{mode},'_CM']);
    if options.lc.graph.degree_centrality
        disp(['Calculating degree centrality for ',finas{mode},' data...'])
        C=ea_deg(X);
        ea_export_cmeasure(C,'deg',finas{mode},options);
        disp('Done.');
    end
    if options.lc.graph.eigenvector_centrality
        disp(['Calculating eigenvector centrality for ',finas{mode},' data...'])
        C=ea_eig(X);
        ea_export_cmeasure(C,'eig',finas{mode},options);
        disp('Done.');
    end
    if options.lc.graph.nodal_efficiency
        disp(['Calculating nodal efficiency for ',finas{mode},' data...'])
        C=ea_eff(X);
        ea_export_cmeasure(C,'eff',finas{mode},options);
        disp('Done.');
    end
end


if options.lc.graph.struc_func_sim
    load([expfolder,'DTI_CM.mat']);
    Y=DTI_CM;
    for mode = find(fs==1)
        load([expfolder,finas{mode},'_CM.mat']);
        X=fMRI_CM;

        disp('Calculating structure-function similarity ...')
        C=ea_sfs(X,Y);

        ea_export_cmeasure(C,'sfs',['DTI_',finas{mode}],options);
        disp('Done.');
    end
end



function C=ea_deg(X)
C=nanmean(X);


function C=ea_eig(X)
X=X-ea_nanmin(X(:));

nnanix=~all(isnan(X));
smallX=X(nnanix,nnanix);

%X(isnan(X))=0;
[smallC,~]=eigs(sparse(smallX));
smallC=abs(smallC(:,1));
C=nan(size(X,1),1);
C(nnanix)=smallC;
%C=C/sum(C);


function C=ea_eff(X)
X(isnan(X))=0;

% based on BCT (c) M. Rubinov (see below)

if islogical(X)
    % based on
    %EFFICIENCY_BIN
    %   Mika Rubinov, U Cambridge
    %   Jonathan Clayden, UCL
    %   2008-2013

    % Modification history:
    % 2008: Original (MR)
    % 2013: Bug fix, enforce zero distance for self-connections (JC)
    % 2013: Local efficiency generalized to directed networks

    n=length(X);                                %number of nodes
    X(1:n+1:end)=0;                             %clear diagonal
    X=double(X~=0);                             %enforce double precision

    C=zeros(n,1);
    for u=1:n
        V=find(X(u,:)|X(:,u).');            %neighbors
        sa=X(u,V)+X(V,u).';                 %symmetrized adjacency vector
        e=distance_inv(X(V,V));             %inverse distance matrix
        se=e+e.';                           %symmetrized inverse distance matrix
        numer=sum(sum((sa.'*sa).*se))/2;    %numerator
        if numer~=0
            denom=sum(sa).^2 - sum(sa.^2);  %denominator
            C(u)=numer/denom;               %local efficiency
        end
    end
else

    % based on
    %EFFICIENCY_WEI
    %
    %   Mika Rubinov, U Cambridge, 2011-2012

    %Modification history
    % 2011: Original (based on efficiency.m and distance_wei.m)
    % 2013: Local efficiency generalized to directed networks

    n=length(X);                                    %number of nodes
    L = X;
    A = X~=0;
    ind = L~=0;
    L(ind) = 1./L(ind);                             %connection-length matrix

    C=zeros(n,1);
    for u=1:n
        V=find(A(u,:)|A(:,u).');                %neighbors
        sw=X(u,V).^(1/3)+X(V,u).^(1/3).';       %symmetrized weights vector
        e=distance_inv_wei(L(V,V));             %inverse distance matrix
        se=e.^(1/3)+e.'.^(1/3);                 %symmetrized inverse distance matrix
        numer=(sum(sum((sw.'*sw).*se)))/2;      %numerator
        if numer~=0
            sa=A(u,V)+A(V,u).';                 %symmetrized adjacency vector
            denom=sum(sa).^2 - sum(sa.^2);      %denominator
            C(u)=numer/denom;                   %local efficiency
        end
    end
end


function C=ea_sfs(X,Y)
X=atanh(X); % Fisher transform fucntional CM.
C=nan(size(X,1),1);
for n=1:size(X,2)
    % mask nan and inf values for both variables equally.
    nanix = isnan(X(:,n)) + isnan(Y(:,n)) + isinf(X(:,n)) + isinf(Y(:,n));
    if ~isempty(X(~nanix,n)) && ~isempty(Y(~nanix,n))
        C(n) = corr(X(~nanix,n),Y(~nanix,n));
    end
end


function ea_export_cmeasure(C,exstr,mode,options)
V=spm_vol([ea_space(options,'labeling'),options.lc.general.parcellation,'.nii']);
X=spm_read_vols(V);
X=round(X);
aID = fopen([ea_space(options,'labeling'),options.lc.general.parcellation,'.txt']);
atlas_lgnd=textscan(aID,'%d %s');
d=length(atlas_lgnd{1}); % how many ROI.
Y=X;
Y(:)=nan;
for node=1:d
    Y(X==node)=C(node);
end

expfolder=[options.root,options.patientname,filesep,'connectomics',filesep,options.lc.general.parcellation,filesep,'graph',filesep];
if ~exist(expfolder,'dir')
    mkdir(expfolder);
end

V.fname=[expfolder,exstr,'_',mode,'.nii'];
V.dt(1) = 64;
spm_write_vol(V,Y);

% also write out .MAT file with values:
save([expfolder,exstr,'_',mode,'.mat'],'C');


function D=distance_inv(A_)
l=1;                                        %path length
Lpath=A_;                                   %matrix of paths l
D=A_;                                       %distance matrix
n_=length(A_);

Idx=true;
while any(Idx(:))
    l=l+1;
    Lpath=Lpath*A_;
    Idx=(Lpath~=0)&(D==0);
    D(Idx)=l;
end

D(~D | eye(n_))=inf;                        %assign inf to disconnected nodes and to diagonal
D=1./D;                                     %invert distance


function D=distance_inv_wei(W_)
n_=length(W_);
D=inf(n_);                                      %distance matrix
D(1:n_+1:end)=0;

for u=1:n_
    S=true(1,n_);                               %distance permanence (true is temporary)
    W1_=W_;
    V=u;
    while 1
        S(V)=0;                                 %distance u->V is now permanent
        W1_(:,V)=0;                             %no in-edges as already shortest
        for v=V
            T=find(W1_(v,:));                   %neighbours of shortest nodes
            D(u,T)=min([D(u,T);D(u,v)+W1_(v,T)]);%smallest of old/new path lengths
        end

        minD=min(D(u,S));
        if isempty(minD)||isinf(minD),          %isempty: all nodes reached;
            break,                              %isinf: some nodes cannot be reached
        end;

        V=find(D(u,:)==minD);
    end
end

D=1./D;                                         %invert distance
D(1:n_+1:end)=0;
