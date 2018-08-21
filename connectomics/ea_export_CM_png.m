function cm=ea_export_CM_png(varargin)

X=varargin{1};
name=varargin{2};
options=varargin{3};

cm=figure('color','w','NumberTitle','off','Name',name);
if nargin==4
    imagesc(X,varargin{4});
else
    imagesc(X);
end
try
caxis([ea_nanmin(X(:))/2,ea_nanmax(X(:))/2]);
end
aID = fopen([ea_space(options,'labeling'),options.lc.general.parcellation,'.txt']);
atlas_lgnd=textscan(aID,'%d %s');
if length(atlas_lgnd{2})<20
    set(gca,'YTickLabel',atlas_lgnd{2}','YTickMode','manual','YTick',[1:length(atlas_lgnd{2})])
end
axis square
colorbar
