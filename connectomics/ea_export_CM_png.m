function cm=ea_export_CM_png(X,name,options)


cm=figure('color','w','NumberTitle','off','Name',name);
imagesc(X);
aID = fopen([options.earoot,'templates',filesep,'labeling',filesep,options.lc.general.parcellation,'.txt']);
atlas_lgnd=textscan(aID,'%d %s');
if length(atlas_lgnd{2})<20
set(gca,'YTickLabel',atlas_lgnd{2}','YTickMode','manual','YTick',[1:length(atlas_lgnd{2})])
end
axis square
colorbar
