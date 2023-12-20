function [tmp] = ea_diode_respecifyslices(tmp,ct,cscale,electrode,activecontact)
if ct.mat(1,1) < 0
    for k = 1:length(tmp)
        tmp{k} = flip(tmp{k},1);
    end
end
if ct.mat(2,2) < 0    
    for k = 1:length(tmp)
        tmp{k} = flip(tmp{k},2);
    end
end

for k = 1:length(tmp)
    sp(k) = subplot(numel(tmp),1,numel(tmp)-(k-1));
    imagesc(tmp{k}')
    view(-180,90)
    hold on
    colormap(gray)
    caxis manual
    caxis(cscale)
    axis equal
    axis off
%     title('Center -1')
end

%% show lead
ax_elec = axes('Position',[0.125 0.125 0.1 0.75]);
hold on
for k = 1:length(electrode.insulation)
    patch('Faces',electrode.insulation(k).faces,'Vertices',electrode.insulation(k).vertices,'Edgecolor','none','Facecolor',[electrode.lead_color electrode.lead_color electrode.lead_color]);
end
for k = 1:length(electrode.contacts)
    patch('Faces',electrode.contacts(k).faces,'Vertices',electrode.contacts(k).vertices,'Edgecolor','none','Facecolor',[electrode.contact_color electrode.contact_color electrode.contact_color]);
end
for k = activecontact
    patch('Faces',electrode.contacts(k).faces,'Vertices',electrode.contacts(k).vertices,'Edgecolor','none','Facecolor','r');
end
view (90,0)
ylim([-1 1])
xlim([-1 1])
zlim([0 max(electrode.insulation(end-1).vertices(:,3))])
axis equal
axis off
end
