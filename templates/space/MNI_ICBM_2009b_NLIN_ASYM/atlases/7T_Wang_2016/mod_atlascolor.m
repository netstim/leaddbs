load atlas_index

pastellize=1.5; % how much pastel.

P=jet;

PH=rgb2hsv(P);
PH(:,2)=PH(:,2)/pastellize;
P=hsv2rgb(PH);
J=jet;
J=rgb2hsv(J);
J(:,2)=J(:,2)/pastellize;
J=hsv2rgb(J);
atlases.colormap=P;

atlases.colors(1:end)=round(linspace(20,50,length(atlases.names))); % use parula colors for others

atlases.colormap(atlases.colors(3),:)=[0.8,0.2,0.2]; % RN


% also map to structures so we dont need to rebuild:
for ostr=1:length(atlases.names)
   atlases.cdat{ostr,1}(:)=atlases.colors(ostr);
   atlases.cdat{ostr,2}(:)=atlases.colors(ostr);
end

save('atlas_index','atlases');