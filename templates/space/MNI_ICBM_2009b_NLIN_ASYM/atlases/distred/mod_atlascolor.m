load atlas_index

pastellize=1; % how much pastel.

P=parula;

PH=rgb2hsv(P);
PH(:,2)=PH(:,2)/pastellize;
P=hsv2rgb(PH);
J=jet;
atlases.colormap=P;

atlases.colors(1:54)=1:54; % use parula colors for others

% define structures 55 - 60 to have colors from JET

atlases.colormap(55+4,:)=J(20,:);
atlases.colors(55)=55+4;
atlases.colormap(56+4,:)=J(40,:);
atlases.colors(56)=56+4;
atlases.colormap(57+4,:)=J(60,:);
atlases.colors(57)=57+4;
atlases.colormap(58+4,:)=J(20,:);
atlases.colors(58)=58+4;
atlases.colormap(59+4,:)=J(30,:);
atlases.colors(59)=59+4;
atlases.colormap(60+4,:)=J(64,:);
atlases.colors(60)=60+4;



% also map to structures so we dont need to rebuild:
for ostr=1:60
   atlases.cdat{ostr,1}(:)=atlases.colors(ostr);
   atlases.cdat{ostr,2}(:)=atlases.colors(ostr);
end

save('atlas_index','atlases');