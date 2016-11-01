load atlas_index

pastellize=3; % how much pastel.

P=jet;

PH=rgb2hsv(P);
PH(:,2)=PH(:,2)/pastellize;
P=hsv2rgb(PH);
J=jet;
J=rgb2hsv(J);
J(:,2)=J(:,2)/pastellize;
J=hsv2rgb(J);
atlases.colormap=P;

atlases.colors(1:55)=1:55; % use parula colors for others

% define structures 55 - 60 to have colors from JET

atlases.colormap(56+3,:)=J(20,:);
atlases.colors(56)=56+3;
atlases.colormap(57+3,:)=J(40,:);
atlases.colors(57)=57+3;
atlases.colormap(58+3,:)=J(60,:);
atlases.colors(58)=58+3;
atlases.colormap(59+3,:)=J(20,:);
atlases.colors(59)=59+3;
atlases.colormap(60+3,:)=J(30,:);
atlases.colors(60)=60+3;
atlases.colormap(61+3,:)=J(64,:);
atlases.colors(61)=61+3;



% also map to structures so we dont need to rebuild:
for ostr=1:60
   atlases.cdat{ostr,1}(:)=atlases.colors(ostr);
   atlases.cdat{ostr,2}(:)=atlases.colors(ostr);
end

save('atlas_index','atlases');