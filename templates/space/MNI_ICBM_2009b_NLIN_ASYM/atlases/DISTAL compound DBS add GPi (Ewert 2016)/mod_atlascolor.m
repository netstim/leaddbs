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

% define structures 31-33 and 62 - 64 to have colors from JET

atlases.colormap(32,:)=J(20,:); % STN assoc
atlases.colors(32)=32;
atlases.colormap(33,:)=J(40,:); % STN limb
atlases.colors(33)=33;
atlases.colormap(34,:)=J(60,:); % STN mot
atlases.colors(34)=34;
atlases.colormap(62,:)=J(20,:); % GPe
atlases.colors(62)=62;
atlases.colormap(63,:)=J(30,:); % GPi
atlases.colors(63)=63;
%atlases.colormap(63,:)=J(30,:); % GPi
atlases.colors(9)=40;
atlases.colormap(64,:)=J(64,:); % RN
atlases.colors(64)=64;


% also map to structures so we dont need to rebuild:
for ostr=1:64
   atlases.cdat{ostr,1}(:)=atlases.colors(ostr);
   atlases.cdat{ostr,2}(:)=atlases.colors(ostr);
end

save('atlas_index','atlases');