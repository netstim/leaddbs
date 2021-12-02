function outS=ea_conflateS_lr(S)

outS=S(1); % rh
for source=1:4
outS.(['Ls',num2str(source)])=S(2).(['Ls',num2str(source)]);
end

outS.amplitude{2}=S(2).amplitude{2};
outS.activecontacts{2}=S(2).activecontacts{2};

