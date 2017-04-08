function c = ea_colormap(map)
%
% Small function to return colormap as a matrix without affecting your scene
% __________________________________________________________________________________
% Copyright (C) 2017 University of Pittsburgh, Brain Modulation Lab
% Ari Kappel

figH = figure('visible','off');

c = colormap(map);

close(figH)

