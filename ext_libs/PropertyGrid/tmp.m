clear all

f = figure( ...
    'MenuBar', 'none', ...
    'Name', 'Property grid demo - Copyright 2010 Levente Hunyadi', ...
    'NumberTitle', 'off', ...
    'Toolbar', 'none');

h = PropertyGrid(f, ...
    'Position', [0.5 0 0.5 1]);

for pt=1:9
struc(pt).U=[0,1,0,0];
struc(pt).Im=[1000,1000,1000,1000];
end

h.Item = struc;        % bind object, discards any previously set properties

% update the type of a property assigned with type autodiscovery
%userproperties = PropertyGridField.GenerateFrom(struc);
%userproperties.FindByName('IntegerMatrix').Type = PropertyType('denserealdouble', 'matrix');
%disp(userproperties.FindByName('IntegerMatrix').Type);
h.Bind(struc, userproperties);

uiwait(f);

