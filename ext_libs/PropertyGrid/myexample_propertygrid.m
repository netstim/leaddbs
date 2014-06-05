% Demonstrates how to use the property pane.
%
% See also: PropertyGrid

% Copyright 2010 Levente Hunyadi
function myexample_propertygrid
patients={'examplepatient1','examplepatient2'};
for p=1:length(patients)
Uproperties(p) =    PropertyGridField('double', pi, ...
        'Category', 'Voltage', ...
        'DisplayName', patients{p}, ...
        'Description', 'Double Matrix.');

Iproperties(p) =    PropertyGridField('double', pi, ...
        'Category', 'Impedance', ...
        'DisplayName', patients{p}, ...
        'Description', 'Double Matrix.');
end




% arrange flat list into a hierarchy based on qualified names
Uproperties = Uproperties.GetHierarchy();
Iproperties = Iproperties.GetHierarchy();

% create figure
f = figure( ...
    'MenuBar', 'none', ...
    'Name', 'Property grid demo - Copyright 2010 Levente Hunyadi', ...
    'NumberTitle', 'off', ...
    'Toolbar', 'none');

% procedural usage
gU = PropertyGrid(f, ...            % add property pane to figure
    'Properties', Uproperties,'Position', [0 0 0.5 1]);     % set properties explicitly;
gI = PropertyGrid(f, ...            % add property pane to figure
    'Properties', Iproperties,'Position', [0.5 0 0.5 1]);     % set properties explicitly;


% declarative usage, bind object to grid
obj = SampleObject;  % a value object

% update the type of a property assigned with type autodiscovery
%userproperties = PropertyGridField.GenerateFrom(obj);
%userproperties.FindByName('IntegerMatrix').Type = PropertyType('denserealdouble', 'matrix');
%disp(userproperties.FindByName('IntegerMatrix').Type);

% wait for figure to close
uiwait(f);

% display all properties and their values on screen
disp('Left-hand property grid');
disp(gU.GetPropertyValues());
keyboard
