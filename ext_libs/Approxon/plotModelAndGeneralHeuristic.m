%% Examples on the "General Heuristic" of 0.2V/mm as proposed in the literature
%
%  Andreas Husch, University of Luxembourg, 2019
%  andreas.husch(at)uni.lu
%
function pointHandle = plotModelAndGeneralHeuristic(ax)
load('astrom_table3_3v.mat'); %#ok<LOAD> % loads PW, D, T for 3V 
%[activation_model_3v, gof] =  fitModel(PW, D, T); %#ok<ASGLU>
%load('activation_model_3v.mat'); %#ok<LOAD>
load('logLinMod.mat'); % workaround for people without curve fitting toolbox

%% Hightlight 0.2V contour in 3D 
if(nargin < 1)
    fig = figure(); %#ok<NASGU>
    ax = axis();
else
    fig = gcf(); %#ok<NASGU>
end
axes(ax);
%h = plot(activation_model_3v, [PW, D], T ,'Xlim', [10 240], 'YLim', [1 8], 'Parent', ax);
x = linspace(10,240, 60);
y = linspace(1,8, 60);
[X,Y] = meshgrid(x,y);
Z = logLinMod(X,Y);

% plot model surface
hs = surface(X,Y,Z,'Parent', ax);%, 'EdgeColor', 'none');
hold(ax, 'on');

% plot general heuristic
[~,hc] = contour3(X,Y,Z, [0.2 0.2], 'r', 'LineWidth', 2, 'Parent', ax);

% plot data by astrom at all where the model was fitted to as smaller blue
hd = plot3(PW,D,T, 'o', 'MarkerSize', 5, ....
    'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'b', 'Parent', ax);

% plot large red point for current parameters
pointHandle = plot3(60,3.5,logLinMod(60,3.5), 'o', 'MarkerSize', 20, ....
    'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'r', 'Parent', ax); % initial choosen parameters
%vd = plot3(60,7.5,0.064, '*r', 'Parent', ax); % validation data point


title('Axonal E-Field Activation Threshold for (PW,D) at 3V', 'Parent', ax);

% Label axes
xlabel( 'PW [\mus]', 'Interpreter', 'Tex', 'Parent', ax);
ylabel( 'D [\mum]', 'Interpreter', 'Tex', 'Parent', ax);
zlabel( 'T [V/mm]', 'Interpreter', 'Tex', 'Parent', ax);
grid(ax, 'on');
view(ax, 157.2, 26.3 );
caxis(ax, [0.0369 2])


legend( [hs; hd ; pointHandle; hc], 'Model', 'Data (Aström et al.)', 'Choosen Parameters', 'Proposed General Heuristics of 0.2V/mm', 'Location', 'NorthEast', 'Interpreter', 'none' );

ax.Color = 'w';
ax.XTick = 30:30:270;

%% Altenative: Highlight "relevant" area:
%zlim([0 1])
%caxis([0.0369 0.5])

