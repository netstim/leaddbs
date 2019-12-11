%% Examples on the "General Heuristic" of 0.2V/mm as proposed in the literature
%
%  Andreas Husch, University of Luxembourg, 2019
%  andreas.husch(at)uni.lu
%
function pointHandle = plotModelAndGeneralHeuristic(ax)
load('astrom_table3_3v.mat'); %#ok<LOAD> % loads PW, D, T for 3V 
%[activation_model_3v, gof] =  fitModel(PW, D, T); %#ok<ASGLU>
load('activation_model_3v.mat'); %#ok<LOAD>

%% Hightlight 0.2V contour in 3D 
if(nargin < 1)
    fig = figure(); %#ok<NASGU>
    ax = axis();
else
    fig = gcf(); %#ok<NASGU>
end
axes(ax);
h = plot(activation_model_3v, [PW, D], T ,'Xlim', [10 240], 'YLim', [1 8], 'Parent', ax);
hold(ax, 'on');
pointHandle = plot3(60,3.5,activation_model_3v(60,3.5), 'o', 'MarkerSize', 20, ....
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

x = linspace(10,240);
y = linspace(1,8);
[X,Y] = meshgrid(x,y);
Z = activation_model_3v(X,Y);
[~,hc] = contour3(X,Y,Z, [0.2 0.2], 'r', 'LineWidth', 2, 'Parent', ax);

legend( [h; pointHandle; hc], 'Model', 'Data (Aström et al.)', 'Choosen Parameters', 'Proposed General Heuristics of 0.2V/mm', 'Location', 'NorthEast', 'Interpreter', 'none' );

ax.Color = 'w';
ax.XTick = 30:30:270;

%% Altenative: Highlight "relevant" area:
%zlim([0 1])
%caxis([0.0369 0.5])

