function c = ea_colorlover(n,plot)
% Version 2.0
% c = colorlover(palette_number,plot);
% c is a matrix with 5 colors in matlab color codes 
% e.g. plot(x,y,'color',c(1,:)) will plot in color 1 from the matrix
% if n is empty a random palette will be chosen
if ~exist('plot','var') 
    plot = 0;
end

if ~exist('n','var')
    n=19;
end
colors{19} = [213,62,79
244,109,67
253,174,97
254,224,139
255,255,191
230,245,152
171,221,164
102,194,165
50,136,189];
colors{1} = [zeros(1,3)+.9;zeros(1,3)+.7;zeros(1,3)+.5;zeros(1,3)+.3;zeros(1,3)+.1];
colors{2} = [255,120,105;255,188,84;244,255,98;167,218,255;255,36,191];
colors{3} = [0.0431,0.2078,0.3765;0.3843,0.6353,0.6588;0.6706,0.2235,0.1176;0.5,0.5,0.5;0.9,0.9,.9];
% Northwest Rain
colors{4} = [252,237,208;62,85,103;207,255,255;124,30,28;105,152,180];
% Ocean Five
colors{5} = [0,160,176;106,74,60;204,51,63;235,104,65;237,201,81];
% It's raining love
colors{6} = [163,169,72;237,185,46;248,89,49;206,24,54;0,153,137];
colors{7} =[2 157 175;227 37 82; 229 213 153; 255 194 25; 249 212 35; 121 183 180];
colors{8} = [249 212 35; 237 229 116; 225 245 196; 173 214 188; 121 183 180];
% flowers
colors{17} = [207,95,83;251,153,82;228,193,135;139,176,135;112,140,115];
% Thought Provoking
colors{18} = [236,208,120;217,91,67;192,41,66;84,36,55;83,119,122];
% Adrift in Dreams
colors{9} = [207,240,158 ; 168,219,168 ; 121,189,154 ; 59,134,134 ; 11,72,107];
% you are beautiful
colors{10} = [53,19,48 ; 66,66,84 ; 100,144,138 ; 232,202,164 ; 204,42,65];
% vintage modern
colors{11} = [140,35,24 ; 94,140,106 ; 136,166,94 ; 191,179,90 ; 242,196,90];

colors{12} = [ea_hex2rgb('#f8B195');ea_hex2rgb('#f67280');ea_hex2rgb('#c06c84');ea_hex2rgb('#6c5b7b');ea_hex2rgb('#355c7d')];

colors{13} = [ea_hex2rgb('#acdbc9');ea_hex2rgb('#dbebc2');ea_hex2rgb('#fdd2b5');ea_hex2rgb('#f7a7a6');ea_hex2rgb('#48b94')];

colors{14} = [ea_hex2rgb('#acdbc9');ea_hex2rgb('#dbebc2');ea_hex2rgb('#fdd2b5');ea_hex2rgb('#f7a7a6');ea_hex2rgb('#48b94')];

colors{15} = [ea_hex2rgb('#ef4566');ea_hex2rgb('#f69a9a');ea_hex2rgb('#f9cdae');ea_hex2rgb('#c8c8a9');ea_hex2rgb('#83ae9b')];

colors{16} = [ea_hex2rgb('#a6206a');ea_hex2rgb('#ec1c4b');ea_hex2rgb('#f16a43');ea_hex2rgb('#f7d969');ea_hex2rgb('#2f9395')];


% random
colors{length(colors)+1} = rand(30,3);

if ~exist('n','var') || isempty(n)
    nn=randperm(numel(colors));
    n=nn(1);
end

if max(colors{n}) >=1
    c = color_converter_inline(squeeze(colors{n}));
else
    c= squeeze(colors{n});
end
if plot
figure;
for a = 1:size(c,1);
    
    x=nan(1,size(c,1));
    x(a) = 1;
    b=bar(x);hold on;
    set(b,'EdgeColor',c(a,:),'FaceColor',c(a,:),'BarWidth',1)
end
figone(3);
title(['Palette Number ' num2str(n)])
set(gca,'YTick',[])
xlim([0.5 size(c,1)+.5])
end

c = repmat(c,[100,1]);
function new_color=color_converter_inline(ci)
new_color = ci./255;