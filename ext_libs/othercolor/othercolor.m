function c = othercolor(n,m)
% OTHERCOLOR alternative colormaps from various sources
%
%   OTHERCOLOR(N,M) returns a M-by-3 matrix contaning a colormap given by
%   the name N.  OTHERCOLOR, by itself, is the same length as the current
%   figure's colormap.  If no figure exists, MATLAB creates one.
%
%   Supported colormaps are stored in a colorData.mat.  And you can easily
%   add your own: they are just standard MATLAB RGB colormap matrices.
%
%   From: http://geography.uoregon.edu/datagraphics/color_scales.htm
%   BrBu_10     BuDOr_18    BuGr_14     BuOr_10     Bu_10       GrMg_16     
%   BrBu_12     BuDRd_12    BuGy_8      BuOr_12     Bu_7        RdYlBu_11b  
%   BuDOr_12    BuDRd_18    BuOrR_14    BuOr_8      Cat_12      StepSeq_25  
%   
%   From: http://www.colorbrewer2.org
%   Accent3     GnBu6       Paired4     PuOr4       RdYlBu10    Set37       
%   Accent4     GnBu7       Paired5     PuOr5       RdYlBu11    Set38       
%   Accent5     GnBu8       Paired6     PuOr6       RdYlBu3     Set39       
%   Accent6     GnBu9       Paired7     PuOr7       RdYlBu4     Spectral10  
%   Accent7     Greens3     Paired8     PuOr8       RdYlBu5     Spectral11  
%   Accent8     Greens4     Paired9     PuOr9       RdYlBu6     Spectral3   
%   Blues3      Greens5     Pastel13    PuRd3       RdYlBu7     Spectral4   
%   Blues4      Greens6     Pastel14    PuRd4       RdYlBu8     Spectral5   
%   Blues5      Greens7     Pastel15    PuRd5       RdYlBu9     Spectral6   
%   Blues6      Greens8     Pastel16    PuRd6       RdYlGn10    Spectral7   
%   Blues7      Greens9     Pastel17    PuRd7       RdYlGn11    Spectral8   
%   Blues8      Greys3      Pastel18    PuRd8       RdYlGn3     Spectral9   
%   Blues9      Greys4      Pastel19    PuRd9       RdYlGn4     YlGn3       
%   BrBG10      Greys5      Pastel23    Purples3    RdYlGn5     YlGn4       
%   BrBG11      Greys6      Pastel24    Purples4    RdYlGn6     YlGn5       
%   BrBG3       Greys7      Pastel25    Purples5    RdYlGn7     YlGn6       
%   BrBG4       Greys8      Pastel26    Purples6    RdYlGn8     YlGn7       
%   BrBG5       Greys9      Pastel27    Purples7    RdYlGn9     YlGn8       
%   BrBG6       OrRd3       Pastel28    Purples8    Reds3       YlGn9       
%   BrBG7       OrRd4       PiYG10      Purples9    Reds4       YlGnBu3     
%   BrBG8       OrRd5       PiYG11      RdBu10      Reds5       YlGnBu4     
%   BrBG9       OrRd6       PiYG3       RdBu11      Reds6       YlGnBu5     
%   BuGn3       OrRd7       PiYG4       RdBu3       Reds7       YlGnBu6     
%   BuGn4       OrRd8       PiYG5       RdBu4       Reds8       YlGnBu7     
%   BuGn5       OrRd9       PiYG6       RdBu5       Reds9       YlGnBu8     
%   BuGn6       Oranges3    PiYG7       RdBu6       Set13       YlGnBu9     
%   BuGn7       Oranges4    PiYG8       RdBu7       Set14       YlOrBr3     
%   BuGn8       Oranges5    PiYG9       RdBu8       Set15       YlOrBr4     
%   BuGn9       Oranges6    PuBu3       RdBu9       Set16       YlOrBr5     
%   BuPu3       Oranges7    PuBu4       RdGy10      Set17       YlOrBr6     
%   BuPu4       Oranges8    PuBu5       RdGy11      Set18       YlOrBr7     
%   BuPu5       Oranges9    PuBu6       RdGy3       Set19       YlOrBr8     
%   BuPu6       PRGn10      PuBu7       RdGy4       Set23       YlOrBr9     
%   BuPu7       PRGn11      PuBu8       RdGy5       Set24       YlOrRd3     
%   BuPu8       PRGn3       PuBu9       RdGy6       Set25       YlOrRd4     
%   BuPu9       PRGn4       PuBuGn3     RdGy7       Set26       YlOrRd5     
%   Dark23      PRGn5       PuBuGn4     RdGy8       Set27       YlOrRd6     
%   Dark24      PRGn6       PuBuGn5     RdGy9       Set28       YlOrRd7     
%   Dark25      PRGn7       PuBuGn6     RdPu3       Set310      YlOrRd8     
%   Dark26      PRGn8       PuBuGn7     RdPu4       Set311      YlOrRd9     
%   Dark27      PRGn9       PuBuGn8     RdPu5       Set312      
%   Dark28      Paired10    PuBuGn9     RdPu6       Set33       
%   GnBu3       Paired11    PuOr10      RdPu7       Set34       
%   GnBu4       Paired12    PuOr11      RdPu8       Set35       
%   GnBu5       Paired3     PuOr3       RdPu9       Set36       
%
%   From Mathematica:
%     MCMYKcolors           MIndexed45            Mdarkterrain          
%     MHTML                 MIndexed46            Mdeepseacolors        
%     MIndexed1             MIndexed47            Mfallcolors           
%     MIndexed10            MIndexed48            Mfruitpunchcolors     
%     MIndexed11            MIndexed49            Mfuchsiatones         
%     MIndexed12            MIndexed5             Mgeologicages         
%     MIndexed13            MIndexed50            Mgraytones            
%     MIndexed14            MIndexed51            Mgrayyellowtones      
%     MIndexed15            MIndexed52            Mgreenbrownterrain    
%     MIndexed16            MIndexed53            Mgreenpinktones       
%     MIndexed17            MIndexed54            Mhypsometrictints     
%     MIndexed18            MIndexed55            Mislandcolors         
%     MIndexed19            MIndexed56            Mlakecolors           
%     MIndexed2             MIndexed57            Mlegacy               
%     MIndexed20            MIndexed58            Mlighttemperaturemap  
%     MIndexed21            MIndexed59            Mlightterrain         
%     MIndexed22            MIndexed6             Mmintcolors           
%     MIndexed23            MIndexed60            Mneoncolors           
%     MIndexed24            MIndexed61            Mpastel               
%     MIndexed25            MIndexed62            Mpearlcolors          
%     MIndexed26            MIndexed7             Mpigeontones          
%     MIndexed27            MIndexed8             Mplumcolors           
%     MIndexed28            MIndexed9             Mrainbow              
%     MIndexed29            Malpinecolors         Mredbluetones         
%     MIndexed3             Maquamarine           Mredgreensplit        
%     MIndexed30            Marmycolors           Mrosecolors           
%     MIndexed31            Matlanticcolors       Mrusttones            
%     MIndexed32            Matoms                Msandyterrain         
%     MIndexed33            Mauroracolors         Msiennatones          
%     MIndexed34            Mavocadocolors        Msolarcolors          
%     MIndexed35            Mbeachcolors          Msouthwestcolors      
%     MIndexed36            Mblackbodyspectrum    Mstarrynightcolors    
%     MIndexed37            Mbluegreenyellow      Msunsetcolors         
%     MIndexed38            Mbrasstones           Mtemperaturemap       
%     MIndexed39            Mbrightbands          Mthermometercolors    
%     MIndexed4             Mbrowncyantones       Mvalentinetones       
%     MIndexed40            Mcandycolors          Mvisiblespectrum      
%     MIndexed41            Mcherrytones          Mwatermeloncolors     
%     MIndexed42            Mcoffeetones          Mwebsafe              
%     MIndexed43            Mdarkbands            
%     MIndexed44            Mdarkrainbow     
%
%   Usage:
%   A typical 3D plot:
%   >> [X,Y,Z] = peaks(30);
%   >> surfc(X,Y,Z)
%   >> colormap(othercolor('RdYlBu_11b'))
%   >> colorbar
%   >> axis([-3 3 -3 3 -10 5])
%
%   To get the list of available colormaps in a cellarray:
%   >> colormapNames = othercolor();
%
%   Iterate through colormaps (enter to move to next, ctrl+c to exit loop)
%   >> l = othercolor; for i=1:length(l), colormap(othercolor(i));pause;end
%
%   Plot the first 50 colormaps
%   >> colors = othercolor();
%   >> l = 50;
%   >> for i=1:l
%   >> subplot(ceil(l/10),10,i);
%   >> c = othercolor(i);
%   >> imagesc(reshape(c,1,size(c,1),size(c,2)));
%   >> title(char(colors(i)),'interpreter','none');
%   >> axis off;
%   >> end
%
%   Author: Joshua Atkins
%   Date: March 1, 2011

types = who('-file','colorData.mat');

% if no colormap is choosen then display available colormaps
if nargin < 1,
    c = types;
else
    % default number of points
    if nargin < 2, m = size(get(gcf,'colormap'),1); end

    % allows numerical indexing
    if isnumeric(n), n = char(types(n)); end
        
    % load color data
    data = load('colorData.mat',n);
    if isempty(fieldnames(data))
        c = [];
    else
        c = interp1(linspace(0,1,size(data.(n),1)),data.(n),linspace(0,1,m),'pchip');
        c(c<0) = 0;
        c(c>1) = 1;
    end
end