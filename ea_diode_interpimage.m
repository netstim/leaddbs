function interVal=ea_diode_interpimage(img,yx,zpad)
% BiLinear interpolation using 4 pixels around the target location with ceil convention
% RGB = 1 for gray scale images.
% img can be a single layer matrix or a RGB layer colored image
% yx =[y_value, x_value]; It can be either horizontal or vertical vector
%
% zpad is a boolean variable. if true, zeros are used for pixel values
% outside of the given img. If false, the nearest edge value is repeated.
%
% Example:
% [m,n]=meshgrid(1:3);img=[m+n]
% --> 2     3     4
%     3     4     5
%     4     5     6
% interpImg(img,[2.4,2.2])
% --> 4.6
%
%                                   Disi A, Sep,16th,2013
%                                   adis@mit.edu

% Update with Bug fix suggested by Sanjuro
% There is bug. When x or y is under 1 or over max size (eg. [2.4, 0]) and with zpad = false, then you try get not integer indexed pixel in line 60: 
% pixVal=img(y0,x0,:);

% should be: 
% pixVal=img(ceil(y0),ceil(x0),:);

if nargin<4,RGB=ndims(img);RGB(RGB<3)=1; end
if nargin<3,zpad=true; end

yx0=floor(yx);
wt=yx-yx0; wtConj=1-wt;
interTop=wtConj(2)*pixLookup(img,yx0(1),yx0(2),zpad,RGB)+wt(2)*pixLookup(img,yx0(1),yx(2),zpad,RGB);
interBtm=wtConj(2)*pixLookup(img,yx(1),yx0(2),zpad,RGB)+wt(2)*pixLookup(img,yx(1),yx(2),zpad,RGB);
interVal=wtConj(1)*interTop+wt(1)*interBtm;

end


function pixVal=pixLookup(img,y,x,zpad,RGB)
% This helper function looks up a pixel value from a given input image
% img is the input image (RGB or Grayscale)
% yx is the coordinate and repEdge tells the condition for pixel values out
% side of img (Use round up convention)
% For grayscale use RGB =1
if nargin<4,RGB=3;end

pixVal=zeros(1,1,RGB); %Initialize the pixel 

if nargin<3
    zpad=true; %pad with black pixels
end

if RGB==3
    [ROW,COL,~]=size(img);
else
    [ROW,COL]=size(img);
end
% If the pixel value is outside of image given
if (x<=0)||(x>COL)||(y<=0)||(y>ROW) 
    if zpad
        pixVal(:)=0;
    else
        y0=y;x0=x;
        y0(y0<1)=1; x0(x0<1)=1;
        y0(y0>ROW)=ROW;x0(x0>COL)=COL;
        % Bug fix suggested by Sanjuro
        pixVal=img(ceil(y0),ceil(x0),:);
    end
else
    pixVal=img(ceil(y),ceil(x),:);
end

end