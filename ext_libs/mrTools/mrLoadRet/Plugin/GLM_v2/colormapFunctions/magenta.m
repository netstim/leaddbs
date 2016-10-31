% magenta.m
%
%        $Id: magenta.m 1950 2010-12-18 10:12:48Z julien $ 
%      usage: colorMap = magenta(numberColors)
%         by: julien besle
%       date: 13/02/11
%    purpose: returns a colorMap from dark to bright MLR magenta
%

function colorMap = magenta(numberColors)

minValue = .2;
hsvBlue = rgb2hsv(color2RGB('magenta'));
colorMap = repmat(hsvBlue,numberColors,1);
colorMap(:,3) = (1:numberColors)/numberColors*(1-minValue)+minValue;
colorMap = hsv2rgb(colorMap);