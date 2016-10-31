% blue.m
%
%        $Id: blue.m 1950 2010-12-18 10:12:48Z julien $ 
%      usage: colorMap = blue(numberColors)
%         by: julien besle
%       date: 13/02/11
%    purpose: returns a colorMap from dark to bright MLR blue
%

function colorMap = blue2(numberColors)

minValue = .2;
hsvBlue = rgb2hsv([0 0.33 1]);
colorMap = repmat(hsvBlue,numberColors,1);
colorMap(:,3) = (1:numberColors)/numberColors*(1-minValue)+minValue;
colorMap = hsv2rgb(colorMap);