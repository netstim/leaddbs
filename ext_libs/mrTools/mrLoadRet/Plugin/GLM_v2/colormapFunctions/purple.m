% purple.m
%
%        $Id: purple.m 1950 2010-12-18 10:12:48Z julien $ 
%      usage: colorMap = purple(numberColors)
%         by: julien besle
%       date: 13/02/11
%    purpose: returns a colorMap from dark to bright MLR purple
%

function colorMap = purple(numberColors)

minValue = .2;
hsvPurple = rgb2hsv(color2RGB('purple'));
colorMap = repmat(hsvPurple,numberColors,1);
colorMap(:,3) = (1:numberColors)/numberColors*(1-minValue)+minValue;
colorMap = hsv2rgb(colorMap);