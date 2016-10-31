%function randomColors : 
%
%      usage: cmap = randomColors(nColors)
%     author: Julien Besle
%       date: 19/04/2010
%       goal: creates a list of non adjacent RGB triplets
%       $Id$

function cmap = randomColors(nColors)
cycle_index = 0;                                            %cycle through the colors in order to avoid adjacent colors for adjacent events
color_cycle = [];                                                                           %
for i_event = 1:nColors                                                      %
    cycle_index = rem(cycle_index+round(nColors/2),nColors);  %
    while ~isempty(find(color_cycle==cycle_index, 1))                                            %
        cycle_index = rem(cycle_index+1, nColors);                           %
    end                                                                                     %
    color_cycle(i_event) = cycle_index;                                                     %    
end                                                                                         %
cmap = hsv(nColors)*.7;   %creates a matrix of number_of_events colors
cmap = cmap(color_cycle+1,:); %different events will have different colors
