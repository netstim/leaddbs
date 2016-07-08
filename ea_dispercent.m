function  ea_dispercent(varargin)
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

percent=round(varargin{1}*100);

if nargin==2
    if strcmp(varargin{2},'end')
        fprintf(1,[varargin{2},':     ']); 
        fprintf('\n')
        fprintf('\n')          
    else
        fprintf(1,[varargin{2},':     ']);       
    end
else
    fprintf(1,[repmat('\b',1,(length(num2str(percent))+1)),'%d','%%'],percent);
end