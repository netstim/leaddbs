%% Abstract Class - configurablem marker interface
%
% A subclass has to implement a config panel to allow 
% graphical confiuration of it's features
%
% Andreas Husch
% Centre Hospitalier de Luxembourg, Dep. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicine
% 2013 - 2017
% mail@andreashusch.de, husch.andreas@chl.lu

classdef configurable < handle
    methods
        panel = getConfigPanel(this, varargin);
    end
end