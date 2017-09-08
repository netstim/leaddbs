%% MetaTrajectory - Maker interface marking all sorts of classes represeting trajectories, and forcing them Copyable
% including non-straight line trajectories
%
% Andreas Husch
% Centre Hospitalier de Luxembourg, Dep. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicine
% 2017
% mail@andreashusch.de, fbernardpi@gmail.com

classdef MetaTrajectory < handle & matlab.mixin.Copyable
    events
        trajectoryChanged; %to be fired if a trajectory paramater is changed
    end
end