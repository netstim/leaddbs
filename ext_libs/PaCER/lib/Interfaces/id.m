%% id -  add a unique id to an object
% e.g. to be used as a key
%
% Florian Bernard, Andreas Husch
% Centre Hospitalier de Luxembourg, Dep. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicne
% 2014 - 2017
% mail@andreashusch.de, fbernardpi@gmail.com

classdef id < handle
    properties ( GetAccess = 'public', SetAccess = 'private' )
        objectId
    end
    
    methods ( Access = 'protected' )
        function this = id()
            this.objectId = id.increment();
        end
    end
    
    methods ( Static, Access = 'private' )
        function result = increment()
            persistent stamp;
            if isempty( stamp )
                stamp = 0;
            end
            stamp = stamp + uint32(1);
            result = stamp;
        end
    end
end