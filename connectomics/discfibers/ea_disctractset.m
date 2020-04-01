classdef ea_disctractset < handle
    % Discriminative fiber set class to handle visualizations of sets of discriminative fibers in lead dbs resultfig / 3D Matlab figures
    % A. Horn

    properties (SetObservable)
        disctracts ea_disctract % set of discriminative tracts associated with this
        app % handle to discfiberexplorer app view
        modeltype % redundancy protocol only 1: t-tests/vta, 2: efields/spearman-R {'T-Tests / VTAs', 'Spearman''s Correlation / Efields'}
        leadgrouppath % redundancy protocol only, path to original lead group project
        M % content of Lead group analysis
    end
    
    properties (Access = private)
    end
    
    methods
        function obj=ea_disctractset(pobj) % class constructor, pobj needs to have M, ispositive and app entries. If it has more, it is treated as a loaded object with all entries.
            if isfield(pobj,'fibcell')
                obj=pobj; % e.g. loaded object from disk supplied.
            else
                % set up data from pobj.M and pobj.app
                obj.M=pobj.M;
                obj.app=pobj.app;
                obj.modeltype = obj.app.ModelSetupDropDown.Value;
                obj.leadgrouppath = pobj.leadgrouppath;
                
                inpobj.tractset=obj;
                inpobj.app=obj.app;
                inpobj.modeltype=obj.modeltype;
                inpobj.M=obj.M;
                obj.disctracts(1,1) = ea_disctract(inpobj);
                
            end
        end
        
        function calculate(obj)
            % calculate all visible tractsets:
            
            for set=1:size(obj.disctracts,1)
                    obj.disctracts(set,1).calculate;
            end
            
        end
        
        function draw(obj)
            for set=1:size(obj.disctracts,1)
                if obj.app.ShowPositiveFibersCheckBox.Value == 1
                    obj.disctracts(set,1).draw;
                end
                if obj.app.ShowNegativeFibersCheckBox.Value == 1
                    obj.disctracts(set,1).draw;
                end
            end
        end
    end
    methods (Static)
        function changeevent(~,event)
        end
    end
end
