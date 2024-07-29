classdef ea_bfigure < handle
    properties
        Vertices = [];
        Faces = [];
        Lines = [];
        SurfaceData = struct('X', [], 'Y', [], 'Z', []);
        StreamtubeData = struct('X', [], 'Y', [], 'Z', [], 'U', [], 'V', [], 'W', []);
        Lights = {};
        CamLights = {};
        RunInBackground = false;  % Default to running in background
        BlenderPath = fullfile(ea_getearoot,'ext_libs','blender','Blender.app','Contents','MacOS','Blender');  % Path to Blender executable
    end
    
    methods
        function obj = ea_bfigure()
            % Constructor
        end
        
        function patch(obj, varargin)
            % Handle patch command
            p = inputParser;
            addParameter(p, 'Vertices', []);
            addParameter(p, 'Faces', []);
            parse(p, varargin{:});
            
            obj.Vertices = p.Results.Vertices;
            obj.Faces = p.Results.Faces;
            
            % Save data and call Blender script
            obj.saveData('patch');
            obj.callBlenderScript('patch');
        end
        
        function plot3(obj, x, y, z)
            % Handle plot3 command
            obj.Lines = [x(:), y(:), z(:)];
            
            % Save data and call Blender script
            obj.saveData('plot3');
            obj.callBlenderScript('plot3');
        end
        
        function surf(obj, X, Y, Z)
            % Handle surf command
            obj.SurfaceData.X = X;
            obj.SurfaceData.Y = Y;
            obj.SurfaceData.Z = Z;
            
            % Save data and call Blender script
            obj.saveData('surf');
            obj.callBlenderScript('surf');
        end
        
        function surface(obj, X, Y, Z)
            % Handle surface command (similar to surf)
            obj.surf(X, Y, Z);
        end
        
        function streamtube(obj, X, Y, Z, U, V, W)
            % Handle streamtube command
            obj.StreamtubeData.X = X;
            obj.StreamtubeData.Y = Y;
            obj.StreamtubeData.Z = Z;
            obj.StreamtubeData.U = U;
            obj.StreamtubeData.V = V;
            obj.StreamtubeData.W = W;
            
            % Save data and call Blender script
            obj.saveData('streamtube');
            obj.callBlenderScript('streamtube');
        end
        
        function light(obj, varargin)
            % Handle light command
            p = inputParser;
            addParameter(p, 'Position', []);
            addParameter(p, 'Type', 'SUN');
            addParameter(p, 'Energy', 3);
            parse(p, varargin{:});
            
            lightData = struct('position', p.Results.Position, 'type', p.Results.Type, 'energy', p.Results.Energy);
            obj.Lights{end+1} = lightData;
            
            % Save light data and call Blender script
            obj.saveData('light');
            obj.callBlenderScript('light');
        end
        
        function camlight(obj, varargin)
            % Handle camlight command
            p = inputParser;
            addParameter(p, 'Position', []);
            addParameter(p, 'Type', 'SUN');
            addParameter(p, 'Energy', 3);
            parse(p, varargin{:});
            
            camlightData = struct('position', p.Results.Position, 'type', p.Results.Type, 'energy', p.Results.Energy);
            obj.CamLights{end+1} = camlightData;
            
            % Save camlight data and call Blender script
            obj.saveData('camlight');
            obj.callBlenderScript('camlight');
        end
        
        function saveData(obj, plotType)
            % Save data to CSV files based on plot type
            switch plotType
                case 'patch'
                    csvwrite('vertices.csv', obj.Vertices);
                    csvwrite('faces.csv', obj.Faces);
                case 'plot3'
                    csvwrite('lines.csv', obj.Lines);
                case 'surf'
                    csvwrite('X.csv', obj.SurfaceData.X);
                    csvwrite('Y.csv', obj.SurfaceData.Y);
                    csvwrite('Z.csv', obj.SurfaceData.Z);
                case 'streamtube'
                    csvwrite('X.csv', obj.StreamtubeData.X);
                    csvwrite('Y.csv', obj.StreamtubeData.Y);
                    csvwrite('Z.csv', obj.StreamtubeData.Z);
                    csvwrite('U.csv', obj.StreamtubeData.U);
                    csvwrite('V.csv', obj.StreamtubeData.V);
                    csvwrite('W.csv', obj.StreamtubeData.W);
                case 'light'
                    jsonStr = jsonencode(obj.Lights);
                    fid = fopen('lights.json', 'w');
                    fwrite(fid, jsonStr, 'char');
                    fclose(fid);
                case 'camlight'
                    jsonStr = jsonencode(obj.CamLights);
                    fid = fopen('camlights.json', 'w');
                    fwrite(fid, jsonStr, 'char');
                    fclose(fid);
            end
        end
        
        function callBlenderScript(obj, plotType)
            % Define the path to the Blender script
            blenderScript = fullfile(pwd, 'ea_blenderplot.py');
            
            % Determine whether to run Blender in background or not
            if obj.RunInBackground
                system([obj.BlenderPath, ' --background --python ', blenderScript, ' -- ', plotType]);
            else
                system([obj.BlenderPath, ' --python ', blenderScript, ' -- ', plotType]);
            end
        end
    end
end