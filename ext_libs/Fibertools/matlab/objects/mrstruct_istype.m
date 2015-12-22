%MRSTRUCT_ISTYPE assess validity or type of mrStruct.
%		truth = mrstruct_istype(mrStruct,typeStr)
%
%		If nargin==1 the function checks if mrStruct is a valid MR-Structure.
%		If typeStr is provided, the given mrStruct is checked for beeing a 
%		mrStruct of the given type (image,volume....).
%     For performance reasons the first argument will not be checked for
%     being a valid mrStruct.
%
%		Function returns 1 for 'true' and 0 for 'false'.
%
%		Examples:
%		truth = mrstruct_istype(myExperimentStruct);
%		truth = mrstruct_istype(myExperimentStruct,'spectrum');
%
%		Harald Fischer
%		1/00
%
%		Frederik Testud
%       20110118: added some lines to read in and recognize as mrStructs
%                         structs which have not the fields 'dim10'  and 'dim11'
%                         added previously
%
%
%		Frederik Testud
%       20140120: added some lines to read in and recognize as mrStructs
%                         structs which have not the field 'dim12' added previously
%
%		PC
%

function truth = mrstruct_istype(mrStruct,typeStr)


%%%%% init and error check
truth   = 0;
if nargin==2,
   if ~isstruct(mrStruct) || ~ischar(typeStr),
      warning('input argument have wrong type');
      return;
   end;
   
elseif nargin==1,
   if ~isstruct(mrStruct),
      return;
   end;

elseif nargin==0,
   warning('insufficient input arguments');
   return;
end
%%%%% End of: init and error check




%%%%% check if input struct is a valid mrStruct (usage 1)
if nargin>=1,
	try
    	tmpStruct      = mrstruct_init;
        mrStructNames  = fieldnames(mrStruct);
        tmpStructNames = fieldnames(tmpStruct);
        
        truth = true;
        for k = 1:length(mrStructNames),
            truth = truth & any(cellfun(@(x) strcmp(mrStructNames{k},x),tmpStructNames));
        end;
        
%         if ~any(strcmp(mrStructNames, 'dim10'))
%             tmpStruct = rmfield(tmpStruct, 'dim10');
%             tmpStruct = rmfield(tmpStruct, 'dim11');
%             tmpStruct = rmfield(tmpStruct, 'dim12');
%         elseif ~any(strcmp(mrStructNames, 'dim11')) 
%             tmpStruct = rmfield(tmpStruct, 'dim11');
%             tmpStruct = rmfield(tmpStruct, 'dim12');
%         elseif ~any(strcmp(mrStructNames, 'dim12'))
%             tmpStruct = rmfield(tmpStruct, 'dim12');
%         end
%         tmpStructNames = fieldnames(tmpStruct);
%         truth = isequal(tmpStructNames, mrStructNames);
    catch
        warning('input structure is not a valid mrStruct');
    end;
    %return;
end;
%%%%% End of: check if input struct is a valid mrStruct (usage 1)

%%%%% check cases (usage 2)
if (nargin>=2 && truth),
	if strcmp(typeStr,'image')
       if strcmp(mrStruct.dim1,'size_x') && strcmp(mrStruct.dim2,'size_y') && ...
          strcmp(mrStruct.dim3,'unused'),
          truth=1;
       else
          truth=0;
       end;
       
       
	elseif strcmp(typeStr,'imageEchos'),
       if strcmp(mrStruct.dim1,'size_x') && strcmp(mrStruct.dim2,'size_y') && ...
          strcmp(mrStruct.dim3,'echos') && strcmp(mrStruct.dim4,'unused'),
          truth=1; 
       else
          truth=0;
       end;
       
       
	elseif strcmp(typeStr,'volume'),
       if strcmp(mrStruct.dim1,'size_x') && strcmp(mrStruct.dim2,'size_y') && ...
          strcmp(mrStruct.dim3,'size_z') && strcmp(mrStruct.dim4,'unused'),
          truth=1; 
       else
          truth=0;
       end;
       
       
	elseif strcmp(typeStr,'volumeEchos'),
       if strcmp(mrStruct.dim1,'size_x') && strcmp(mrStruct.dim2,'size_y') && ...
          strcmp(mrStruct.dim3,'size_z') && strcmp(mrStruct.dim4,'echos') && ...
          strcmp(mrStruct.dim5,'unused'),   
          truth=1; 
       else
          truth=0;
       end;
       
       
	elseif strcmp(typeStr,'series2D'),
       if strcmp(mrStruct.dim1,'size_x') && strcmp(mrStruct.dim2,'size_y') && ...
          strcmp(mrStruct.dim3,'size_t') && strcmp(mrStruct.dim4,'unused'),
          truth=1; 
       else
          truth=0;
       end;
       
       
	elseif strcmp(typeStr,'series2DEchos'),
       if strcmp(mrStruct.dim1,'size_x') && strcmp(mrStruct.dim2,'size_y') && ...
          strcmp(mrStruct.dim3,'echos') && strcmp(mrStruct.dim4,'size_t') && ...
          strcmp(mrStruct.dim5,'unused'),
          truth=1; 
       else
          truth=0;
       end;
	
       
          
	elseif strcmp(typeStr,'series3D'),
       if strcmp(mrStruct.dim1,'size_x') && strcmp(mrStruct.dim2,'size_y') && ...
          strcmp(mrStruct.dim3,'size_z') && strcmp(mrStruct.dim4,'size_t') && ...
          strcmp(mrStruct.dim5,'unused'),
          truth=1; 
       else
          truth=0;
       end;
       
       
	elseif strcmp(typeStr,'series3DEchos'),
       if strcmp(mrStruct.dim1,'size_x') && strcmp(mrStruct.dim2,'size_y') && ...
          strcmp(mrStruct.dim3,'size_z') && strcmp(mrStruct.dim4,'echos') && ...
          strcmp(mrStruct.dim5,'size_t') && strcmp(mrStruct.dim6,'unused'),
          truth=1; 
       else
          truth=0;
       end;
       
       
	elseif strcmp(typeStr,'spectrum'),
       if strcmp(mrStruct.dim1,'spectral') && strcmp(mrStruct.dim2,'unused'),
          truth=1; 
       else
          truth=0;
       end;
       
       
	elseif strcmp(typeStr,'spectrum1D'),
       if strcmp(mrStruct.dim1,'size_x') && strcmp(mrStruct.dim2,'spectral') && ...
          strcmp(mrStruct.dim3,'unused'),
          truth=1; 
       else
          truth=0;
       end;
       
       
	elseif strcmp(typeStr,'spectrum2D')
       if strcmp(mrStruct.dim1,'size_x') && strcmp(mrStruct.dim2,'size_y') && ...
          strcmp(mrStruct.dim3,'spectral') && strcmp(mrStruct.dim4,'unused')
          truth=1; 
       else
          truth=0;
       end
       
       
	elseif strcmp(typeStr,'spectrum3D')
       if strcmp(mrStruct.dim1,'size_x') && strcmp(mrStruct.dim2,'size_y') && ...
          strcmp(mrStruct.dim3,'size_z') && strcmp(mrStruct.dim4,'spectral') && ...
          strcmp(mrStruct.dim5,'unused'),
          truth=1; 
       else
          truth=0;
       end
       
       elseif strcmp(typeStr,'DiffusionEchos2D'),
       truth = 0;
       if strcmp(mrStruct.dim1,'size_x') && strcmp(mrStruct.dim2,'size_y') && ...
         strcmp(mrStruct.dim3,'unused') && strcmp(mrStruct.dim4,'unused') && ...
         strcmp(mrStruct.dim5,'unused'),
         if isfield(mrStruct,'dim11'),
            if strcmp(mrStruct.dim11,'diffusion'),
                truth = 1;
            end;
         end;
       end
       
       elseif strcmp(typeStr,'DiffusionEchos3D'),
       truth = 0;
       if strcmp(mrStruct.dim1,'size_x') && strcmp(mrStruct.dim2,'size_y') && ...
          strcmp(mrStruct.dim3,'size_z') && strcmp(mrStruct.dim4,'unused') && ...
          strcmp(mrStruct.dim5,'unused') ,
          if isfield(mrStruct,'dim11'),
              if strcmp(mrStruct.dim11,'diffusion'),
                  truth = 1;
              end;
          end;      
       end

       elseif strcmp(typeStr,'Blade'),
       truth = 0;
       if strcmp(mrStruct.dim1,'size_x') && strcmp(mrStruct.dim2,'size_y') && ...
          strcmp(mrStruct.dim3,'size_z') && strcmp(mrStruct.dim4,'unused') && ...
          strcmp(mrStruct.dim5,'unused') ,
          if isfield(mrStruct,'dim12'),
              if strcmp(mrStruct.dim12,'ice_ind4'),
                  truth = 1;
              end;
          end;      
       end
       
%        elseif strcmp(typeStr,'Segment'),
%        truth = 0;
%        if strcmp(mrStruct.dim1,'size_x') && strcmp(mrStruct.dim2,'size_y') && ...
%           strcmp(mrStruct.dim3,'size_z') && strcmp(mrStruct.dim4,'unused') && ...
%           strcmp(mrStruct.dim5,'unused') ,
%           if isfield(mrStruct,'dim13'),
%               if strcmp(mrStruct.dim13,'segment'),
%                   truth = 1;
%               end;
%           end;      
%        end
       
       else
       warning('input type not recognized');
       truth=0;
       end
end
%%%%% End of: check cases (usage 2)
% Copyright (c) May 15th, 2001 by University of Freiburg, Dept. of Radiology, Section of Medical Physics