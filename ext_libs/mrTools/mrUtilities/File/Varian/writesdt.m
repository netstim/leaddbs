% writesdt
%
%      usage: writesdt(filename,d);
%         by: justin gardner
%       date: 07/24/03
%    purpose: write sdt/spr format. d should be a structure like the one returned
%             by readsdt. Otherwise it can be an image array. If the data are complex
%             and edt/epr file will be saved
%    history:  11/26/2008 modified by Pei, just added at line 117 about endian part 
%                 endian:ieee-le or ieee-be don't have a space between field and data.
%                 this writesdt only can be used after readsdt, d should be
%                 a structure
%
function writesdt(filename,d)

% init some variables
t0 = clock;
MAXFIELDSIZE = 100;

% check arguments
if (nargin ~= 2)
  help writesdt;
  return
end

% if we are just getting image data then convert it into a stuct
if ~isstruct(d) && isnumeric(d)
  data = d;
  d = [];
  d.data = data;
  d.dataType = 'REAL';
  d.endian = 'ieee-be';
end

if ~isreal(d.data)
  disp(sprintf('(writesdt) Data is complex. Saving as an edt/epm file'));
  dataFileExt = 'edt';
  paramFileExt = 'epr';
  d.dataType = 'COMPLEX';
else
  dataFileExt = 'sdt';
  paramFileExt = 'spr';
end

% check to see if this looks like data
if ~isfield(d,'data')
  disp('(writesdt) No data in input structure');
  help writesdt;
  return
end

% okay, now we have the filenames of the two files
dataFilename = setext(filename,dataFileExt);
paramFilename = setext(filename,paramFileExt);

% check to see if the files already exists
fdata = fopen(dataFilename,'r');
if (fdata ~= -1)
  fclose(fdata);
  overwrite = input(sprintf('File %s already exists. Overwrite? (y/n) ',dataFilename),'s');  
  if (overwrite ~= 'y')
    return
  end
else
  fparam = fopen(paramFilename,'r');
  if (fparam ~= -1)
    fclose(fparam);
    overwrite = input(sprintf('File %s already exists. Overwrite? (y/n) ',paramFilename),'s');  
    if (overwrite ~= 'y')
      return
    end
  end
end

% get wordisze
datatypes = strvcat('BYTE','WORD','LWORD','REAL','COMPLEX');
datasizes = [           1      2       4      4         8];
matlabtypes = {'int8','int16','int32','float32','float32'};

if (isfield(d,'dataType'))
  whichtype = strmatch(d.dataType,datatypes);
  if (isempty(whichtype))
    disp(sprintf('(writesdt) Unrecogonized dataType %s. Setting to REAL.',d.dataType));
    whichtype = 4;
    d.dataType = 'REAL';
  end
  % set datasize
  d.wordsize = datasizes(whichtype);
else
  disp('(writesdt) Could not find dataType field. Setting to REAL.');
  d.dataType = 'REAL';
  whichtype = 4;
  d.wordsize = datasizes(whichtype);
end

% if dim does not exist then create it
if (~isfield(d,'dim'))
  d.dim = size(d.data);
  d.numDim = length(d.dim);
end

% if numDim does not exist then create it
if (~isfield(d,'numDim'))
  d.numDim = length(d.dim);
end

% if fidName does not exist
if (~isfield(d,'fidName'))
  if (isfield(d,'filepath'))
    d.fidName = d.filepath;
  else
    d.fidName = 'unknown';
  end
end

% calculate data size
d.filelen = prod(d.dim)*d.wordsize;

% open the data file
fdata = fopen(dataFilename,'w',d.endian);

if (fdata == -1)
  disp(sprintf('(writesdt) **** Could not open data file %s ****',dataFilename));
  return
end

%reshape data
data = reshape(d.data,prod(d.dim),1);

% if complex then we need to split into real and imaginary part and save
if strcmp(d.dataType,'COMPLEX')
  % split into real and imaginary
  realData = real(data);
  imagData = imag(data);
  % now put them back to gether interleaved
  data = [realData';imagData'];
  data = data(:);
  % we are writing out twice as much data
end

% write the data
count = fwrite(fdata,data,matlabtypes{whichtype});

% check that we wrote everything out
if strcmp(d.dataType,'COMPLEX')
  if count ~= 2*(d.filelen/datasizes(whichtype))
    disp(sprintf('(writesdt) **** Only wrote %i out of %i bytes ****',count,d.filelen));
  end
else
  if count ~= (d.filelen/datasizes(whichtype))
    disp(sprintf('(writesdt) **** Only wrote %i out of %i bytes ****',count,d.filelen));
  end
end

% close the file. we are done.
fclose(fdata);


% open the parameter file 
fparam = fopen(paramFilename,'w');
if (fparam == -1)
  disp(sprintf('(writesdt) **** Could not open parameter file %s ****',paramFilename));
  return
end

% try to save all parameters
allfields = fieldnames(d);
for i = 1:length(allfields)
  fielddata = eval(sprintf('d.%s',allfields{i}));
  % if the field is a string, then write that out
  if (ischar(fielddata))
    % fidName and sdtOrient don't have a space between field and data.
    %if (strcmp(allfields{i},'fidName') || strcmp(allfields{i},'sdtOrient') || strcmp(allfields{i},'endian'))
    %  fprintf(fparam,sprintf('%s:%s\n',allfields{i},fielddata));
    %else
    fprintf(fparam,sprintf('%s: %s\n',allfields{i},fielddata));
    %end
  % numeric fields
  elseif ((~iscell(fielddata)) && (~isstruct(fielddata)))
    % only save arrays (not larger matrices
    if ((length(size(fielddata)) == 2) && (min(size(fielddata)) == 1))
      % don't save if the array is too big    
      if (max(size(fielddata)) <= MAXFIELDSIZE)
	fprintf(fparam,sprintf('%s:',allfields{i}));
	% if the numbers are not integer
	if (sum(floor(fielddata)~=fielddata))
	  % than write as floating point
	  for j = 1:length(fielddata)
	    fprintf(fparam,sprintf(' %.6f',fielddata(j)));
	  end
	else
	  % otherwise write as ints
	  for j = 1:length(fielddata)
	    fprintf(fparam,sprintf(' %i',fielddata(j)));
	  end
	end	  
	fprintf(fparam,sprintf('\n'));
      end
    end
  end
end

fclose(fparam);

disp(sprintf('(writesdt) Took %0.2f seconds',etime(clock,t0)));