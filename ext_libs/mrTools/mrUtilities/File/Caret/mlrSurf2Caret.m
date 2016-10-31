% mlrSurf2Caret.m
%
%        $Id: mlrSurf2Caret.m 1529 2009-09-24 02:31:03Z justin $ 
%      usage: mlrSurf2Caret()
%         by: justin gardner
%       date: 08/19/08
%    purpose: Convert surf files from surfRelax to Caret
%
function retval = mlrSurf2Caret(filename)

% check arguments
if ~any(nargin == [0 1])
  help mlrSurf2Caret
  return
end

if ieNotDefined('filename'),filename = '';end
params = mrSurfViewer(filename);

if isempty(params),return,end

% load the anatomy header to look for talinfo
header = mlrImageHeaderLoad(params.anatomy);

% make sure we have AC point
if isempty(header.talInfo)
  warndlg('(mlrSurf2Caret) Please set AC point (no need to set any other tal points)');
  % get the AC point
  header.talInfo = talairach(params.anatomy);
  if isempty(header.talInfo),return,end
end

% AC point print out
disp(sprintf('(mlrSurf2Caret) Using [%i %i %i] as AC point',header.talInfo.AC(1),header.talInfo.AC(2),header.talInfo.AC(3)));

% load the outer surface
outerSurf = loadSurfOFF(params.outerSurface);
if isempty(outerSurf),return,end

% load the inner surface
innerSurf = loadSurfOFF(params.innerSurface);
if isempty(innerSurf),return,end

% make the fiducial surface (i.e. the one centered at the
% midpoint between the outer and inner surface)
surf = outerSurf;
surf.vtcs = (outerSurf.vtcs+innerSurf.vtcs)/2;

% and shift the points so AC is center
surf.vtcs(:,1) = surf.vtcs(:,1)-header.talInfo.AC(1);
surf.vtcs(:,2) = surf.vtcs(:,2)-header.talInfo.AC(2);
surf.vtcs(:,3) = surf.vtcs(:,3)-header.talInfo.AC(3);

% save as caret file
saveCaretFile(filename,surf,params);

%%%%%%%%%%%%%%%%%%%%%%%
%%   saveCaretFile   %%
%%%%%%%%%%%%%%%%%%%%%%%
function saveCaretFile(filename,surf,params)

surfFilename = sprintf('%s_fiducial',stripext(filename));

% save in ascii
encoding = 'ASCII';

% convert into structure that is called for by caret_save
% first make the coord file
M.encoding{1} = encoding;
M.num_nodes = surf.Nvtcs;
M.index = (0:(M.num_nodes-1))';
M.data = surf.vtcs;
% and save
caret_save(setext(surfFilename,'coord'),M);

% now make the topo file
clear M;
M.encoding{1} = encoding;
M.num_tiles = surf.Ntris;
M.data = surf.tris;
% and save
caret_save(setext(surfFilename,'topo'),M);

% now make the spec file
specFilename = setext(filename,'spec');
fid=fopen(specFilename,'w','ieee-be');
if isempty(fid)
  disp(sprintf('(mlrSurf2Caret) Could not create spec file %s',specFilename));
  return
end

% write the spec file
fprintf(fid,'BeginHeader\n');
fprintf(fid,'Encoding ASCII\n');
fprintf(fid,'EndHeader\n');
fprintf(fid,'tag-version 1\n');
fprintf(fid,'volume_anatomy_file %s\n',params.anatomy);
fprintf(fid,'FIDUCIALcoord_file %s\n',setext(surfFilename,'coord'));
fprintf(fid,'CLOSEDtopo_file %s\n',setext(surfFilename,'topo'));

% and close
fclose(fid);

% Following function is available for download from:
% http://www.bangor.ac.uk/~pss412/imaging/surface_stats.htm
% modified slightly to add tag-version 1 and number of nodes in topo file
%%%%%%%%%%%%%%%%%%%%
%%   caret_save   %%
% Jörn Diedrichsen %
%%%%%%%%%%%%%%%%%%%%
function caret_save(filename,M)
% function caret_save(filename,M)
% Save a M-file structure as a binary or ascii metric, coord, topo, pait file 
% M.encoding should be {'ASCII','BINARY'}
% --------------------------------------------------------------
% v.1.0 Joern Diedrichsen 04/11/28
if nargin==1
    M=filename;
    filename='';
end;
if isempty(filename)
    [F,P]=uiputfile('*.*','Save Structure as');
    filename = [P,F];
end
fid=fopen(filename,'w','ieee-be');
if (fid==-1)
    fprintf('Error opening file %s\n',filename);
    return;
end;
linefeed=sprintf('\n');

% Figure out the type of file we are loading: 
s=strfind(filename,'.');
type=filename(s(end)+1:end);

% --------------------------------------------------------------------
% Now switch the format and make the header depending on the type of file 
switch (type)
    % --------------------------------------------------------------------
    % Coordinate file
    case 'coord'
        if (~isfield(M,'header'))
            M.header{1}='BeginHeader';
            M.header{2}=['encoding ' M.encoding{1}];
	    M.header{3} = 'configuration_id FIDUCIAL';
            M.header{4}='EndHeader';
        end;
        for line=1:length(M.header) 
            fprintf(fid,'%s\n',M.header{line});
        end;
        % Format string: 
        if (strcmp(M.encoding,'ASCII'))
            fprintf(fid,'%d\n',M.num_nodes);
            M.data=[M.index M.data];
            format_str=['%i %3.3f %3.3f %3.3f\n'];
            fprintf(fid,format_str,M.data');
        else
            fwrite(fid,M.num_nodes,'int32');
            fwrite(fid,M.data','float32');
        end;
    % --------------------------------------------------------------------
    % Coordinate file
    case 'topo'
        if (~isfield(M,'header'))
        M.header{1}='BeginHeader';
        M.header{2}=['encoding ' M.encoding{1}];
        M.header{3}='filetype topo';
        M.header{4}='perimeter_id CLOSED';
        M.header{5}='EndHeader';
        M.header{6}='tag-version 1';
        end;
        for line=1:length(M.header) 
            fprintf(fid,'%s\n',M.header{line});
        end;
        M.data=M.data-1;
	fprintf(fid,'%i\n',M.num_tiles);
        if (strcmp(M.encoding,'ASCII'))
            format_str=['%d %d %d\n'];
            fprintf(fid,format_str,M.data');
        else
            fwrite(fid,M.num_tiles,'int32');
            fwrite(fid,M.data','int32');
        end;


    % --------------------------------------------------------------
    % Metric file 
    case 'metric'
        M.header={};
        M.header{1}='BeginHeader';
        M.header{2}=['encoding ' M.encoding{1}];
        M.header{3}='EndHeader';
        M.header{4}='tag-version 2';
        M.header{5}=sprintf('tag-number-of-nodes %d',size(M.data,1));
        M.header{6}=sprintf('tag-number-of-columns %d',size(M.data,2));
        M.header{7}=sprintf('tag-title');
        where=7;
        for i=1:length(M.column_name)
            M.header{where+i}=sprintf('tag-column-name %d %s',i-1,M.column_name{i});
        end;
        where=where+length(M.column_name);
        for i=1:length(M.column_name)
            M.header{where+i}=sprintf('tag-column-color-mapping %d %8.8f %8.8f',i-1,M.column_color_mapping(i,1),M.column_color_mapping(i,2));
        end;
        where=where+length(M.column_name);
        M.header{where+1}='tag-BEGIN-DATA';
        
        for line=1:length(M.header) 
            fprintf(fid,'%s\n',M.header{line});
        end;
        % Format string: 
        if (strcmp(M.encoding,'ASCII'))
            M.data=[M.index M.data];
            format_str=['%i' repmat([' %6.6f'],1,size(M.data,2)-1) '\n'];
            fprintf(fid,format_str,M.data');
        else
            fwrite(fid,M.data','float32');
        end;
end;
fclose(fid);


