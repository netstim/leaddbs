function genFDfromFTR(dtdname,ftrname,oversampling,selection,outname)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       genFDfromFTR(dtdname,ftrname,oversampling,selection,outname)
%   
%  generates fiber densities (colored and uncolred) and endpoint densities
%
%    dtdname - name of corresponding dtd (leaving blank causes to open a uigetfile)
%    ftrname - name of ftr (leaving blank causes to open a uigetfile)
%    oversampling - oversampling factor wrt to original resolution  (optionl)
%    selection - which bundle (empty for all) (optional)
%    outname - name prefix of resultrung files (optional)
%
%    2012 Marco Reisert UKLFR
%

if nargin < 4,
    selection = [];
end;

if nargin < 3,
    oversampling = 1;
end;

if nargin < 2,
    ftrname = [];
end;

if nargin < 1,
    dtdname = [];
end;



dtd = dtdstruct_read(dtdname);
mr = dtd.b0_image_struc;

[ftr err ftrname] = ftrstruct_read(ftrname);


if nargin < 5,
    [pp nn ] = fileparts(ftrname);
    outname = fullfile(pp,nn);
end;



display('generating FDmaps');
[fdrgb fd ep] = ftr2FDmaps(ftr,size(mr.dataAy),selection,oversampling);


display('saving FDmaps');

mr.edges(1:3,1:3) = mr.edges(1:3,1:3)/oversampling;

    
fdname = [outname '_fd'];
epname = [outname '_epd'];
fdrgbname = [outname '_fdrgb'];

mr.dataAy = fd; mrstruct_write(mr,fdname);
mr.dataAy = ep; mrstruct_write(mr,epname);
mr.dataAy = fdrgb; mrstruct_write(mr,fdrgbname);





    
