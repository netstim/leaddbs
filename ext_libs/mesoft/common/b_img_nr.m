%	b_img_nr.m 	reads parts of a sequence of Bruker images
% 			eg all images of one slice or of one echo time
%
%	usage: X = b_img_nr(fname,sz,nrep,rep, nsl,sl,ne,e,fmt);
%
%	input: 	fname   : string with filename of the binary Bruker file (2dseq)
%          	sz	: vector [m,n] size of the images in pixels, 
%			  		m:PHASE, n:READ;
%				nr    : total number of repetitions (time frames);
%				rep   : number of desired repetition(time frames);
%				nsl	: total number of slices;
%				sl	: number of desired slice;
%				ne	: total number of echos;
%				e	: number of desired echo;
%				fmt	: data format (default: 'long' = 32bit);
%	output: X	: Tensor containing images;
%			  size(X) = m n rep
%
%	possible abbreviations:
%			X = b_img_all(fname,sz,nr);	
%			X = b_img_all(fname,sz,rep,nsl,sl);	
%			X = b_img_all(fname,sz,rep,ne,e);	
%			X = b_img_all(fname,sz,rep,nsl,sl,ne,e);	
%		  		  
%	see also b_img_all, b_raw, b_raw_all, b_raw_epi 
% 	written by Claudia Oesterle;
%	tested on 02, 5.12.97;	


function X = b_img_nr(fname,sz,nr,rep,nsl,sl,ne,e,fmt)

% errors and tests
if (nargin==3) 
   ne = 1; e = 1; nsl = 1; sl = 1; end;
if (nargin==5) 
   ne = 1; e = 1; end;
if (nargin<3) 
   error('not enough input variables, specify at least fname, sz and rep'); end;
if (nargin<8),	fmt = 'long'; end;
if strcmp(fmt,'long'), o=4; 
elseif strcmp(fmt,'short'), o=2; end;
% end of errors and tests

% preallocation of X
X = zeros([sz,length(rep)]);
% end of preallocation of X

% read of the data
offs = o*sz(1)*sz(2)*(e-1)+o*sz(1)*sz(2)*(sl-1);
for r=1:(min(rep)-1);
   offs = offs+o*sz(1)*sz(2)*ne*nsl;
end;      


fid=fopen(fname,'rb','ieee-be');
if fid>=0
   tmpIndx =0;
   for r=rep;
      tmpIndx = tmpIndx + 1;
      fseek(fid,offs,-1); 
      tmp=fread(fid,[sz(2),sz(1)],fmt);
      X(:,:,tmpIndx) = tmp.';
      offs = offs+o*sz(1)*sz(2)*ne*nsl;
   end;      
else
   disp([fname,' could not be opened, check path'])
end;
fclose(fid);
X = squeeze(X);
% end of read of the data	

