%function M= open_brooker(fName, res, gradNo, vox)
%
%   function returns mrstruct
%   example result= open_brooker('data/vox2x2x4', [128 128 29], 66, [2 2 4])

function M= open_brooker(fName, res, gradNo, vox)

if ~exist('fName')
    fName= 'draft_data'
end
if ~exist('res')
    res= [128 128 29];
end
if ~exist('gradNo')
    gradNo= 66;
end

if ~exist('vox')
    vox= [2 2 4];
end

M= mrstruct_init('series3D');
M.vox= vox;

M.dataAy= zeros(res(1), res(2), res(3), gradNo);
for gradI=1:1:gradNo
%    strcat('GradientNo', num2str(gradI))
    for sliceI=1:1:res(3)
% 	usage: X = b_img_nr(fname,sz,nrep,rep, nsl,sl,ne,e,fmt);
% 
%  	input: 	fname   : string with filename of the binary Bruker file (2dseq)
%            	sz	: vector [m,n] size of the images in pixels, 
%  			  		m:PHASE, n:READ;
%  				nr    : total number of repetitions (time frames);
%  				rep   : number of desired repetition(time frames);
%  				nsl	: total number of slices;
%  				sl	: number of desired slice;
%  				ne	: total number of echos;
%  				e	: number of desired echo;
%  				fmt	: data format (default: 'long' = 32bit);
%  	output: X	: Tensor containing images;
        M.dataAy(:, :, sliceI, gradI)= b_img_nr(fName, [res(1), res(2)], gradNo, gradI, res(3), sliceI, 1, 1, 'long');
    end
end