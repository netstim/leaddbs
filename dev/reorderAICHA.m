function reorderAICHA

f=fopen('AICHA (Joliot 2015).txt')
fn=fopen('AICHA reordered (Joliot 2015).txt','w');

A=textscan(f,'%d %s')
nii=ea_load_nii('AICHA (Joliot 2015).nii')
nunii=nii;
nunii.img(:)=0;

for reg=1:length(A{2})
    
 
if isodd(reg)
nureg(reg)=round(reg/2);
else
    nureg(reg)=192+round(reg/2);
end
end
fnureg=[1:2:384,2:2:384];

cnt=1;
for reg=1:384
    
    nunii.img(nii.img==reg)=nureg(reg);
    fprintf(fn,'%d %s\n',reg,A{2}{fnureg(reg)});
    cnt=cnt+1;
end

nunii.fname='AICHA reordered (Joliot 2015).nii';
spm_write_vol(nunii,nunii.img);





function iso = isodd(x)
%ISODD True for odd numbers.
%   ISODD(X) is 1 for the elements of X that are odd and 0 for even elements.
%
%   Class support for input X:
%      float: double, single
%      integer: [u]int8, [u]int16, [u]int32, [u]int64
%
%   An error is raised if X is a float and is:
%   1) too large or small to determine if it is odd or even (this includes -Inf
%      and +Inf), or
%   2) contains a fraction (e.g. isodd(1.5) raises an error), or
%   3) NaN (not-a-number).
%
%   See also ISEVEN.

% By Ulf Carlberg 29-MAY-2011.

if isa(x, 'float') && isreal(x)
    ax = abs(double(x));
    if isa(x, 'single') && any(ax(:)>=2^24)
        % Also x=Inf or x=-Inf is will be taken care of here.
        error('isodd:InputOutOfRange', ...
            'The maximum value of abs(X) allowed is 2^24-1 for singles.');
    end
    if any(ax(:)>=2^53)
        % Also x=Inf or x=-Inf is will be taken care of here.
        error('isodd:InputOutOfRange', ...
            'The maximum value of abs(X) allowed is 2^53-1 for doubles.');
    end
    if any(isnan(x(:)))
        error('isodd:InputNaN', 'No entries of X are allowed to be NaN.');
    end
    if any(floor(ax(:))~=ax(:))
        error('isodd:InputNotInt', 'All entries of X must be integers.');
    end
    iso = logical(bitget(ax, 1));
elseif isa(x, 'integer')
    if isa(x, 'uint64')
        % MOD can not handle 64 bit numbers on some Matlab versions, but BITGET
        % works on unsigned 64-bit integers.
        iso = logical(bitget(x, 1));
    elseif isa(x, 'int64')
        % MOD can not handle 64 bit numbers on some Matlab versions, but ABS
        % works, and BITGET works on unsigned 64-bit integers.
        ax = uint64(abs(x));
        ax(x==-2^63) = 0;  % "the weird number" is even
        iso = logical(bitget(ax, 1));
    else
        iso = logical(mod(x, 2));
    end
else
    error('isodd:UnsupportedClass', 'Unsupported class.');
end