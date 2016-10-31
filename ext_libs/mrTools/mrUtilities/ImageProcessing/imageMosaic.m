function mosaic = imageMosaic(volume1, volume2, N)
% imgMosaic - generates a checkerboard mosaic from the two input images
%             with blocks of size NxN (default 10)
%
%       mosaic = regMosaic(volume1, volume2, <N>)
%
% Oscar Nestares - 5/99
% jb 10/2010 - made into a function
%            - extended to 3D volumes

[Ny Nx Nz] = size(volume1);

if nargin<3
   N=round(max(Ny,Nx)/15);
end

mosaic = zeros(Ny,Nx,Nz);

for k = 1:Nz
   slice1 = volume1(:,:,k);
   slice2 = volume2(:,:,k);
    
   slice1(isnan(slice1)) = 0;
   slice2(isnan(slice2)) = 0;

   basic = [ones(N) zeros(N); zeros(N) ones(N)];
   check = repmat(basic, ceil(Ny/N), ceil(Nx/N));
   check = check(1:Ny,1:Nx);

   mosaic(:,:,k) = check.*slice1 + (~check).*slice2;

% puts a slightly different contrast in one image than in the other
% mosaic = check.*slice1*0.9 + (~check).*slice2;

end
