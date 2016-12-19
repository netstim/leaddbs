function Z = imresize3D_diff(A, dim, method)

% for other methods than cubic consult 'doc interp3'

if nargin<=2
    method = 'cubic';
end

[Nr, Nc, Ns] = size(A);
Mr = dim(1);
Mc = dim(2);
Ms = dim(3);

[X, Y, Z] = meshgrid(1:Nc, 1:Nr, 1:Ns);

% Ir = [1:Nr/Mr:Nr] + (Nr - (Mr-1)*Nr/Mr - 1)/2;
% Ic = [1:Nc/Mc:Nc] + (Nc - (Mc-1)*Nc/Mc - 1)/2;
% Is = [1:Ns/Ms:Ns] + (Ns - (Ms-1)*Ns/Ms - 1)/2;
Ir = linspace(1,Nr,Mr);
Ic = linspace(1,Nc,Mc);
Is = linspace(1,Ns,Ms);

[Xi,Yi,Zi] = meshgrid(Ic,Ir,Is);

Z = interp3(X,Y,Z,A,Xi,Yi,Zi,method,0); % values lieing outside of the grid will be set to zero
