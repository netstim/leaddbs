function [ADC,FA,VectorF,DifT]=DTI(DTIdata,parameters)
% This function will perform DTI calculations on a certain
% DTI dataset.
%
% [ADC,FA,VectorF,DifT]=DTI(DTIdata,parameters)
%
% DTIdata:  A struct containing DTIdata(i).VoxelData, DTIdata(i).Gradient
%           and DTIdata(i).Bvalue of all the DTI datasets
% parameters: A struct containing, parameters.BackgroundTreshold and
%           parameters.WhiteMatterExtractionThreshold, parameters.textdisplay
%
% ADC: A 3D matrix with the Apparent Diffuse Coefficient (ADC)
% FA: A 3D matrix with the fractional anistropy (FA)
% VectorF: A 4D matrix with the (main) fiber direction in each pixel
% DifT: A 3D matrix with all  Diffusion tensors [Dxx,Dxy,Dxz,Dyy,Dyz,Dzz]
%
% Example,
%  see DTI_example.m
%
% This function is written by D.Kroon University of Twente (August 2008)

if(parameters.textdisplay), disp('Start DTI function'); pause(0.1); end
    
% Make a 4D matrix to store the different gradient voxel volumes
% (The constant 6 is just a minimum number of volumedatasets,
% this matrix will automaticaly expand if more volumedatasets are used.)
S=zeros([size(DTIdata(1).VoxelData) 6],'single');
% Make a 3D matrix to store the zero gradient voxel volume(s)
S0=zeros(size(DTIdata(1).VoxelData),'single');
% Make a matrix to store the gradients
H=zeros(6,3);
% Make a vector to store the different B-values (timing)
Bvalue=zeros(6,1);

% Read the input data (DTIdata), and seperate the zero gradient
% and other gradients into different matrices.
if(parameters.textdisplay), disp('Separate gradient and none gradient datasets'); pause(0.1); end
voxel0=0; voxelg=0;

for i=1:length(DTIdata)
    if(nnz(DTIdata(i).Gradient(:)==[0;0;0])==3)
        voxel0=voxel0+1;
        S0=S0+single(DTIdata(i).VoxelData);
    else
        voxelg=voxelg+1;
        S(:,:,:,voxelg)=single(DTIdata(i).VoxelData);
        H(voxelg,:)=single(DTIdata(i).Gradient);
        Bvalue(voxelg) = single(DTIdata(i).Bvalue);
    end
end
% The zero gradient matrix is the mean of all zero gradient datasets
S0=S0/voxel0;

% Free some memory
clear DTIdata;

if(parameters.textdisplay), disp('Create the b matrices'); pause(0.1); end
% Create the b matrices
% (http://www.meteoreservice.com/PDFs/Mattiello97.pdf)
b=zeros([3 3 size(H,1)]);
for i=1:size(H,1),
    b(:,:,i)=Bvalue(i)*H(i,:)'*H(i,:);
end

if(parameters.textdisplay), disp('Voxel intensity to absorption conversion'); pause(0.1); end
% Convert measurement intentensity into absorption (log)
Slog=zeros(size(S),'single');
for i=1:size(H,1),
    Slog(:,:,:,i)=log((S(:,:,:,i)./S0)+eps);
end

if(parameters.textdisplay), disp('Create B matrix vector [Bxx,2*Bxy,2*Bxz,Byy,2*Byz,Bzz]'); pause(0.1); end
% Sort all b matrices in to a vector Bv=[Bxx,2*Bxy,2*Bxz,Byy,2*Byz,Bzz];
Bv=squeeze([b(1,1,:),2*b(1,2,:),2*b(1,3,:),b(2,2,:),2*b(2,3,:),b(3,3,:)])';

% Create a matrix to store the Diffusion tensor of every voxel
% [Dxx,Dxy,Dxz,Dyy,Dyz,Dzz]
DifT=zeros([size(S0) 6],'single');
% Create a matrix to store the eigen values of every voxel
Y=zeros([size(S0) 3],'single');
% Create a matrix to store the fractional anistropy (FA)
FA=zeros(size(S0),'single');
% Create a matrix to store the Apparent Diffuse Coefficient (ADC)
ADC=zeros(size(S0),'single');
% Create a maxtrix to store the (main) fiber direction in each pixel
VectorF=zeros([size(S0) 3],'single');

if(parameters.textdisplay), disp('Calculate Diffusion tensor, eigenvalues, and other parameters of each voxel'); pause(0.1); end
% Loop through all voxel coordinates
for x=1:size(S0,1),
    for y=1:size(S0,2),
        for z=1:size(S0,3),
            
            % Only process a pixel if it isn't background
            %if(S0(x,y,z)>parameters.BackgroundTreshold)
                
                % Calculate the Diffusion tensor [Dxx,Dxy,Dxz,Dyy,Dyz,Dzz],
                % with a simple matrix inverse.
                Z=-squeeze(Slog(x,y,z,:));
                M=Bv\Z;
               
                % The DiffusionTensor (Remember it is a symetric matrix,
                % thus for instance Dxy == Dyx)
                DiffusionTensor=[M(1) M(2) M(3); M(2) M(4) M(5); M(3) M(5) M(6)];

                % Calculate the eigenvalues and vectors, and sort the 
                % eigenvalues from small to large
                [EigenVectors,D]=eig(DiffusionTensor); EigenValues=diag(D);
                [t,index]=sort(EigenValues); 
                EigenValues=EigenValues(index); EigenVectors=EigenVectors(:,index);
                EigenValues_old=EigenValues;
                
                % Regulating of the eigen values (negative eigenvalues are
                % due to noise and other non-idealities of MRI)
                if((EigenValues(1)<0)&&(EigenValues(2)<0)&&(EigenValues(3)<0)), EigenValues=abs(EigenValues);end
                if(EigenValues(1)<=0), EigenValues(1)=eps; end
                if(EigenValues(2)<=0), EigenValues(2)=eps; end
                
                % Apparent Diffuse Coefficient
                ADCv=(EigenValues(1)+EigenValues(2)+EigenValues(3))/3;
                
                % Fractional Anistropy (2 different definitions exist)
                % First FA definition:
                %FAv=(1/sqrt(2))*( sqrt((EigenValues(1)-EigenValues(2)).^2+(EigenValues(2)-EigenValues(3)).^2+(EigenValues(1)-EigenValues(3)).^2)./sqrt(EigenValues(1).^2+EigenValues(2).^2+EigenValues(3).^2) );
                % Second FA definition:
                FAv=sqrt(1.5)*( sqrt((EigenValues(1)-ADCv).^2+(EigenValues(2)-ADCv).^2+(EigenValues(3)-ADCv).^2)./sqrt(EigenValues(1).^2+EigenValues(2).^2+EigenValues(3).^2) );
                
                % Store the results of this pixel in the volume matrices
                ADC(x,y,z)=ADCv;
                Y(x,y,z,:)=EigenValues;
                DifT(x,y,z,:)=[DiffusionTensor(1:3) DiffusionTensor(5:6) DiffusionTensor(9)];
                % Only store the FA and fiber vector of a voxel, if it exceed an anistropy treshold
               FA(x,y,z)=FAv;
                if(FAv>parameters.WhiteMatterExtractionThreshold)
                    VectorF(x,y,z,:)=EigenVectors(:,end)*EigenValues_old(end);
                end
            %end
        end
    end
end
if(parameters.textdisplay), disp('DTI function Finished'); pause(0.1); end





