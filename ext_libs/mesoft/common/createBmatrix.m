% createBmatrix

% DE_schema12_b0_1000
% ika 030124
% script to calculate DTI fron Siem DTI images
% last update 021029
% example : DE_schema12_b0_1000('J.Rick-030124-1.5T','E:\zDataTemp\Rick_Jochen_20030124', DTI, 49, 40)
% explanation: Res =DE_schema12_b0_12(FileName,OutputDir, DTI, N_slices, tresholdAbs, DE_schema, b_factor)
% if b_factor is undefinedm, accumed to be =1000; %
% if b_factor=BMatrix with size 3x3xNumberOfDESteps, this BMatrix will be used for calculations 
% ika 031117 added option for exact b-matrix calculation  
% % 040726 - correction for b-matrix - to keep identical with fiber viever_tool
%
%   Wriiten by:
%   ------------------------------------------------------------
%   Dr.Kamil Il'yasov
%   Section of Medical Physics, Dept. of Diagnostic Radiology
%   University Medical Center of Freiburg
%   Hugstetter Str. 55, D-79106 Freiburg, Germany
%   Tel.: +49 761 270 5235
%   Fax : +49 761 270 3831
%   E_mail: Kamil.Ilyasov@uniklinik-freiburg.de
%   ------------------------------------------------------------
%% ika041121 added Liambda1, 2,3 in save DTD struct
%


function res =createBmatrix(DE_schema, b_factor)



DE_12 = [0 0 0; 1,0, .5; 0,.5, 1; .5 1 0; 1 .5 0; 0 1 .5; .5 0 1; 1 0 -.5; 0 -.5 1; -.5 1 0; 1 -.5 0; 0 1 -.5; -.5 0 1];
DE_6 = [0 0 0; 1 0 1; -1 0 1; 0 1 1 ; 0 1 -1; 1 1 0 ; -1 1 0];

% 		%%%NB!!!!!%%% moved 10 lines lower - done 
% 		DE_12(:,2) = -DE_12(:,2); % ika 030114  % 31.10.03 - why I have done that????
% 		DE_6(:,2) = -DE_6(:,2);  % 31.10.03 - accume that has to be done similat to line before???
% 		%%%NB!!!!!%%%
switch DT_Encoding_Schema
    case  'DE12', 
        DE=DE_12;
    case  'DE06', 
        DE =DE_6;         ...
        otherwise
        'unknown diffusion ensoding schema, only Siemens 6 & 12 diffusion encoding directions supported'
        return   
end
DE_size = size(DE);

DE(:,2) = - DE(:,2);% 31.10.03 - why I have done that????


%%%%%%%%%%%%%%%%555


DE_12 = [0 0 0; 1,0, .5; 0,.5, 1; .5 1 0; 1 .5 0; 0 1 .5; .5 0 1; 1 0 -.5; 0 -.5 1; -.5 1 0; 1 -.5 0; 0 1 -.5; -.5 0 1];
DE_6 = [0 0 0; 1 0 1; -1 0 1; 0 1 1 ; 0 1 -1; 1 1 0 ; -1 1 0];

%%%NB!!!!!%%%
DE_12(:,2) = -DE_12(:,2); % ika 030114  % 31.10.03 - why I have done that????
DE_6(:,2) = -DE_6(:,2);  % 31.10.03 - accume that has to be done similat to line before???
%%%NB!!!!!%%%

% if (DE_schema == 'DE12')
%     DE=DE_12;
% end
% 
%     if (DE_schema == 'DE06') % fixed bug on 040319 ika
%     DE =DE_6;
% else
%     'unknown diffusion ensoding schema, only Siemens 6 & 12 diffusion encoding directions supported'
%     break
% end


switch DE_schema
    case  'DE12', 
        DE=DE_12;
    case  'DE06', 
        DE =DE_6;         ...
        otherwise
        'unknown diffusion encoding schema, only Siemens 6 & 12 diffusion encoding directions supported'
        return   
end


DE_size = size(DE);
b_factorSize = size(b_factor)
if (length(b_factor)>1 )
    if (length(b_factorSize)==3 && b_factorSize(1)==3 && b_factorSize(2)==3) 
        B_tensor =b_factor;
        %%%NB!!!!!%%%
        B_tensor(2,1,:) = -B_tensor(2,1,:);
        B_tensor(1,2,:) = -B_tensor(1,2,:);
        B_tensor(2,3,:) = -B_tensor(2,3,:);
        B_tensor(3,2,:) = -B_tensor(3,3,:);
        %%%NB!!!!!%%%
    elseif(length(b_factorSize)==2 && b_factorSize(1)==1 ) 
        if b_factor(1)==0
            B_tensor(:,:,1)= zeros(3,3);
            b_start = 2;
        else
            b_start = 1;
        end
        i_count = b_start;
        for j= b_start:length(b_factor)
            for i = 2:DE_size(1)
                tmp =DE(i,:)'*DE(i,:);
                B_tensor(:,:,i_count)= (b_factor(j)/trace(tmp))*tmp;
                i_count = i_count +1:
            end
        end   
    else
        'b-factor has unsupported format'
        ret = 0;
        return
    end
else
    for i = 2:DE_size(1)
        tmp =DE(i,:)'*DE(i,:);
        B_tensor(:,:,i)= (b_factor/trace(tmp))*tmp;
    end
end;


N_DWIperSlice

b_factorSize = size(b_factor)
if (length(b_factor)>1 )
    if (length(b_factorSize)==3 && b_factorSize(1)==3 && b_factorSize(2)==3) 
        N_DWIperSlice = length(b_factor);
    elseif(length(b_factorSize)==2 && b_factorSize(1)==1 ) 
        if b_factor(1)==0
            N_DWIperSlice = 1+ (N_DWIperSlice-1)*(length(b_factor) -1)
        else
            N_DWIperSlice = 1+ (N_DWIperSlice-1)*(length(b_factor) -1)
        end
    else
        'b-factor has unsupported format'
        ret = 0;
        return
    end
end;
