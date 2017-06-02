function mat=ea_antsmat2mat(afftransform,fixed)
% followed instructions from
% https://www.neuro.polymtl.ca/tips_and_tricks/how_to_use_ants and ITK code
% in ComputeOffset() itkMatrixOffsetTransformBase.hxx

offset=zeros(3,1);
translation=afftransform(end-2:end);

mat=reshape(afftransform,[3,4]);

for i=1:3
    offset(i)=translation(i)+fixed(i);
    for j=1:3
       offset(i)=offset(i)-(mat(j,i) * fixed(j)); 
    end
end
 % convert RAS to LPS (ITK uses RAS)
offset(3)=-offset(3);
mat(:,4)=offset;

mat=[mat;[0,0,0,1]];

 % convert RAS to LPS (ITK uses RAS)
mat=mat.*...
    [1 1 -1 1
    1 1 -1 1
    -1 -1 1 1
    1 1 1 1];
    %% original code in itkMatrixOffsetTransformBase
% {
%   const MatrixType & matrix = this->GetMatrix();
%   
%   OffsetType offset;
%   for(unsigned int i=0; i<NOutputDimensions; i++)
%     {
%     offset[i] = m_Translation[i] + m_Center[i];
%     for(unsigned int j=0; j<NInputDimensions; j++)
%       {
%       offset[i] -= matrix[i][j] * m_Center[j];
%       }
%     }
% 
%   m_Offset = offset;
% }