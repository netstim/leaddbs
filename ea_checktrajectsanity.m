function res=checktrajectsanity(trajvector)

load('trajvectors');
trajvectors=[trajvectors;trajvector];
save('trajvectors','trajvectors');
res=1;


