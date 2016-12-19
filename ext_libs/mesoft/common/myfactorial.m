function res= myfactorial(in)

res= in;
for i= 1:numel(in)
    res(i)= factorial(in(i));    
end