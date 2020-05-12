function fibs=ea_discfibers_addjitter(fibs,jit)

for fib=1:length(fibs)
    fibs{fib}=fibs{fib}+randn(size(fibs{fib}))*jit;
end

