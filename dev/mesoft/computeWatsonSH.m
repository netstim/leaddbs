%%

%syms k t 
%/ ((3.1416^(1/2)*erfi((k)^(1/2)))/(2*(k)^(1/2))))
%int(exp(k*t^2)*orthpoly::legendre(12,1.0*t),t,0,1) 

for l = 2:2:14,

    x = evalin(symengine,sprintf('simplify(int(exp(k*t^2)*orthpoly::legendre(%i,1.0*t),t=0..1) / ((PI^(1/2)*erfi((k)^(1/2)))/(2*(k)^(1/2))));',l))



    y = strrep(char(x),'erfi','tmp1');
    y = strrep(y,'erfc','tmp2');
    y = strrep(y,'erf','myerf');
    y = strrep(y,'tmp1','erfi');
    y = strrep(y,'tmp2','erfc');

    y = strrep(strrep(strrep(y,'*','.*'),'^','.^'),'/','./');

    f = @(k) real(eval(y));

    k = 0.1:0.1:200;
    plot(k,(f(k)));
    drawnow;
    l
    lut(l/2,:) = f(k);
end;

watsonlut.lut = lut;
watsonlut.k = k;

save watsonlut -struct watsonlut

    