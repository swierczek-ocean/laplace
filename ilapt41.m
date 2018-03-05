function f = ilapt41(T,a)
f = exp(0.5*a.*T).*(sqrt(3).*sin(sqrt(3)*a.*T./2)+...
    cos(sqrt(3)*a.*T./2)-exp(-3*a.*T./2))./(3*a);
end

