function f = ilapt45(T,a)
f = (exp(a.*T)+2.*exp(-a.*T./2).*cos(sqrt(3)*a.*T./2))./3;
end

