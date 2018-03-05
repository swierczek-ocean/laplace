function f = ilapt48(T,a)
f = (sin(a.*T).*cosh(a.*T)+cos(a.*T).*sinh(a.*T))./(2*a);
end

