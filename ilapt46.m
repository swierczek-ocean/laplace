function f = ilapt46(T,a)
f = (sin(a.*T).*cosh(a.*T)-cos(a.*T).*sinh(a.*T))./(4*a^3);
end

