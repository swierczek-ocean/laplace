function f = ilapt15(T,a)
f = (a.*T.*cosh(a.*T) - sinh(a.*T))./(2*a^3);
end

