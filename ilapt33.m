function f = ilapt33(T,a)
f = (3.*T.*sinh(a.*T)+a.*T.^2.*cosh(a.*T))./(8*a);
end

