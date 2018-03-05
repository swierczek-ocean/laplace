function f = ilapt30(T,a)
f = ((3+a^2.*T.^2).*sinh(a.*T)-3*a.*T.*cosh(a.*T))./(8*a^5);
end

