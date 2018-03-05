function f = ilapt20(T,a)
f = ((3-a^2.*T.^2).*sin(a.*T)-3*a.*T.*cos(a.*T))./(8*a^5);
end

