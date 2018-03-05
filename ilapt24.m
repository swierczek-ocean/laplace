function f = ilapt24(T,a)
f = ((3-a^2.*T.^2).*sin(a.*T)+5*a.*T.*cos(a.*T))./(8*a);
end

