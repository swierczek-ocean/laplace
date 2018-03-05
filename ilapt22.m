function f = ilapt22(T,a)
f = ((1+a^2.*T.^2).*sin(a.*T) - a.*T.*cos(a.*T))./(8*a^3);
end


