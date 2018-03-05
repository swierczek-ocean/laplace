function f = ilapt21(T,a)
f = (T.*sin(a.*T)-a.*T.^2.*cos(a.*T))./(8*a^3);
end

