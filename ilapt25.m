function f = ilapt25(T,a)
f = ((8-a^2.*T.^2).*cos(a.*T)-7*a.*T.*sin(a.*T))./8;
end

