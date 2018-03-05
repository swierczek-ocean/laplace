function f = ilapt35(T,a)
f = ((8+a^2.*T.^2).*cosh(a.*T)+7*a.*T.*sinh(a.*T))./8;
end

