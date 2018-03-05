function f = ilapt66(T,a)
f = T.*besselj(1,a.*T)./a;
end

