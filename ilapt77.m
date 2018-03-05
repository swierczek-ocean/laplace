function f = ilapt77(T,a)
f = (T./a).*besselj(1,2.*sqrt(a.*T));
end

