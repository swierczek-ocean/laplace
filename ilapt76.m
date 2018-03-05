function f = ilapt76(T,a)
f = sqrt(T./a).*besselj(0,2.*sqrt(a.*T));
end

