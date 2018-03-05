function f = ilapt68(T,a)
f = besselj(0,a.*T)-a.*T.*besselj(1,a.*T);
end

