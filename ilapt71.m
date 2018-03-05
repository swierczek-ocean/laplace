function f = ilapt71(T,a)
f = besseli(0,a.*T)+a.*T.*besseli(1,a.*T);
end

