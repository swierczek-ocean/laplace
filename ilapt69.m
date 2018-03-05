function f = ilapt69(T,a)
f = T.*besseli(1,a.*T)./a;
end

