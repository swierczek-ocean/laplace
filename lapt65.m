function F = lapt65(s,a)
F = (s - sqrt(s.^2-a^2)).^4./sqrt(s.^2-a^2);
end

