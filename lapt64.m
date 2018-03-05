function F = lapt64(s,a)
F = (s - sqrt(s.^2-a^2)).^3./sqrt(s.^2-a^2);
end

