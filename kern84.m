function F = kern84(u,a,t,n)
F = (u.^n).*exp(-u.^2./(4*a^2.*t)).*besselj(2*n,2.*sqrt(u));
end

