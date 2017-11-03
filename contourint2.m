function Output = contourint2(minx,maxx,miny,maxy,dx,fun)
tic();

sum = 0;

X1 = maxx:-dx:minx;
X1 = X1+1i*maxy;
X2 = maxy:-dx:miny;
X2 = 1i.*X2+minx;
X3 = minx:dx:maxx;
X3 = X3+1i*miny;
X4 = miny:dx:maxy;
X4 = +1i.*X4+maxx;

X = [X1,X2,X3,X4];
% X = transpose(X)
% X = transpose(X);

[n,m] = size(X);

for i=1:m-1
   deltax = X(1,i+1)-X(1,i);
   z1 = X(i);
   z2 = X(i+1);
   sum = sum + 0.5*deltax*(fun(z1) + fun(z2));    
end
   deltax = X(1,1)-X(1,m);
   z1 = X(1);
   z2 = X(n);
   Output = sum + 0.5*deltax*(fun(z1) + fun(z2)); 
   
toc()
end

