function mobius(m,n,radius,b,sd,R)
tic();

X = unifrnd(-10,10,n,1) + 1i.*unifrnd(-10,10,n,1);
sigma1 = max(real(X));
sigma = sigma1 + radius;

Contour1 = sigma + 1i.*normrnd(0,sd,m,1);
Contour2 = sigma + 1i.*unifrnd(-R,R,m,1);

MX = (X-sigma-b)./(X-sigma+b);
MC1 = (Contour1-sigma-b)./(Contour1-sigma+b);
MC2 = (Contour2-sigma-b)./(Contour2-sigma+b);

figure
plot(X,'o','MarkerSize',5,'Color','red')
hold on
plot(Contour1,'o','MarkerSize',5,'Color','blue')
plot(Contour2,'o','MarkerSize',5,'Color','green')

figure
plot(MX,'o','MarkerSize',5,'Color','red')
hold on
plot(MC1,'o','MarkerSize',5,'Color','blue')
plot(MC2,'o','MarkerSize',5,'Color','green')
axis([-1.1 1.2 -1.1 1.2])

toc()
end

