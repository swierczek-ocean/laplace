function  mobius2(n,b,sigma)
tic();

X = linspace(0,2*pi,n);
X = X(2:end-1);

W = exp(1i.*X);

figure
plot(W,'o','MarkerSize',6,'Color','red')
axis([-1.1 1.2 -1.1 1.2])
title('Uniform Points from Unit Circle')
xlabel('Real Part')
ylabel('Imaginary Part')
print('UnitCircle','-djpeg')

Z = (sigma-b).*W-sigma-b;
Z = Z./(W-1);

figure
plot(Z,'o','MarkerSize',6,'Color','red')
title(['Inverse Mobius Transform with b = ', num2str(b),' and sigma = ', num2str(sigma)])
xlabel('Real Part')
ylabel('Imaginary Part')
print(['InverseMobius',num2str(b)],'-djpeg')

toc()
end

