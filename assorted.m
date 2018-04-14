a=3;
b=5;
eps=0.1;



[X,Y] = meshgrid(-10:0.1:10,-10:0.1:10);
Z = imag(master_laplace_fcn(X+1i.*Y,a,b,233,eps));
surf(X,Y,Z)
xlabel('x')
ylabel('iy')
zlabel('real part')