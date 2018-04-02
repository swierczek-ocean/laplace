a=1;
b=2;
eps=0.1;

t = 0.1:pi/10:25;

X = master_inverse_laplace_fcn(t,a,b,84,eps);


plot(t,X)


X