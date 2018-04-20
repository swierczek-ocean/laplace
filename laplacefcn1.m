function F = laplacefcn1(s)
global a, global b, global eps, global ll, global hshift
F = master_laplace_fcn(s+hshift,a,b,ll,eps);
end

