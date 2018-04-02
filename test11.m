%% preliminaries
% clc
close all
clear
format long 
tic()
set(0,'DefaultFigureVisible','on')
colors
%%

%% setting up time values and shifts
a = 10;
b = 2;
eps = 0.5;
Sings = makesings(a,b);
t = 1.5:pi/20:10;
ub = 1000;
sw = 2;
start = 97;
num_test = 97;
error = zeros(num_test-start+1,2);
%%

%% tests
for jj=start:num_test
    True = master_inverse_laplace_fcn(t,a,b,jj,eps);
    fun = @(x)master_laplace_fcn(x,a,b,jj,eps);
    NAB = nabilt(fun,t,ub,jj,sw,Sings);
    error(jj-start+1,:) = [LaplacePlot(True,NAB,t,jj,ub),jj];
end
%%

%% error summary
fprintf('max relative error = %g percent\n',100*max(error(:,1)))
fprintf('mean relative error = %g percent\n',100*mean(error(:,1)))
%%

error

toc()