%% preliminaries
% clc
close all
clear
format long 
tic()
set(0,'DefaultFigureVisible','off')
colors
%%

%% setting up time values and shifts
a = 1;
b = 2;
t = 3;
[n,m] = size(t);
ub = 5;
start = 11;
num_test = 20;
error = zeros(num_test-start+1,1);
%%

%% tests
for jj=start:num_test
    True = master_inverse_laplace_fcn(t,a,b,jj);
    fun = @(x)master_laplace_fcn(x,a,b,jj);
    NAB = nabilt(fun,t,ub);
    error(jj,1) = LaplacePlot(True,NAB,t,jj,ub);
end
%%

%% error summary
% fprintf('max relative error = %g\n',max(error))
% fprintf('mean relative error = %g\n',mean(error))
fprintf('max relative error = %g percent\n',100*max(error))
fprintf('mean relative error = %g percent\n',100*mean(error))
%%

toc()