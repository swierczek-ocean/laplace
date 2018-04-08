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
a = 5;
b = 2;
eps = 0.5;
Sings = makesings(a,b);
t = 1.5:pi/20:10;
ub = 100;
sw = 2;
start = 125;
num_test = 1;
error = zeros(num_test,2);
%%

%% tests
for jj=1:num_test
    ll = jj+start-1;
    True = master_inverse_laplace_fcn(t,a,b,ll,eps);
    fun = @(x)master_laplace_fcn(x,a,b,ll,eps);
    NAB = nabilt(fun,t,ub,ll,sw,Sings);
    error(jj,:) = [LaplacePlot(True,NAB,t,ll,ub),ll];
end
%%

%% error summary
fprintf('max relative error = %g percent\n',100*max(error(:,1)))
fprintf('mean relative error = %g percent\n',100*mean(error(:,1)))
%%

error

toc()