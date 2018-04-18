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
a = 5;
b = 3;
eps = 0.5;
Sings = makesings(a,b);
% t = 1.1:0.125:7.1;
t = 1.15:0.25:15.15;
length = size(t,2);
ub = 100;
sw = 2;
start = 111;
num_test = 42;
error = zeros(num_test,length+1);
%%

%% tests
for jj=1:num_test
    ll = jj+start-1;
    True = master_inverse_laplace_fcn(t,a,b,ll,eps);
    fun = @(x)master_laplace_fcn(x,a,b,ll,eps);
    NAB = nabilt(fun,t,ub,ll,sw,Sings);
    error(jj,:) = [LaplacePlot(True,NAB,t,ll,ub,a,b),ll];
end
%%

%% error summary
fprintf('max relative error = %g percent\n',100*max(reshape(error(:,1:end-1),num_test*length,1)))
fprintf('mean relative error = %g percent\n',100*mean(reshape(error(:,1:end-1),num_test*length,1)))
%%

[mean(error(:,1:end),2),error(:,end)]

save error

toc()