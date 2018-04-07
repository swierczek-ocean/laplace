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
a = 10;
b = 2;
eps = 0.5;
Sings = makesings(a,b);
t = 1.5:0.75:75;
ub = 100;
sw = 2;
start = 1;
num_test = 100;
error9 = zeros(num_test,2);
%%

%% tests
parfor jj=1:num_test
    ll = jj+start;
    True = master_inverse_laplace_fcn(t,a,b,ll,eps);
    fun = @(x)master_laplace_fcn(x,a,b,ll,eps);
    NAB = nabilt(fun,t,ub,ll,sw,Sings);
    error9(jj,:) = [LaplacePlot(True,NAB,t,ll,ub),ll];
end
%%

%% error summary
fprintf('max relative error = %g percent\n',100*max(error9(:,1)))
fprintf('mean relative error = %g percent\n',100*mean(error9(:,1)))
%%

save error

toc()