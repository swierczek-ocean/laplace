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
t = 1.5:9.5;
ub = 100;
start = 1;
num_test = 100;
error = zeros(num_test-start+1,2);
%counter = 1;
%%

%% tests
parfor jj=start:num_test
    True = master_inverse_laplace_fcn(t,a,b,jj,eps);
    fun = @(x)master_laplace_fcn(x,a,b,jj,eps);
    NAB = nabilt(fun,t,ub);
    error(jj,:) = [LaplacePlot(True,NAB,t,jj,ub),jj];
    %counter = counter + 1;
end
%%

%% error summary
fprintf('max relative error = %g percent\n',100*max(error(:,1)))
fprintf('mean relative error = %g percent\n',100*mean(error(:,1)))
%%

save error

toc()