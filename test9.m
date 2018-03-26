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
t = 4;
ub = 100;
start = 89;
num_test = 90;
error = zeros(num_test-start+1,2);
counter = 1;
%%

%% tests
for jj=start:num_test
    True = master_inverse_laplace_fcn(t,a,b,jj);
    fun = @(x)master_laplace_fcn(x,a,b,jj);
    NAB = nabilt(fun,t,ub);
    error(counter,:) = [LaplacePlot(True,NAB,t,jj,ub),jj];
    counter = counter + 1;
end
%%

%% error summary
fprintf('max relative error = %g percent\n',100*max(error(:,1)))
fprintf('mean relative error = %g percent\n',100*mean(error(:,1)))
%%

error

toc()