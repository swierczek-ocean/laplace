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
a = 6;
b = 2;
eps = 0.5;
Sings = makesings(a,b);
t = 1.1:0.25:25.1;
length = size(t,2);
ub = 100;
sw = 2;
start = 1;
num_test = 100;
error1 = zeros(num_test,length+1);
%%

%% tests
parfor jj=1:num_test
    ll = jj+start-1;
    True = master_inverse_laplace_fcn(t,a,b,ll,eps);
    fun = @(x)master_laplace_fcn(x,a,b,ll,eps);
    NAB = nabilt(fun,t,ub,ll,sw,Sings);
    error1(jj,:) = [LaplacePlot(True,NAB,t,ll,ub,a,b),ll];
end
%%

%% error summary
fprintf('max relative error = %g percent\n',100*max(reshape(error1(:,1:end-1),num_test*length,1)))
fprintf('mean relative error = %g percent\n',100*mean(reshape(error1(:,1:end-1),num_test*length,1)))
%%

error_summary1 = [mean(error1(:,1:end),2),error1(:,end)];

save error1
save error_summary1

toc()