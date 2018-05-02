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
global a, global b, global eps
a = 3;
b = 1;
eps = 0.5;
Sings = makesings(a,b);
t = 1.15:0.25:15.15;
timelength = size(t,2);
ub = 100;
sw = 1;
start = 3;
num_test = 1;
errorW = zeros(num_test,timelength+1);
Weeks = zeros(1,timelength);
RET = 0.02*ones(1,timelength);
%%

%% tests
for jj=1:num_test
    global ll
    ll = jj+start-1;
    True = master_inverse_laplace_fcn(t,a,b,ll,eps);
    fun2 = 'laplacefcn(s)';
    Weeks = WeeksMethod(fun2,t,RET,sw);
    errorW(jj,:) = [LaplacePlot2(True,real(Weeks),t,ll,a,b,sw),ll];
end
%%

%% error summary
fprintf('max relative error = %g percent\n',100*max(reshape(errorW(:,1:end-1),num_test*timelength,1)))
fprintf('mean relative error = %g percent\n',100*mean(reshape(errorW(:,1:end-1),num_test*timelength,1)))
%%

[mean(errorW(:,1:end),2),errorW(:,end)]

save errorW

toc()