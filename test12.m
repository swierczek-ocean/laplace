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
global a, global b, global eps
a = 3;
b = 7;
eps = 0.5;
Sings = makesings(a,b);
% t = 1.1:0.125:7.1;
t = 1.15:0.5:25.15;
length = size(t,2);
ub = 100;
sw = 5;
start = 1;
num_test = 100;
errorW = zeros(num_test,length+1);
Weeks = zeros(1,length);
%%

%% tests
for jj=1:num_test
    global ll
    ll = jj+start-1;
    True = master_inverse_laplace_fcn(t,a,b,ll,eps);
    fun2 = 'laplacefcn(s)';
    % Weeks = WeeksMethod(fun,t,0.0000001,sw);
    for kk=1:length
        Weeks(kk) = wfnWeeksCoreSigmab(fun2,t(kk),512,0.5,Sings(ll)+20);
        % Weeks(kk) = WeeksMethod(fun2,t(kk),0.0000001,sw);
    end
    errorW(jj,:) = [LaplacePlot2(True,real(Weeks),t,ll),ll];
end
%%

%% error summary
fprintf('max relative error = %g percent\n',100*max(reshape(errorW(:,1:end-1),num_test*length,1)))
fprintf('mean relative error = %g percent\n',100*mean(reshape(errorW(:,1:end-1),num_test*length,1)))
%%

[mean(errorW(:,1:end),2),errorW(:,end)]

toc()