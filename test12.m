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
a = 1;
b = 1;
eps = 0.5;
Sings = makesings(a,b);
t = 1.15:0.5:15.15;
length = size(t,2);
ub = 100;
sw = 5;
start = 1;
num_test = 5;
errorW = zeros(num_test,length+1);
Weeks = zeros(1,length);
%%

%% tests
for jj=1:num_test
    global ll
    ll = jj+start-1;
    True = master_inverse_laplace_fcn(t,a,b,ll,eps);
    fun2 = 'laplacefcn(s)';
%    Weeks = WeeksMethod(fun2,t,0.0000001,sw);
    for kk=1:length
        Weeks(kk) = wfnWeeksCoreSigmab(fun2,t(kk),512,0.5,Sings(ll)+1/t(kk));
        % Weeks(kk) = WeeksMethod(fun2,t(kk),0.0000001,sw);
    end
    errorW(jj,:) = [LaplacePlot2(True,real(Weeks),t,ll,a,b,sw),ll];
end
%%

%% error summary
fprintf('max relative error = %g percent\n',100*max(reshape(errorW(:,1:end-1),num_test*length,1)))
fprintf('mean relative error = %g percent\n',100*mean(reshape(errorW(:,1:end-1),num_test*length,1)))
%%

[mean(errorW(:,1:end),2),errorW(:,end)]

toc()