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
b = 1;
eps = 0.5;
Sings = makesings(a,b);
t1 = 1.2:0.5:15.2;
t2 = 1.3:0.5:15.3;
t3 = 1.1:0.5:15.1;
t = 1.1:0.1:15.2;
timelength = size(t1,2);
ub = 100;
sw = 0;
start = 130;
num_test = 10;
NAB = zeros(1,timelength);
Weeks = zeros(1,timelength);
SWeeks = zeros(1,timelength);
RET = 0.02*ones(1,timelength);
%%

%% tests
for jj=1:num_test
    global ll
    ll = jj+start-1;
    fun = @(x)master_laplace_fcn(x,a,b,ll,eps);
    NAB = nabilt(fun,t3,ub,ll,sw,Sings);
    global hshift
    hshift = Sings(ll,1);
    True = master_inverse_laplace_fcn(t,a,b,ll,eps);
    fun2 = 'laplacefcn1(s)';
    SWeeks = WeeksMethod(fun2,t2,RET,sw);
    SWeeks = exp(hshift.*t2).*SWeeks;
    fun3 = 'laplacefcn(s)';
    Weeks = WeeksMethod(fun3,t1,RET,sw);
    LaplacePlot3(True,NAB,Weeks,SWeeks,t,t1,t2,t3,ll,a,b,eps)
end
%%

toc()