clc
close all
clear
format long 

t = 1:3;

fun = @(s)1./(s).^2;

fun6 = @(t)(t>0);

tf = fun6(t)

f = nabilt(fun,t)

plot(t,f)
hold on
plot(t,tf)
hold off





