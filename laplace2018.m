clc
close all
clear


fun = @(x)1/(x^2-9)^2;

syms x

poles(fun,x)

fun2 = @(x)1/fun(x);

fzero(fun2,2.999)

vpasolve(fun2(x) == 0, x)














