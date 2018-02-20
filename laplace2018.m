clc
close all
clear

t = 1:100;

fun = @(s)1/(s-9)^2;

fun3 = @(t)t.*exp(9.*t);

ans = fun3(t);

f = nabilt(fun,t);

plot(t,f)
hold on
plot(t,ans)
hold off





