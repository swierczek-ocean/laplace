
%% preliminaries
clc
close all
clear
format long 
tic()
set(0,'DefaultFigureVisible','on')
colors
%%

%% setting up time values and shifts
a = 1;
t = 3;
[n,m] = size(t);
%%

%% test 1
True = ilapt1(t,a);
fun = @(x)lapt1(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,1)
%%

%% test 2
True = ilapt2(t,a);
fun = @(x)lapt2(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,2)
%%

%% test 3
True = ilapt3(t,a);
fun = @(x)lapt3(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,3)
%%

%% test 4
True = ilapt4(t,a);
fun = @(x)lapt4(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,4)
%%

%% test 5
True = ilapt5(t,a);
fun = @(x)lapt5(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,5)
%%

%% test 6
True = ilapt6(t,a);
fun = @(x)lapt6(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,6)
%%

%% test 7
True = ilapt7(t,a);
fun = @(x)lapt7(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,7)
%%

%% test 8
True = ilapt8(t,a);
fun = @(x)lapt8(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,8)
%%

%% test 9
True = ilapt9(t,a);
fun = @(x)lapt9(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,9)
%%

%% test 10
True = ilapt10(t,a);
fun = @(x)lapt10(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,10)
%%

%% test 11
True = ilapt11(t,a);
fun = @(x)lapt11(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,11)
%%

%% test 12
True = ilapt12(t,a);
fun = @(x)lapt12(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,12)
%%

%% test 13
True = ilapt13(t,a);
fun = @(x)lapt13(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,13)
%%

%% test 14
True = ilapt14(t,a);
fun = @(x)lapt14(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,14)
%%

%% test 15
True = ilapt15(t,a);
fun = @(x)lapt15(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,15)
%%

%% test 16
True = ilapt16(t,a);
fun = @(x)lapt16(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,16)
%%

%% test 17
True = ilapt17(t,a);
fun = @(x)lapt17(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,17)
%%

%% test 18
True = ilapt18(t,a);
fun = @(x)lapt18(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,18)
%%

%% test 19
True = ilapt19(t,a);
fun = @(x)lapt19(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,19)
%%

%% test 20
True = ilapt20(t,a);
fun = @(x)lapt20(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,20)
%%

%% test 21
True = ilapt21(t,a);
fun = @(x)lapt21(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,21)
%%

%% test 22
True = ilapt22(t,a);
fun = @(x)lapt22(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,22)
%%

%% test 23
True = ilapt23(t,a);
fun = @(x)lapt23(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,23)
%%

%% test 24
True = ilapt24(t,a);
fun = @(x)lapt24(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,24)
%%

%% test 25
True = ilapt25(t,a);
fun = @(x)lapt25(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,25)
%%

%% test 26
True = ilapt26(t,a);
fun = @(x)lapt26(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,26)
%%

%% test 27
True = ilapt27(t,a);
fun = @(x)lapt27(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,27)
%%

%% test 28
True = ilapt28(t,a);
fun = @(x)lapt28(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,28)
%%

%% test 29
True = ilapt29(t,a);
fun = @(x)lapt29(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,29)
%%

%% test 30
True = ilapt30(t,a);
fun = @(x)lapt30(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,30)
%%

%% test 31
True = ilapt31(t,a);
fun = @(x)lapt31(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,31)
%%

%% test 32
True = ilapt32(t,a);
fun = @(x)lapt32(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,32)
%%

%% test 33
True = ilapt33(t,a);
fun = @(x)lapt33(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,33)
%%

%% test 34
True = ilapt34(t,a);
fun = @(x)lapt34(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,34)
%%

%% test 35
True = ilapt35(t,a);
fun = @(x)lapt35(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,35)
%%

%% test 36
True = ilapt36(t,a);
fun = @(x)lapt36(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,36)
%%

%% test 37
True = ilapt37(t,a);
fun = @(x)lapt37(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,37)
%%

%% test 38
True = ilapt38(t,a);
fun = @(x)lapt38(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,38)
%%

%% test 39
True = ilapt39(t,a);
fun = @(x)lapt39(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,39)
%%

%% test 40
True = ilapt40(t,a);
fun = @(x)lapt40(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,40)
%%

%% test 41
True = ilapt41(t,a);
fun = @(x)lapt41(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,41)
%%

%% test 42
True = ilapt42(t,a);
fun = @(x)lapt42(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,42)
%%

%% test 43
True = ilapt43(t,a);
fun = @(x)lapt43(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,43)
%%

%% test 44
True = ilapt44(t,a);
fun = @(x)lapt44(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,44)
%%

%% test 45
True = ilapt45(t,a);
fun = @(x)lapt45(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,45)
%%

%% test 46
True = ilapt46(t,a);
fun = @(x)lapt46(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,46)
%%

%% test 47
True = ilapt47(t,a);
fun = @(x)lapt47(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,47)
%%

%% test 48
True = ilapt48(t,a);
fun = @(x)lapt48(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,48)
%%

%% test 49
True = ilapt49(t,a);
fun = @(x)lapt49(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,49)
%%

%% test 50
True = ilapt50(t,a);
fun = @(x)lapt50(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,50)
%%

%% test 51
True = ilapt51(t,a);
fun = @(x)lapt51(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,51)
%%

%% test 52
True = ilapt52(t,a);
fun = @(x)lapt52(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,52)
%%

%% test 53
True = ilapt53(t,a);
fun = @(x)lapt53(x,a);
NAB = nabilt(fun,t);
LaplacePlot(True,NAB,t,53)
%%




toc()


















