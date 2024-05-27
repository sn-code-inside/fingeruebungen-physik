clc
clear 
close all 
format long g

addpath('IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

p = 10;
x = linspace(-10,10,500);
y = x.^2/2/p;
d = p/2;
y1 = x.^2/2/p + d;
figure(1)
plot(x,y);
hold on
plot(x,y1);