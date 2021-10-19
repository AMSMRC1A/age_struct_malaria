close all
clear all
clc

x = linspace(0,1,100);
figure_setups;
plot(x,x./(5*x+1))