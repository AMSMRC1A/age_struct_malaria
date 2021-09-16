clear all
close all

%%
g = @(x) x;
f = @(a) integral(g, 0,a);
integral(f, 0, 1,'ArrayValued',true)

%%
xmax = @(a) a;
integral2(@(a,x) g(x), 0,1,0,xmax)