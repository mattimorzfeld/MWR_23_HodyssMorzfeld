%%
clear 
close all
clc


%% Simulation parameters
%% ------------------------------------------
dt = 0.05;
T = .2*1e5;
t = 0:dt:T;
Steps = length(t);
%% ------------------------------------------

%% Model parameters
%% ------------------------------------------
F = 8;
n = 500;
%% ------------------------------------------

%% Run the model
%% ------------------------------------------
xo = randn(n,1);
X = model(xo,dt,Steps,F);
% plot(t,yC(1,:))

FileName = strcat('LongSim_n',num2str(n),'.mat');
save(FileName,'X')
%% ------------------------------------------