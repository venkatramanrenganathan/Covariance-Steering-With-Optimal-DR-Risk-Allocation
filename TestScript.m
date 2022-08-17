clear all; close all; clc;
N = 1000;
PS = loadProblemSetup(N);
V = zeros(PS.N*PS.nu,1);
K = zeros(PS.N*PS.nu,(PS.N+1)*PS.nx);
PostProcess(1,10,V,K,PS);