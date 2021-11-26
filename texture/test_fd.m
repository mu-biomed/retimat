%% Test FD
close all;clc;clearvars;

I = imread('../data/D33.gif');
FD = compute_FD(I, true);
