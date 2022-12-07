clc
clear
close all

load iddata1 z1
nx = 4;
sys = n4sid(z1, nx)
compare(z1, sys)
box on
