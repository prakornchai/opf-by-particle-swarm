clear;clc
% Optimal Power Flow
format short
Lf_ieee30bus_bygauss
%initial value
%Voltage of Gen
busdata(1,3)=1.05;
busdata(2,3)=1.0462;
busdata(5,3)=1.0459;
busdata(8,3)=1.04165;
busdata(11,3)=0.95231;
busdata(13,3)=1.05;
%Transformer Tap
linedata(11,6)=1.01;
linedata(12,6)=0.98;
linedata(15,6)=1.01;
linedata(36,6)=1.02;
%Power generator
busdata(2,7)=49.209;
busdata(5,7)=21.5135;
busdata(8,7)=22.648;
busdata(11,7)=10.4146;
busdata(13,7)=12;
Lfybus;                            % form the bus admittance matrix
% Lfnewton;                % Load flow solution 

Lfnewton_Lid
Busout;              % Prints the power flow solution on the screen
Lineflow;          % Computes and displays the line flow and losses
Lmax