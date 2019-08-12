function FT = find_cost(p)
Lf_ieee30bus_bygauss1;
%Power generator
busdata(2,7)=p(1);
busdata(5,7)=p(2);
busdata(8,7)=p(3);
busdata(11,7)=p(4);
busdata(13,7)=p(5);
%Voltage of Gen
busdata(1,3)=p(6);
busdata(2,3)=p(7);
busdata(5,3)=p(8);
busdata(8,3)=p(9);
busdata(11,3)=p(10);
busdata(13,3)=p(11);
%Transformer Tap
linedata(11,6)=p(12);
linedata(12,6)=p(13);
linedata(15,6)=p(14);
linedata(36,6)=p(15);


Lfybus;                            % form the bus admittance matrix
ss=3;
switch ss
    case 1
            lfgauss;                % Load flow solution by Gauss-Seidel method
    case 2
            decouple              % Load flow solution by fast decoupled method
    case 3
            Lfnewton             % Load flow solution by Newton-Raphson method
end
%busout;              % Prints the power flow solution on the screen
Lineflow;          % Computes and displays the line flow and losses
pp=[Pg(1);Pg(2);Pg(5);Pg(8);Pg(11);Pg(13)];

%find cost
F = zeros(6,1);
c = [0;0;0;0;0;0];
b = [2;1.75;1;3.25;3;3];
a = [0.00375;0.0175;0.0625;0.00834;0.02500;0.02500];
for i=1:6,
    F(i)=(a(i)*pp(i)^2)+b(i)*pp(i)+c(i);
end
FT = sum(F);
% VAR constrain
ppp=[Pg(1);Qg(1);Qg(2);Qg(5);Qg(8);Qg(11);Qg(13)];
if ppp(1)<50|ppp(1)>200|ppp(2)<-20|ppp(2)>200|ppp(3)<-20|ppp(3)>100|ppp(4)<-15|ppp(4)>80|ppp(5)<-15|ppp(5)>60|ppp(6)<-10|ppp(6)>50|ppp(7)<-15|ppp(7)>60
    FT=FT+10000;
end
%Line flow security

%line loading constrain
%Trm_loss=Transmissionloss;
%if Trm_loss>100
    %FT=FT+50*(Trm_loss-100);
%end

MVA_line=S_line;
if MVA_line(1)>130|MVA_line(2)>130|MVA_line(4)>65|MVA_line(8)>130|MVA_line(5)>130|MVA_line(6)>65|MVA_line(11)>90|...
   MVA_line(14)>70|MVA_line(17)>130|MVA_line(18)>32|MVA_line(19)>65|MVA_line(20)>32|MVA_line(27)>65|MVA_line(28)>65|MVA_line(12)>65|...
   MVA_line(37)>65|MVA_line(38)>32|MVA_line(39)>32|MVA_line(40)>32|MVA_line(43)>16|MVA_line(49)>16|MVA_line(46)>16|MVA_line(53)>16|...
   MVA_line(55)>32|MVA_line(31)>32|MVA_line(32)>32|MVA_line(33)>32|MVA_line(34)>32|MVA_line(59)>32|MVA_line(47)>16|MVA_line(62)>16|...
   MVA_line(64)>16|MVA_line(67)>16|MVA_line(69)>16|MVA_line(70)>16|MVA_line(76)>65|MVA_line(74)>16|MVA_line(75)>16|MVA_line(80)>16|MVA_line(25)>32|MVA_line(21)>32
                              FT=FT+10000;
                             % disp('Over LineFlow')
end
Vm_max=1.1;
Vm_min=0.9;
kv=10000;
for kk=1:15
    if Vm(kk)>Vm_max
        FT=FT+kv*(Vm(kk)-Vm_max)^2;
    elseif Vm(kk)<Vm_min
        FT=FT+kv*(Vm_min-Vm(kk))^2;
    end
end
