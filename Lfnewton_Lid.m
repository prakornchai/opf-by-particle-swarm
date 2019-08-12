%   Power flow solution by Newton-Raphson method
%   Copyright (c) 1998 by  H. Saadat
%   Revision 1 (Aug. 99) To include two or more parallel lines


ns=0; ng=0; Vm=0; delta=0; yload=0; deltad=0;
nbus = length(busdata(:,1));
kb=[];Vm=[]; delta=[]; Pd=[]; Qd=[]; Pg=[]; Qg=[]; Qmin=[]; Qmax=[];  % Added (6-8-00)
Pk=[]; P=[]; Qk=[]; Q=[]; S=[]; V=[]; % Added (6-8-00)
for k=1:nbus
n=busdata(k,1);
kb(n)=busdata(k,2); Vm(n)=busdata(k,3); delta(n)=busdata(k, 4);
Pd(n)=busdata(k,5); Qd(n)=busdata(k,6); Pg(n)=busdata(k,7); Qg(n) = busdata(k,8);
Qmin(n)=busdata(k, 9); Qmax(n)=busdata(k, 10);
Qsh(n)=busdata(k, 11);
    if Vm(n) <= 0  Vm(n) = 1.0; V(n) = 1 + j*0;
    else delta(n) = pi/180*delta(n);
         V(n) = Vm(n)*(cos(delta(n)) + j*sin(delta(n)));
         P(n)=(Pg(n)-Pd(n))/basemva;
         Q(n)=(Qg(n)-Qd(n)+ Qsh(n))/basemva;
         S(n) = P(n) + j*Q(n);
    end
end
for k=1:nbus
if kb(k) == 1, ns = ns+1; else, end
if kb(k) == 2 ng = ng+1; else, end
ngs(k) = ng;
nss(k) = ns;
end
Ym=abs(Ybus); t = angle(Ybus);
m=2*nbus-ng-2*ns;
maxerror = 1; converge=1;
iter = 0;
%%%% added for parallel lines (Aug. 99)
mline=ones(nbr,1);
for k=1:nbr
      for m=k+1:nbr
         if((nl(k)==nl(m)) & (nr(k)==nr(m)));
            mline(m)=2;
         elseif ((nl(k)==nr(m)) & (nr(k)==nl(m)));
         mline(m)=2;
         else, end
      end
   end
%%%   end of statements for parallel lines (Aug. 99)   

% Start of iterations
clear A  DC   J  DX
    while maxerror >= accuracy & iter <= maxiter % Test for max. power mismatch
for ii=1:m
for k=1:m
   A(ii,k)=0;      %Initializing Jacobian matrix
end, end
iter = iter+1;
for n=1:nbus
nn=n-nss(n);
lm=nbus+n-ngs(n)-nss(n)-ns;
J11=0; J22=0; J33=0; J44=0;
   for ii=1:nbr
   	if mline(ii)==1   % Added to include parallel lines (Aug. 99)
      	if nl(ii) == n | nr(ii) == n
         	if nl(ii) == n ,  l = nr(ii); end
         	if nr(ii) == n , l = nl(ii); end
         J11=J11+ Vm(n)*Vm(l)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));
         J33=J33+ Vm(n)*Vm(l)*Ym(n,l)*cos(t(n,l)- delta(n) + delta(l));
        		if kb(n)~=1
        		J22=J22+ Vm(l)*Ym(n,l)*cos(t(n,l)- delta(n) + delta(l));
        		J44=J44+ Vm(l)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));
        		else, end
        		if kb(n) ~= 1  & kb(l) ~=1
        			lk = nbus+l-ngs(l)-nss(l)-ns;
        			ll = l -nss(l);
      			% off diagonalelements of J1
        			A(nn, ll) =-Vm(n)*Vm(l)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));
              		if kb(l) == 0  % off diagonal elements of J2
              		A(nn, lk) =Vm(n)*Ym(n,l)*cos(t(n,l)- delta(n) + delta(l));end
              		if kb(n) == 0  % off diagonal elements of J3
              		A(lm, ll) =-Vm(n)*Vm(l)*Ym(n,l)*cos(t(n,l)- delta(n)+delta(l)); end
              		if kb(n) == 0 & kb(l) == 0  % off diagonal elements of  J4
              		A(lm, lk) =-Vm(n)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));end
        	   else end
     		else , end
      else, end   
   end
   Pk = Vm(n)^2*Ym(n,n)*cos(t(n,n))+J33;
   Qk = -Vm(n)^2*Ym(n,n)*sin(t(n,n))-J11;
   if kb(n) == 1 P(n)=Pk; Q(n) = Qk; end   % Swing bus P
     if kb(n) == 2  Q(n)=Qk;
         if Qmax(n) ~= 0
           Qgc = Q(n)*basemva + Qd(n) - Qsh(n);
           if iter <= 7                  % Between the 2th & 6th iterations
              if iter > 2                % the Mvar of generator buses are
                if Qgc  < Qmin(n),       % tested. If not within limits Vm(n)
                Vm(n) = Vm(n) + 0.01;    % is changed in steps of 0.01 pu to
                elseif Qgc  > Qmax(n),   % bring the generator Mvar within
                Vm(n) = Vm(n) - 0.01;end % the specified limits.
              else, end
           else,end
         else,end
     end
   if kb(n) ~= 1
     A(nn,nn) = J11;  %diagonal elements of J1
     DC(nn) = P(n)-Pk;
   end
   if kb(n) == 0
     A(nn,lm) = 2*Vm(n)*Ym(n,n)*cos(t(n,n))+J22;  %diagonal elements of J2
     A(lm,nn)= J33;        %diagonal elements of J3
     A(lm,lm) =-2*Vm(n)*Ym(n,n)*sin(t(n,n))-J44;  %diagonal of elements of J4
     DC(lm) = Q(n)-Qk;
   end
end
DX=A\DC';
for n=1:nbus
  nn=n-nss(n);
  lm=nbus+n-ngs(n)-nss(n)-ns;
    if kb(n) ~= 1
    delta(n) = delta(n)+DX(nn); end
    if kb(n) == 0
    Vm(n)=Vm(n)+DX(lm); end
 end
  maxerror=max(abs(DC));
     if iter == maxiter & maxerror > accuracy 
  % fprintf('\nWARNING: Iterative solution did not converged after ')
   %fprintf('%g', iter), fprintf(' iterations.\n\n')
   %fprintf('Press Enter to terminate the iterations and print the results \n')
   %converge = 0; pause, 
     else, end
   
end

if converge ~= 1
   tech= ('                      ITERATIVE SOLUTION DID NOT CONVERGE'); else, 
   tech=('                   Power Flow Solution by Newton-Raphson Method');
end   
V = Vm.*cos(delta)+j*Vm.*sin(delta);
deltad=180/pi*delta;
i=sqrt(-1);
k=0;
for n = 1:nbus
     if kb(n) == 1
     k=k+1;
     S(n)= P(n)+j*Q(n);
     Pg(n) = P(n)*basemva + Pd(n);
     Qg(n) = Q(n)*basemva + Qd(n) - Qsh(n);
     Pgg(k)=Pg(n);
     Qgg(k)=Qg(n);     %june 97
     elseif  kb(n) ==2
     k=k+1;
     S(n)=P(n)+j*Q(n);
     Qg(n) = Q(n)*basemva + Qd(n) - Qsh(n);
     Pgg(k)=Pg(n);
     Qgg(k)=Qg(n);  % June 1997
  end
yload(n) = (Pd(n)- j*Qd(n)+j*Qsh(n))/(basemva*Vm(n)^2);
end
busdata(:,3)=Vm'; busdata(:,4)=deltad';
Pgt = sum(Pg);  Qgt = sum(Qg); Pdt = sum(Pd); Qdt = sum(Qd); Qsht = sum(Qsh);

%clear A DC DX  J11 J22 J33 J44 Qk delta lk ll lm
%clear A DC DX  J11 J22 J33  Qk delta lk ll lm
linedt= linedata;                     % calling the line data for the system
busdt= busdata;                    % calling the bus data for the system
nb= max(busdt(:,1));   % gives the total number of buses in the system
type =busdt(:,2); % identify the type of bus i.e. ref., generator, and load
pv = find(type==2 | type==1); % identify the PV bus  
npv = length(pv); % gives the number of PV buses
pq = find(type==0);                         % identify the PQ bus 
npq = length(pq);                              % gives the number of PQ buses

for m=1:npq
   for n= 1:npq
      YLL (m,n) = Ybus(pq(m), pq(n));
   end
end
for m=1:npq,
  for n = 1:npv,
       YLG(m,n)= Ybus(pq(m), pv(n));
  end
end
F = -(inv(YLL))*YLG;
for jj=1:npq
    Ld=0;
 for i=1:npv
    Ld=Ld+F(jj,i)*V(pv(i)); 
 end   
    L(jj)=abs(1-Ld/V(pq(jj)));
end
Lmax=max(L);











