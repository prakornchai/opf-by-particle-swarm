%   Fast Decoupled Power Flow Solution
%   Copyright (c)1998 by Hadi Saadat.
%   Revision 1 (Aug. 99 )  Modified to include two or more parallel lines
ns=0; Vm=0; delta=0; yload=0; deltad=0;
nbus = length(busdata(:,1));
kb=[];Vm=[]; delta=[]; Pd=[]; Qd=[]; Pg=[]; Qg=[]; Qmin=[]; Qmax=[]; % Added (6-8-00)
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
if kb(n) == 1, ns = ns+1; else, end
nss(n) = ns;
end
Ym = abs(Ybus); t = angle(Ybus);
ii=0;
for ib=1:nbus
     if kb(ib) == 0 | kb(ib) == 2
     ii = ii+1;
      jj=0;
         for jb=1:nbus
             if kb(jb) == 0 | kb(jb) == 2
             jj = jj+1;
             B1(ii,jj)=imag(Ybus(ib,jb));
             else,end
         end
     else, end
end

ii=0;
for ib=1:nbus
     if kb(ib) == 0
     ii = ii+1;
      jj=0;
         for jb=1:nbus
             if kb(jb) == 0
             jj = jj+1;
             B2(ii,jj)=imag(Ybus(ib,jb));
             else,end
         end
     else, end
end
B1inv=inv(B1); B2inv = inv(B2);

maxerror = 1; converge = 1; 
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
while maxerror >= accuracy & iter <= maxiter % Test for max. power mismatch
iter = iter+1;
id=0; iv=0;
for n=1:nbus
nn=n-nss(n);
J11=0;  J33=0;
	for ii=1:nbr
		if mline(ii)==1   % Added to include parallel lines (Aug. 99) 
     		if nl(ii) == n | nr(ii) == n
        			if nl(ii) == n,  l = nr(ii); end
        			if nr(ii) == n,  l = nl(ii); end
        		J11=J11+ Vm(n)*Vm(l)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));
        		J33=J33+ Vm(n)*Vm(l)*Ym(n,l)*cos(t(n,l)- delta(n) + delta(l));
     		else , end
      else, end
  end
   Pk = Vm(n)^2*Ym(n,n)*cos(t(n,n))+J33;
   Qk = -Vm(n)^2*Ym(n,n)*sin(t(n,n))-J11;
   if kb(n) == 1 P(n)=Pk; Q(n) = Qk; end   % Swing bus P
     if kb(n) == 2  Q(n)=Qk;
         Qgc = Q(n)*basemva + Qd(n) - Qsh(n);
         if Qmax(n) ~= 0
           if iter <= 20                 % Between the 1th & 6th iterations
              if iter >= 10              % the Mvar of generator buses are
                if Qgc  < Qmin(n),       % tested. If not within limits Vm(n)
                Vm(n) = Vm(n) + 0.005;   % is changed in steps of 0.05 pu to
                elseif Qgc > Qmax(n),    % bring the generator Mvar within
                Vm(n) = Vm(n) - 0.005;end % the specified limits.
              else, end
           else,end
         else,end
     end
   if kb(n) ~= 1
   id = id+1;
     DP(id) = P(n)-Pk;
     DPV(id) = (P(n)-Pk)/Vm(n);
   end
   if kb(n) == 0
   iv=iv+1;
     DQ(iv) = Q(n)-Qk;
     DQV(iv) = (Q(n)-Qk)/Vm(n);
   end
end
Dd=-B1\DPV';
DV=-B2\DQV';
id=0;iv=0;
  for n=1:nbus
    if kb(n) ~= 1
    id = id+1;
    delta(n) = delta(n)+Dd(id); end
    if kb(n) == 0
    iv = iv+1;
    Vm(n)=Vm(n)+DV(iv); end
  end
    maxerror=max(max(abs(DP)),max(abs(DQ)));
   if iter == maxiter & maxerror > accuracy
   %fprintf('\nWARNING: Iterative solution did not converged after ')
  % fprintf('%g', iter), fprintf(' iterations.\n\n')
   %fprintf('Press Enter to terminate the iterations and print the results \n')
  % converge = 0; pause, 
   else, end
   
end
if converge ~= 1
   tech= ('                      ITERATIVE SOLUTION DID NOT CONVERGE'); else, 
   tech=('                   Power Flow Solution by Fast Decoupled Method');
end   
k=0;
V = Vm.*cos(delta)+j*Vm.*sin(delta);
deltad=180/pi*delta;
clear A; clear DC; clear DX
i=sqrt(-1);
for n = 1:nbus
     if kb(n) == 1
     S(n)=P(n)+j*Q(n);
     Pg(n) = P(n)*basemva + Pd(n);
     Qg(n) = Q(n)*basemva + Qd(n) - Qsh(n);
     k=k+1;
     Pgg(k)=Pg(n);
     elseif  kb(n) ==2
     S(n)=P(n)+j*Q(n);
     Qg(n) = Q(n)*basemva + Qd(n) - Qsh(n);
     k=k+1;
     Pgg(k)=Pg(n);
     end
yload(n) = (Pd(n)- j*Qd(n)+j*Qsh(n))/(basemva*Vm(n)^2);
end
busdata(:,3)=Vm'; busdata(:,4)=deltad';
Pgt = sum(Pg);  Qgt = sum(Qg); Pdt = sum(Pd); Qdt = sum(Qd); Qsht = sum(Qsh);
clear Pk Qk  DP DQ J11 J33 B1 B1inv B2 B2inv DPV  DQV Dd delta ib id ii iv jb jj
