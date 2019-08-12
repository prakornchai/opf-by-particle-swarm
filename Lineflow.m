%  This program is used in conjunction with lfgauss or lf Newton
%  for the computation of line flow and line losses.
%
%  Copyright (c) 1998 H. Saadat
SLT = 0;
%fprintf('\n')
%fprintf('                           Line Flow and Losses \n\n')
%fprintf('     --Line--  Power at bus & line flow    --Line loss--  Transformer\n')
%fprintf('     from  to    MW      Mvar     MVA       MW      Mvar      tap\n')
ddd=0;
for n = 1:nbus
busprt = 0;
   for L = 1:nbr;
       if busprt == 0
       %fprintf('   \n'), fprintf('%6g', n), fprintf('      %9.3f', P(n)*basemva)
       %fprintf('%9.3f', Q(n)*basemva), fprintf('%9.3f\n', abs(S(n)*basemva))

       busprt = 1;
       else, end
       if nl(L)==n      
       k = nr(L);
       In = (V(n) - a(L)*V(k))*y(L)/a(L)^2 + Bc(L)/a(L)^2*V(n);
       Ik = (V(k) - V(n)/a(L))*y(L) + Bc(L)*V(k);
       Snk = V(n)*conj(In)*basemva;
       Skn = V(k)*conj(Ik)*basemva;
       SL  = Snk + Skn;
       SLT = SLT + SL;
       elseif nr(L)==n  k = nl(L);
       In = (V(n) - V(k)/a(L))*y(L) + Bc(L)*V(n);
       Ik = (V(k) - a(L)*V(n))*y(L)/a(L)^2 + Bc(L)/a(L)^2*V(k);
       Snk = V(n)*conj(In)*basemva;
       Skn = V(k)*conj(Ik)*basemva;
       SL  = Snk + Skn;
       SLT = SLT + SL;
       else, end
         if nl(L)==n | nr(L)==n
         ddd=ddd+1;
         S_line(ddd)=abs(Snk);
         %fprintf('%12g', k),
         %fprintf('%9.3f', real(Snk)), fprintf('%9.3f', imag(Snk))
         %fprintf('%9.3f', abs(Snk)),
         %fprintf('%9.3f', real(SL)),
             %if nl(L) ==n & a(L) ~= 1
             %fprintf('%9.3f', imag(SL)), fprintf('%9.3f\n', a(L))
             %else, fprintf('%9.3f\n', imag(SL))
             %end
         else, end
  end
end
SLT = SLT/2;
%fprintf('   \n'), fprintf('    Total loss                         ')
%fprintf('%9.3f', real(SLT)), fprintf('%9.3f\n', imag(SLT))
Transmissionloss=real(SLT);
clear Ik In SL SLT Skn Snk
