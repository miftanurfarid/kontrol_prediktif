function [y,T]=bubblet(nc,P,x)
%subprogram untuk menghitung bubble point dari suatu campuran
%algoritma perhitungan mengikuti buku "Introduction to Chemical
%Engineering Thermodynamics" karya Smith & van Ness
Tsat=Tsat01(P);%P dalam atm.
T=0;
for i=1:nc
   T=T+x(i)*Tsat(i);
end
Psat=psat(T);%T dalam Kelvin
gamma=uniq01(x,T);
tol=0.00001;
e=[1 1];
%break
while abs(e)>tol
   xx=0;
   for i=1: nc
      xx=xx+x(i)*gamma(i)*Psat(i)/Psat(1);
   end
   P1=P/xx;
   %Psat(1)=P/(x(1)*gamma(1)+x(2)*gamma(2)*Psat(2)/Psat(1));
   Tsat=Tsat01(P1);
   T=Tsat(1);
   Psat=psat(T);
   gamma=uniq01(x,T);
   Told=T;
   xx=0;
   for i=1:nc
      xx=xx+x(i)*gamma(i)*Psat(i)/Psat(1);
   end
   P1=P/xx;
   Tsat=Tsat01(P1);
   T=Tsat(1);
   e=T-Told;
end
for i=1:nc
   y(i)=x(i)*gamma(i)*Psat(i)/P;
end
%T keluar sudah dalam Kelvin.
%FMLS




