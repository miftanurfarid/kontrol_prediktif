function [x,T]=dewT01(nc,P,y)
%subprogram untuk menghitung dew point dari suatu campuran
%algoritma perhitungan mengikuti buku "Introduction to Chemical
%Engineering Thermodynamics" karya Smith & van Ness
Tsat=Tsat01(P);%P dalam atm.
T=0;
for i=1:nc
   T=T+y(i)*Tsat(i);
end
Psat=psat(T);%T dalam Kelvin
gamma=[1 1];
tol=0.00001;
e=[1 1];
while abs(e)>tol
   xx=0;
   for i=1:nc
      xx=xx+y(i)/gamma(i)*Psat(1)/Psat(i);
   end
   P1=P*xx;
   Tsat=Tsat01(P1);
   T=Tsat(1);
   Told=T;
   Psat=psat(T);   
   for i=1:nc;
      x(i)=y(i)*P/(gamma(i)*Psat(i));
   end
   gamma=uniq01(x,T);%t dalam Kelvin.
   xx=0;
   for i=1:nc
      xx=xx+y(i)/gamma(i)*Psat(1)/Psat(i);
   end
   P1=P*xx;
   Tsat=Tsat01(P1);
   T=Tsat(1);
   e=(Told-T)/Told;
end
%keluaran T sudah dalam Kelvin.
%17 des 2001 20:25



