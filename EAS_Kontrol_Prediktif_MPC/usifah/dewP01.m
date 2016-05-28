function [x,P]=dewP01(nc,T,y)
%subprogram untuk menghitung bubble point dari suatu campuran
%algoritma perhitungan mengikuti buku "Introduction to Chemical
%Engineering Thermodynamics" karya Smith & van Ness
t1=T+273.15;
Psat=psat(T);%t dalam celsius.
gamma=[1 1];
tol=[0.00001 0.00001];
e=[1 1];
while abs(e)>tol
   Pk=0;
for i=1:nc
   Pk=Pk+y(i)/(gamma(i)*Psat(i));
end
P=1/Pk;
for i=1:nc;
   x(i)=y(i)*P/(gamma(i)*Psat(i));
end
gammaold=gamma;
   gamma=uniq01(x,t1)%t dalam Kelvin.
for i=1:nc;
   e(i)=(gammaold(i)-gamma(i))/gammaold(i);
end
end
