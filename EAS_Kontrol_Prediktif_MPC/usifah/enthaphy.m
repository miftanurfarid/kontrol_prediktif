function [Hliq,Hvap]=enthaphy(T,x,y)
%program untuk menghitung enthalpy liquid and vapor
%data-data didapat dari smith van ness
nc=2;
Al(1)=13.431; Bl(1)=-.051280; Cl(1)=1.31113E-4;
Al(2)=8.7120; Bl(2)=0.00125; Cl(2)=-1.8E-7;
Hlt(1)=8957.176413; Hlt(2)=10507.89879;%dalam cal/gmol 
Tc(1)=512.6; Tc(2)=647.1;
Tb(1)=337.85; Tb(2)=373.15;
R=1.987;%dalam cal/molK
Tref=298.15;
Hliq=0; Hvap=0;
det1=T-Tref;
det2=(T^2-Tref^2)/2;
det3=(T^3-Tref^3)/3;
for i=1:nc
    Hl(i)=Al(i)*det1+Bl(i)*det2+Cl(i)*det3;
    Hl1(i)=x(i)*Hl(i)*R;
    Hliq=Hliq+Hl1(i);
    Hltn(i)=Hlt(i)*((1-T/Tc(i))/(1-Tref/Tc(i)))^0.38;
    Hvap=Hvap+y(i)*(R*Hl(i)+Hltn(i));
end

