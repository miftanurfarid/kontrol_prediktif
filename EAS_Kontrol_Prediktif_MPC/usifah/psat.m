function Psat=psat(T)
%fmls
%sistem methanol-air
%Input T dalam Kelvin.
%Data antoine T dalam Celsius dan P dalam kPa.
A=[16.5938 16.262];
B=[3644.3 3799.89];
C=[239.76 226.35];
T=T-273.15;
Psat=exp(A-B./(T+C));%Psat dalam kPa.
Psat=Psat/101.33;%konversi menjadi atm.