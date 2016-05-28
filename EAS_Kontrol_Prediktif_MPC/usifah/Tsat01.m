function Tsat=Tsat01(P)

%fmls
%sistem methanol-air
%data-data antoine untuk T dalam Celsius dan P dalam kPa
A=[16.5938 16.262];
B=[3644.3 3799.89];
C=[239.76 226.35];
P=P*101.33;
Tsat = B./(A-log(P))-C;
Tsat = Tsat + 273.15;   
%14-11-01 3:27
%Tsat(i) sudah dalam Kelvin