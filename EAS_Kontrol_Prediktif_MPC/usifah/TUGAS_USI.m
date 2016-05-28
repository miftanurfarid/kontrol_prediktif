%----------------------------------------------------------------------------------
%-------------------         >>>  INITIALIZATIONS  <<<        ---------------------
%----------------------------------------------------------------------------------

clear all
close all
clc
load data_train_atas;
load data_train_bawah;
usi_initial
eval(['load ' nnfileA]);                % Load JST
eval(['load ' nnfileB]);                % Load JST
load metanol;

IAEA=0;
IAEB=0;
%__________________________ KOLOM DISTILASI _____________________________________
nc      = 2 ;                   % Jumlah komponen feed
nt      = length ( stage ) ;    % Jumlah tray dalam kolom distilasi
nf      = 5 ;                   % Nomor tray feed
Rset    = Reflux ;              % Laju aliran refluks awal
F       = 45000 ;               % Laju aliran feed 
beta    = 0.1 ;                 % 6 detik /  60  = 0.1 menit
xff     = [0.5 0.5] ;           % Komposisi komponen feed
xbset   = [0.01 0.99] ;         % Komposisi komponen bottom awal
xdset   = [0.99 0.01] ;         % Komposisi komponen distilat awal
mb      = 450000 ;              % Holdup pada bottom ( plate ke-0 )
md      = 450000 ;              % Holdup pada kondensor
Qc      = Qcon ;                % Panas kelur kondensor
Qro     = Qreb ;                % Panas masuk reboiler awal
delta   = 0.02 ;                % Interval waktu kolom distilasi
Dset    = 22500 ;               % Laju aliran distilat awal
Bset    = 22500 ;               % Laju aliran bottom awal
Lo      = 1000 * Liqr ;          
T       = Temper ;              % Temperatur tiap tray
Ptot    = 1.0 ;                 % Tekanan total
xset    = datax;                 
xbb     = xbset;                  
xb      = xbb;                     
xdd     = xdset;
xd      = xdd;
D       = Dset;
Bot     = Bset;
Qr      = Qro;
Reb     = Qro;
R       = Rset;
L       = Lo;
for i = 1 : length(xset)
   xset ( 2 , i ) = 1 - xset (1,i); 
end
for n = 1 : nt
    if n <= nf
        mo ( n ) = 6e3 ;
    else
        mo ( n ) = 2e3;
    end
m = mo;
    for i = 1:nc
        x(i,n)=xset(i,n);
        mp(i,n)=m(n)*x(i,n);
        mpx(i,n)=mp(i,n);
    end
end
% Reboiler
[ ybb , Tb ]        = bubblet ( nc, Ptot, xbb);
[ hlbold , hvbold ] = enthaphy( Tb, xbb, ybb );
% Tray feed
[ yff , Tf ]        = bubblet ( nc, Ptot, xff);
[ hlfold , hvfold ] = enthaphy( Tf, xff, yff );
% Tray ke-n
for n = 1:nt
   for j = 1:nc
      xx(j) = x(j,n);
   end
      [ yy, T(n) ] = bubblet ( nc, Ptot, xx);
      [ hlold(n), hvold(n) ] = enthaphy ( T(n), xx, yy );
   for j = 1:nc
      x ( j, n ) = xx(j);
      y ( j, n ) = yy(j);
   end
end
% Kondensor
[ ydd , Td ]        = bubblet ( nc, Ptot, xdd );
[ hldold , hvdold]  = enthaphy ( Td, xdd, ydd );

%________________________________ DISTILAT _____________________________
% >>>>>>>>>>>>>>>>>>>>>>>>   MENENTUKAN STRUKTUR REGRESSION   <<<<<<<<<<<<<<<<<<<<<<   
naA      = NNA(1);              % Jumlah y lampau yang dipakai
nbA      = NNA(2);              % Jumlah u lampau yang dipakai
nkA      = 1;                   % Time delay sistem
dA       = 1;                   % Time delay in addition to the usual 1
N1A      = dA;                  % N1<>d tidak dipakai 
inputsA  = naA+sum(nbA);        % Jumlah input total JST
outputsA = 1;                   % Jumlah output sistem (1 jika SISO)
phiA     = zeros(inputsA,1);    % Inisialisasi vektor regression


% >>>>>>>>>>>>>>>>>    MENENTUKAN STRUKTUR MODEL JST     <<<<<<<<<<<<<<<<<<<
hiddenA   = length(NetDefA(1,:));        % Jumlah hidden node
L_hiddenA = find(NetDefA(1,:)=='L')';    % Lokasi hidden node dg fungsi aktifasi linear
H_hiddenA = find(NetDefA(1,:)=='H')';    % Lokasi hidden node dg fungsi aktifasi tangent hiperbolik
y1A       = zeros(hiddenA,N2A-N1A+1);    % Output hidden layer
yhatA     = zeros(outputsA,1);           % Output JST


%>>>>>>>>>>>>>>>>>>>>>>>        INISIALISASI VARIABEL        <<<<<<<<<<<<<<<<<<<<<<
% Inisialisasi sinyal lampau
maxlengthA = 5;                          
ref_oldA   = repmat(y_0A,maxlengthA,1);    
y_oldA     = repmat(y_0A,N2A,1);
yhat_oldA  = repmat(y_0A,N2A,1);
u_oldA     = repmat(u_0A,maxlengthA,1);


dfA         = ones(hiddenA,1);
d2fA        = ones(hiddenA,1);
dUtilde_dUA = eye(NuA);
dUtilde_dUA(1:NuA-1,2:NuA)=dUtilde_dUA(1:NuA-1,2:NuA)-eye(NuA-1);
d2UA        = rhoA'*dUtilde_dUA*dUtilde_dUA';
dY_dUA      = zeros(NuA,N2A-N1A+1);             % Inisialisasi matriks turunan parsial
hj_matA     = zeros(hiddenA*NuA,N2A-dA+1);
d2Y_dU2A    = zeros(NuA*NuA,N2A-dA+1);

uA          = u_0A;                             % Sinyal kontrol sampai t<=0
upA         = uA(ones(NuA,1));                  % Inisialisasisinyal kontrol masa depan
upminA      = upA;
index1A     = 1:NuA:NuA*(NuA-1)+1;   
index2A     = 1:hiddenA:hiddenA*(hiddenA-1)+1; 
index3A     = 1:(NuA+1):NuA^2;      

tA          = -TsA;
delta2A     = deltaA*deltaA;
fighandle   = progress;



%>>>>>>>>>>>>>>>>    MENENTUKAN SET POINT        <<<<<<<<<<<<<<<<<<
  eval(['refA = ' reftyA ';']);
  refA=refA(:);
  iA=length(refA);
  if iA>samplesA+N2A,
    refA=refA(1:samplesA+N2A);
  else
    refA=[refA;refA(iA)*ones(samplesA+N2A-iA,1)];
  end


% Inisialisasi vektor yang digunakan untuk penyimpanan data lampau
ref_dataA    = refA(1:samplesA);
u_dataA      = zeros(samplesA,1);
y_dataA      = zeros(samplesA,1);
yhat_dataA   = zeros(samplesA,1);
t_dataA      = zeros(samplesA,1);
e_dataA      = zeros(samplesA,1);

% Vektor data yang digunakan dalam algoritma kontrol
% u_vec : Elemen pertama adalah u(t-d+N1-nb+1), elemen terakhir adalah u(t-d+N2)
% Index waktu : tiu = d-N1+nb 
u_vecA = repmat(uA,N2A-N1A+nbA,1);
tiuA   = dA-N1A+nbA;
upiA   = [1:NuA-1 NuA(ones(1,N2A-dA-NuA+2))];    
uviA   = [tiuA:N2A-N1A+nbA];

% y_vec: Elemen pertama adalah y(t-na), elemen terakhir adalah y(t-1)
% Index waktu : tiy=na+1
y_vecA = repmat(y_0A,naA,1);
tiyA   = naA+1;
% yhat_vec: Elemen pertama adalah yhat(t), elemen terakhir adalah yhat(t+N2)
% Index waktu : tiyh=1
yhat_vecA = repmat(y_0A,N2A+1,1);
tiyhA     = 1;

%________________________________ BOTTOM ______________________________
% >>>>>>>>>>>>>>>>>>>>>>>>   MENENTUKAN STRUKTUR REGRESSION   <<<<<<<<<<<<<<<<<<<<<<   
naB      = NNB(1);              % Jumlah y lampau yang dipakai
nbB      = NNB(2);              % Jumlah u lampau yang dipakai
nkB      = 1;                   % Time delay sistem
dB       = 1;                   % Time delay in addition to the usual 1
N1B      = dB;                  % N1<>d tidak dipakai 
inputsB  = naB+sum(nbB);        % Jumlah input total JST
outputsB = 1;                   % Jumlah output sistem (1 jika SISO)
phiB     = zeros(inputsB,1);    % Inisialisasi vektor regression 


% >>>>>>>>>>>>>>>>>    MENENTUKAN STRUKTUR MODEL JST     <<<<<<<<<<<<<<<<<<<
hiddenB   = length(NetDefB(1,:));        % Jumlah hidden node
L_hiddenB = find(NetDefB(1,:)=='L')';    % Lokasi hidden node dg fungsi aktifasi linier
H_hiddenB = find(NetDefB(1,:)=='H')';    % Lokasi hidden node dg fungsi aktifasi tangen hiperbolik
y1B       = zeros(hiddenB,N2B-N1B+1);    % Output hidden layer
yhatB     = zeros(outputsB,1);           % Output JST


%>>>>>>>>>>>>>>>>>>>>>>>        INITIALIZE VARIABLES        <<<<<<<<<<<<<<<<<<<<<<
% Initialization of past signals
maxlengthB = 5;                          
ref_oldB   = repmat(y_0B,maxlengthB,1);    
y_oldB     = repmat(y_0B,N2B,1);
yhat_oldB  = repmat(y_0B,N2B,1);
u_oldB     = repmat(u_0B,maxlengthB,1);


% Miscellaneous initializations
dfB         = ones(hiddenB,1);
d2fB        = ones(hiddenB,1);
dUtilde_dUB = eye(NuB);
dUtilde_dUB(1:NuB-1,2:NuB) = dUtilde_dUB(1:NuB-1,2:NuB)-eye(NuB-1);
d2UB        = rhoB'*dUtilde_dUB*dUtilde_dUB';
dY_dUB      = zeros(NuB,N2B-N1B+1);             % Inisialisasi matriks turunan parsial
hj_matB     = zeros(hiddenB*NuB,N2B-dB+1);
d2Y_dU2B    = zeros(NuB*NuB,N2B-dB+1);

uB          = u_0B;                 % Sinyal kontrol sampai t<=0
upB         = uB(ones(NuB,1));      % Inisialisasi sinyal kontrol masa depan
upminB      = upB;
index1B     = 1:NuB:NuB*(NuB-1)+1;  % A useful vector
index2B     = 1:hiddenB:hiddenB*(hiddenB-1)+1; %Another useful vector
index3B     = 1:(NuB+1):NuB^2;      % A third useful vector

tB          = -TsB;
delta2B     = deltaB*deltaB;



%>>>>>>>>>>>>>>>>     MENENTUKAN SET POINT     <<<<<<<<<<<<<<<<<<
  eval(['refB = ' reftyB ';']);
  refB=refB(:);
  iB=length(refB);
  if iB>samplesB+N2B,
    refB=refB(1:samplesB+N2B);
  else
    refB=[refB;refB(iB)*ones(samplesB+N2B-iB,1)];
  end

% Inisialisasi vektor untuk menyimpan data masa lalu
ref_dataB    = refB(1:samplesB);
u_dataB      = zeros(samplesB,1);
y_dataB      = zeros(samplesB,1);
yhat_dataB   = zeros(samplesB,1);
t_dataB      = zeros(samplesB,1);
e_dataB      = zeros(samplesB,1);

% Vektor data yang dipakai di algoritma kontrol
% u_vec : Elemen pertama adalah u(t-d+N1-nb+1), elemen terakhir adalah is u(t-d+N2)
% Index waktu  : tiu=d-N1+nb 
u_vecB = repmat(uB,N2B-N1B+nbB,1);
tiuB   = dB-N1B+nbB;
upiB   = [1:NuB-1 NuB(ones(1,N2B-dB-NuB+2))];    
uviB   = [tiuB:N2B-N1B+nbB];

% y_vec : Elemen pertama adalah is y(t-na), elemen terakhir adalah y(t-1)
% Index waktu : tiy=na+1
y_vecB = repmat(y_0B,naB,1);
tiyB   = naB+1;

% yhat_vec : Elemen pertama adalah yhat(t), elemen terakhir adalah is yhat(t+N2)
% Index waktu : tiyh=1
yhat_vecB = repmat(y_0B,N2B+1,1);
tiyhB     = 1;

%------------------------------------------------------------------------------
%-------------------         >>>   MAIN LOOP   <<<           ------------------
%------------------------------------------------------------------------------
c = fix(clock);
fprintf('Pengendalian dimulai pukul %2i.%2i.%2i\n\n',c(4),c(5),c(6));
for iB = 1:samplesB,
    tB = tB +delta;
    
  %>>>>>>>>>>>>>>>>>>>>>>>>    PREDIKSI OUTPUT DARI PLANT  <<<<<<<<<<<<<<<<<<<<<<<
  phiA             = [y_vecA(naA:-1:1);u_oldA(dA:dA+nbA-1)];
  h1A              = W1A(:,1:inputsA)*phiA + W1A(:,inputsA+1);  
  y1A(H_hiddenA,1) = pmntanh(h1A(H_hiddenA)); 
  y1A(L_hiddenA,1) = h1A(L_hiddenA);
  yhatA            = W2A(:,1:hiddenA)*y1A(:,1) + W2A(:,hiddenA+1);
  yhat_vecA(tiyhA) = yhatA;

%  uji disturbance
%  if iA>=100
%      yhat_vecA(tiyhA)=yhatA-(0.0000001*yhatA);
%  end
 
   %>>>>>>>>>>>>>>>>>>>>>>>>    PREDIKSI OUTPUT PLANT  <<<<<<<<<<<<<<<<<<<<<<<
   phiB             = [y_vecB(naB:-1:1);u_oldB(dB:dB+nbB-1)];
   h1B              = W1B(:,1:inputsB)*phiB + W1B(:,inputsB+1);  
   y1B(H_hiddenB,1) = pmntanh(h1B(H_hiddenB)); 
   y1B(L_hiddenB,1) = h1B(L_hiddenB);
   yhatB            = W2B(:,1:hiddenB)*y1B(:,1) + W2B(:,hiddenB+1);
   yhat_vecB(tiyhB) = yhatB;

%  uji disturbance
%  if iB>=200
%      yhat_vecB(tiyhB)=yhatB-(0.001*yhatB);
%  end

%________________________________ KOLOM DISTILASI __________________________
   %>>>>>>>>>>>>>>>>>>>>>>>>    MEMBACA OUTPUT DARI PLANT     <<<<<<<<<<<<<<<<<<<<<<<
   R  = Reflux;
   Qr = uB;
   L = Lo+(m-mo)/beta;      
   [ ybb, Tb ]  = bubblet(nc,Ptot,xbb);
   [ hlb, hvb ] = enthaphy(Tb,xbb,ybb);
   [ yff, Tf ]  = bubblet(nc,Ptot,xff);
   [ hlf, hvf ] = enthaphy(Tf,xff,yff);
   for n = 1:nt
       for j = 1:nc
           xx(j) = x(j,n);
       end
       [ yy, T(n) ]    = bubblet(nc,Ptot,xx);
       [ hl(n), hv(n)] = enthaphy(T(n),xx,yy);
       for j = 1:nc
           x(j,n)=xx(j);
           y(j,n)=yy(j);
       end
   end
   [ ydd, Td ] = bubblet(nc,Ptot,xdd);
   [ hld, hvd] = enthaphy(Td,xdd,ydd);
   Vb  = (L(1)*(hl(1)-hlb)+Qr-mb*(hlb-hlbold)/delta)/(hvb-hlb);
   Bot = L(1)-Vb;
   if Bot < 0.0
       disp(['L1 = ',num2str(L(1))]);
       disp(['VB = ',num2str(Vb)]);
       disp('NEGATIVE BOTTOM RATE ENCOUNTERED');
       pause
   end
   yvap = VRATE1(nt,nf,L,hl,hlold,hv,hlb,Qr,Vb,F,R,hvb,hlf,hld,delta,m);
   n    = nt;
   hvn  = hv(n);
   yvapn= yvap(n);
   Qc   = yvapn*(hvn-hld)-md*(hld-hldold)/delta;
   D    = yvapn-R;
   dm   = DERIV1(nt,nf,L,yvap,Vb,F,R);
   [dxb,dxm,dxd] = DERIV2(nt,nf,L,yvap,x,y,xbb,ybb,xdd,xff,Vb,Bot,R,D,F,mb,md);
   for n=1:nt
       m(n)=m(n)+dm(n)*delta;
   end
   for i=1:nc
       xbb(i)=xbb(i)+dxb(i)*delta;
       if xbb(i)<=0.0
           xbb(i)=0.001;
       end
       xdd(i)=xdd(i)+dxd(i)*delta;
       if xdd(i)<=0.0
           xdd(i)=0.001;
       end
   end
   for n=1:nt
       for i=1:nc
           mpx(i,n)=mpx(i,n)+dxm(i,n)*delta;
           x(i,n)=mpx(i,n)/m(n);
           if x(i,n)<=0.0
               x(i,n)=0.001;
           end
       end
   end
   sum1=0.0;
   sum2=0.0;
   for i = 1:nc
       sum1=sum1+xbb(i);
       sum2=sum2+xdd(i);
   end
   for i=1:nc
       xbb(i)=xbb(i)/sum1;
       erxb(i)=xbset(i)-xbb(i);
       xdd(i)=xdd(i)/sum2;
       erxd(i)=xdset(i)-xdd(i);
   end
   for n=1:nt
       sum3=0.0;
       for i=1:nc
           sum3=sum3+x(i,n);
       end
       for i=1:nc
           x(i,n)=x(i,n)/sum3;
       end
   end
   hlbold=hlb;
   hldold=hld;
   for n=1:nt
       hlold(n)=hl(n);
   end
   for i=1:nc
       xxb(i,1)=xbb(i);
       xxd(i,1)=xdd(i);
   end
   yB=yhatB;     
   yA=yhatA;
    
%________________________________ DISTILAT _________________________
  %>>>>>>>>>>>>>>>>>>  MENCARI SINYAL KONTROL BARU DENGAN OPTIMASI  <<<<<<<<<<<<<<<<<
  upmin0A   = upminA([2:NuA NuA]);
  einitvalA = eval(initvalA);     
  for trA=1:length(einitvalA),
      upA=upmin0A;                 % Harga awal untuk pencarian u baru  
      upA(NuA)=einitvalA(trA);
      u_vecA(uviA) = upA(upiA);  
      dwA = 1;                     
      lambdaA = 0.1;               % Inisialisasi parameter Levenberg-Marquardt 
      mA = 1; 
  
  
      %>>>>>>>>>>>>>>> COMPUTE PREDICTIONS FROM TIME t+N1 TO t+N2 <<<<<<<<<<<<<<<<
      for kA=N1A:N2A,
          %----- Determine prediction yhat(t+k) -----
          phiA              = [yhat_vecA(tiyhA+kA-1:-1:tiyhA+kA-min(kA,naA)) ; ...
                  y_vecA(tiyA-1:-1:tiyA-max(naA-kA,0)) ; u_vecA(tiuA-dA+kA:-1:tiuA-dA+1-nbA+kA)];
          h1A               = W1A(:,1:inputsA)*phiA + W1A(:,inputsA+1);  
          y1A(H_hiddenA,kA-N1A+1) = pmntanh(h1A(H_hiddenA)); 
          y1A(L_hiddenA,kA-N1A+1) = h1A(L_hiddenA);
          yhat_vecA(tiyhA+kA) = W2A(:,1:hiddenA)*y1A(:,kA-N1A+1) + W2A(:,hiddenA+1);
      end
      %>>>>>>>>>>>>>>>>>>>>>>    EVALUATE CRITERION    <<<<<<<<<<<<<<<<<<<<<<
      duvecA = u_vecA(tiuA:tiuA+NuA-1)-u_vecA(tiuA-1:tiuA+NuA-2);
      evecA  = refA(iB+N1A:iB+N2A) - yhat_vecA(tiyhA+N1A:tiyhA+N2A);
      JA = evecA'*evecA + rhoA*(duvecA'*duvecA);
% if iB>=200
%     JA=JA+2*JA;
% end
  
      while mA<=maxiterA,
          if dwA == 1,
              %>>>>>>>>>>>>>>>>>>>    DETERMINE dy/du <<<<<<<<<<<<<<<<<<<<<
              for kA=N1A:N2A
                  dfA(H_hiddenA)  = (1-y1A(H_hiddenA,kA-N1A+1).*y1A(H_hiddenA,kA-N1A+1));
                  d2fA(H_hiddenA) = -2*y1A(H_hiddenA,kA-N1A+1).*dfA(H_hiddenA);
                  for lA=0:min(kA-dA,NuA-2)
                      imax1A = min(kA-dA-lA,naA);
                      if lA>=kA-dA-nbA+1,
                          if imax1A>=1,
                              hj_vecA = W1A(:,1:imax1A)*dY_dUA(lA+1,kA-N1A:-1:kA-imax1A-N1A+1)' + W1A(:,naA+kA-dA-lA+1);
                          else
                              hj_vecA=W1A(:,naA+kA-dA-lA+1);
                          end
                      else
                          hj_vecA = W1A(:,1:imax1A)*dY_dUA(lA+1,kA-N1A:-1:kA-imax1A-N1A+1)';
                      end
                      hj_matA(index2A(lA+1):index2A(lA+1)+hiddenA-1,kA-dA+1) = hj_vecA;
                      dY_dUA(lA+1,kA-N1A+1)  = W2A(1:hiddenA)*(dfA.*hj_vecA);
                  end
                  if kA>=NuA
                      lA=NuA-1;
                      imax1A = min(kA-dA-lA,naA);
                      imax2A = min(kA-dA-NuA+2,nbA);
                      if imax2A>1,
                          if kA==NuA,
                              hj_vecA = sum(W1A(:,naA+1:naA+imax2A)')';
                          else
                              hj_vecA = W1A(:,1:imax1A)*dY_dUA(lA+1,kA-N1A:-1:kA-imax1A-N1A+1)'...
                                  + sum(W1A(:,naA+1:naA+imax2A)')';
                          end
                      else
                          if kA==NuA,
                              hj_vecA = W1A(:,naA+1:naA+imax2A);
                          else
                              hj_vecA = W1A(:,1:imax1A)*dY_dUA(lA+1,kA-N1A:-1:kA-imax1A-N1A+1)'...
                                  + W1A(:,naA+1:naA+imax2A);
                          end
                      end
                      hj_matA(index2A(lA+1):index2A(lA+1)+hiddenA-1,kA-dA+1) = hj_vecA;
                      dY_dUA(lA+1,kA-N1A+1)  = W2A(1:hiddenA)*(dfA.*hj_vecA);
                  end 
              end
              %>>>>>>>>>>>>>>>>>>>    DETERMINE d^2 y / du^2 <<<<<<<<<<<<<<<<<<<<<
              for kA=dA:N2A,
                  for lA=0:NuA-1,
                      imax1A = min(kA-dA-lA,naA);
                      for pA=0:lA,
                          if imax1A>=1,
                               tmpvecA = W1A(:,1:imax1A)*d2Y_dU2A(index1A(lA+1)+pA,kA-dA:-1:kA-imax1A-dA+1)';
                           else
                               tmpvecA = zeros(hiddenA,1);
                           end         
                           ddA = d2fA.*hj_matA(index2A(lA+1):index2A(lA+1)+hiddenA-1,kA-dA+1)...
                               .*hj_matA(index2A(pA+1):index2A(pA+1)+hiddenA-1,kA-dA+1);
                           ddA = ddA + dfA.*tmpvecA;
                           d2Y_dU2A(index1A(lA+1)+pA,kA-dA+1) = W2A(1:hiddenA)*ddA;
                       end
                   end
               end
               %>>>>>>>>>>>>>>>>>>>    DETERMINE Hessian matrix <<<<<<<<<<<<<<<<<<<<<
               HA = dY_dUA*dY_dUA' + d2UA;
               for lA=1:NuA,
                   for pA=1:lA,
                       HA(lA,pA) = HA(lA,pA) - d2Y_dU2A(index1A(lA)+pA-1,:)*evecA;
                       HA(pA,lA) = HA(lA,pA);              
                   end
               end
               % -- Gradient --
               dJduA   = -dY_dUA*evecA + rhoA*(dUtilde_dUA*duvecA);
               dwA = 0;   
           end
           %>>>>>>>>>>>>>>>>>>>>>   DETERMINE SEARCH DIRECTION   <<<<<<<<<<<<<<<<<<<<<
           % -- Add lambda to the diagonal --
           HA(index3A) = HA(index3A) + lambdaA;
           % -- Check if H is positive definite --
           pA=1;
           while pA~=0,
               [chA,pA] = chol(HA);
               if pA~=0,         % Hessian + diagonal not positive definite
                   lambdaA = lambdaA*4;
                   HA(index3A) = HA(index3A) + lambdaA;
               end
           end
           % -- Determine search direction --
           fA = -(chA\(chA'\dJduA));
           % -- Compute 'apriori' iterate --
           up_newA= upA + fA;
           u_vecA(uviA) = up_newA(upiA);          % Insert updated controls
           
           %>>>>>>>>>>>>> COMPUTE PREDICTIONS FROM TIME t+N1 TO t+N2  <<<<<<<<<<<<<<<< 
           for kA=N1A:N2A,
               %----- Determine prediction yhat(t+k) -----
               phiA              = [yhat_vecA(tiyhA+kA-1:-1:tiyhA+kA-min(kA,naA)) ; ...
                       y_vecA(tiyA-1:-1:tiyA-max(naA-kA,0)) ; u_vecA(tiuA-dA+kA:-1:tiuA-dA+1-nbA+kA)];
               h1A               = W1A(:,1:inputsA)*phiA + W1A(:,inputsA+1);  
               y1A(H_hiddenA,kA-N1A+1) = pmntanh(h1A(H_hiddenA)); 
               y1A(L_hiddenA,kA-N1A+1) = h1A(L_hiddenA);
               yhat_vecA(tiyhA+kA) = W2A(:,1:hiddenA)*y1A(:,kA-N1A+1) + W2A(:,hiddenA+1);
           end
            
           %>>>>>>>>>>>>>>>>>>>>>>>>    EVALUATE CRITERION   <<<<<<<<<<<<<<<<<<<<<<<<<<  
           duvecA = u_vecA(tiuA:tiuA+NuA-1)-u_vecA(tiuA-1:tiuA+NuA-2);
           evec_newA  = refA(iB+N1A:iB+N2A) - yhat_vecA(tiyhA+N1A:tiyhA+N2A);
           J_newA = evec_newA'*evec_newA + rhoA*(duvecA'*duvecA);


           %>>>>>>>>>>>>>>>>>>>>>>>       UPDATE lambda      <<<<<<<<<<<<<<<<<<<<<<<<<<
           LA = -fA'*dJduA +0.5*lambdaA*fA'*fA;
     
           % Reduce lambda if SSE has decreased 'sufficiently'
           if (JA - J_newA) > (0.75*LA),
               lambdaA = lambdaA/2;
               % Increase lambda if SSE has increased 'sufficiently'
           elseif (JA-J_newA) <= (0.25*LA),
               lambdaA = 2*lambdaA;
           end
           
           %>>>>>>>>>>>>>>>>>>>>>   UPDATES FOR NEXT ITERATION   <<<<<<<<<<<<<<<<<<<<<< 
           % Update only if criterion has decreased
           if J_newA < JA,
               evecA = evec_newA;
               duA = up_newA-upA;                
               upA = up_newA;
               JA = J_newA;
               dwA = 1;
               mA = mA + 1;
               % Check if stop condition is satisfied 
               if sqrt(duA'*duA)<delta2A, break, end 
           end
           if lambdaA>1e3, break, end       
       end
       
       %>>>>>>>>>>>>>>>>>>>>>>>     SELECT BEST MINIMUM     <<<<<<<<<<<<<<<<<<<<<<<<<
       if trA==1,
           Jmin_oldA = JA;
           upminA = upA;
       else
           if JA<Jmin_oldA,
               upminA = upA;
           end
       end
   end
   
   
   %>>>>>>>>>>>>>>>>>>>>>>     CALCULATE CONTROL SIGANL     <<<<<<<<<<<<<<<<<<<<<<
   eA = refA(iB) - yA;
   
   % Predictive Controller
   uA=upminA(1);
  
   % Make sure control input is within limits
   if uA>ulim_maxA,
       uA=ulim_maxA;
   elseif uA<ulim_minA
       uA=ulim_minA;
   end
   upminA(1) = uA;

   %>>>>>>>>>>>>>>>>>>>       STORE DATA IN DATA VECTORS      <<<<<<<<<<<<<<<<<<<
   u_dataA(iB)       = uA;
   y_dataA(iB)       = yA;
   yhat_dataA(iB)    = yhat_vecA(tiyhA);
   t_dataA(iB)       = tA;
   J_dataA(iB)       = JA;
   e_dataA(iB)       = eA;
   %>>>>>>>>>>>>>>>>>>>>>>>>>>       TIME UPDATES        <<<<<<<<<<<<<<<<<<<<<<<<<
   y_oldA    = [yA;y_oldA(1:end-1)];
   u_oldA    = [uA;u_oldA(1:end-1)];
   u_vecA(uviA) = upminA(upiA);
   u_vecA    = [u_vecA(2:length(u_vecA)) ; upminA(length(upA))];
   y_vecA    = [y_vecA(2:length(y_vecA)) ; yA];
   yhat_vecA(1:length(yhat_vecA)-1) = yhat_vecA(2:length(yhat_vecA));

eA   = (((max(xdA)-min(xdA)))*(eA+1)/2)+min(xdA);%descaling
IAEA = IAEA+abs(eA)*delta;




%__________________________________ BOTTOM _________________________
  %>>>>>>>>>>>>>>>>>>  MENCARI SINYAL KONTROL MELALUI OPTIMASI  <<<<<<<<<<<<<<<<<
  upmin0B   = upminB([2:NuB NuB]);
  einitvalB = eval(initvalB);     
  for trB=1:length(einitvalB),
      upB=upmin0B;                  % Nilai awal untuk pencarian u baru  
      upB(NuB)=einitvalB(trB);
      u_vecB(uviB) = upB(upiB);  
      dwB = 1;                    
      lambdaB = 0.1;                % Inisislisasi parameter Levenberg-Marquardt 
      mB = 1; 
  
  
  %>>>>>>>>>>>>>>> MENENTUKAN PREDIKSI DARI t+N1 TO t+N2 <<<<<<<<<<<<<<<<
      for kB=N1B:N2B,
          %----- Menentukan prediksi yhat(t+k) -----
          phiB                    = [yhat_vecB(tiyhB+kB-1:-1:tiyhB+kB-min(kB,naB)) ; ...
                  y_vecB(tiyB-1:-1:tiyB-max(naB-kB,0)) ; u_vecB(tiuB-dB+kB:-1:tiuB-dB+1-nbB+kB)];
          h1B                     = W1B(:,1:inputsB)*phiB + W1B(:,inputsB+1);  
          y1B(H_hiddenB,kB-N1B+1) = pmntanh(h1B(H_hiddenB)); 
          y1B(L_hiddenB,kB-N1B+1) = h1B(L_hiddenB);
          yhat_vecB(tiyhB+kB)     = W2B(:,1:hiddenB)*y1B(:,kB-N1B+1) + W2B(:,hiddenB+1);
      end
  
  
  %>>>>>>>>>>>>>>>>>>>>>>    EVALUASI KRITERIA    <<<<<<<<<<<<<<<<<<<<<<
      duvecB = u_vecB(tiuB:tiuB+NuB-1)-u_vecB(tiuB-1:tiuB+NuB-2);
      evecB  = refB(iB+N1B:iB+N2B) - yhat_vecB(tiyhB+N1B:tiyhB+N2B);
      JB     = evecB'*evecB + rhoB*(duvecB'*duvecB);
  
      while mB<=maxiterB,
          if dwB == 1,
              %>>>>>>>>>>>>>>>>>>>    DETERMINE dy/du <<<<<<<<<<<<<<<<<<<<<
              for kB=N1B:N2B
                  dfB(H_hiddenB)  = (1-y1B(H_hiddenB,kB-N1B+1).*y1B(H_hiddenB,kB-N1B+1));
                  d2fB(H_hiddenB) = -2*y1B(H_hiddenB,kB-N1B+1).*dfB(H_hiddenB);
                  for lB=0:min(kB-dB,NuB-2)
                      imax1B = min(kB-dB-lB,naB);
                      if lB>=kB-dB-nbB+1,
                          if imax1B>=1,
                              hj_vecB = W1B(:,1:imax1B)*dY_dUB(lB+1,kB-N1B:-1:kB-imax1B-N1B+1)' + W1B(:,naB+kB-dB-lB+1);
                          else
                              hj_vecB=W1B(:,naB+kB-dB-lB+1);
                          end
                      else
                          hj_vecB = W1B(:,1:imax1B)*dY_dUB(lB+1,kB-N1B:-1:kB-imax1B-N1B+1)';
                      end
                      hj_matB(index2B(lB+1):index2B(lB+1)+hiddenB-1,kB-dB+1) = hj_vecB;
                      dY_dUB(lB+1,kB-N1B+1)  = W2B(1:hiddenB)*(dfB.*hj_vecB);
                  end
                  if kB>=NuB
                      lB=NuB-1;
                      imax1B = min(kB-dB-lB,naB);
                      imax2B = min(kB-dB-NuB+2,nbB);
                      if imax2B>1,
                          if kB==NuB,
                              hj_vecB = sum(W1B(:,naB+1:naB+imax2B)')';
                          else
                              hj_vecB = W1B(:,1:imax1B)*dY_dUB(lB+1,kB-N1B:-1:kB-imax1B-N1B+1)'...
                                  + sum(W1B(:,naB+1:naB+imax2B)')';
                          end
                      else
                          if kB==NuB,
                              hj_vecB = W1B(:,naB+1:naB+imax2B);
                          else
                              hj_vecB = W1B(:,1:imax1B)*dY_dUB(lB+1,kB-N1B:-1:kB-imax1B-N1B+1)'...
                                  + W1B(:,naB+1:naB+imax2B);
                          end
                      end
                      hj_matB(index2B(lB+1):index2B(lB+1)+hiddenB-1,kB-dB+1) = hj_vecB;
                      dY_dUB(lB+1,kB-N1B+1)  = W2B(1:hiddenB)*(dfB.*hj_vecB);
                  end 
              end
   
   
              %>>>>>>>>>>>>>>>>>>>    MENENTUKAN d^2 y / du^2 <<<<<<<<<<<<<<<<<<<<<
              for kB=dB:N2B,
                  for lB=0:NuB-1,
                      imax1B = min(kB-dB-lB,naB);
                      for pB=0:lB,
                          if imax1B>=1,
                              tmpvecB = W1B(:,1:imax1B)*d2Y_dU2B(index1B(lB+1)+pB,kB-dB:-1:kB-imax1B-dB+1)';
                          else
                              tmpvecB = zeros(hiddenB,1);
                          end         
                          ddB = d2fB.*hj_matB(index2B(lB+1):index2B(lB+1)+hiddenB-1,kB-dB+1)...
                              .*hj_matB(index2B(pB+1):index2B(pB+1)+hiddenB-1,kB-dB+1);
                          ddB = ddB + dfB.*tmpvecB;
                          d2Y_dU2B(index1B(lB+1)+pB,kB-dB+1) = W2B(1:hiddenB)*ddB;
                      end
                  end
              end
   
   
              %>>>>>>>>>>>>>>>>>>>    DETERMINE Hessian matrix <<<<<<<<<<<<<<<<<<<<<
              HB = dY_dUB*dY_dUB' + d2UB;
              for lB=1:NuB,
                  for pB=1:lB,
                      HB(lB,pB) = HB(lB,pB) - d2Y_dU2B(index1B(lB)+pB-1,:)*evecB;
                      HB(pB,lB) = HB(lB,pB);             
                  end
              end
              % -- Gradient --
              dJduB   = -dY_dUB*evecB + rhoB*(dUtilde_dUB*duvecB);
              dwB = 0;   
          end
          %>>>>>>>>>>>>>>>>>>>>>   DETERMINE SEARCH DIRECTION   <<<<<<<<<<<<<<<<<<<<<
          % -- Menambahkan lambda ke diagonal --
          HB(index3B) = HB(index3B) + lambdaB;
    
          % -- Check if H is positive definite --
          pB=1;
          while pB~=0,
              [chB,pB] = chol(HB);
              if pB~=0,         % Hessian + diagonal tidak berharga positif
                  lambdaB = lambdaB*4;
                  HB(index3B) = HB(index3B) + lambdaB;
              end
          end
          % -- Determine search direction --
          fB = -(chB\(chB'\dJduB));

          % -- Compute 'apriori' iterate --
          up_newB = upB + fB;
          u_vecB(uviB) = up_newB(upiB);          % Menyisipkan sinyal kontrol ter-update
          %>>>>>>>>>>>>> COMPUTE PREDICTIONS FROM TIME t+N1 TO t+N2  <<<<<<<<<<<<<<<< 
          for kB=N1B:N2B,
              %----- Menentukan prediksi yhat(t+k) -----
              phiB              = [yhat_vecB(tiyhB+kB-1:-1:tiyhB+kB-min(kB,naB)) ; ...
                      y_vecB(tiyB-1:-1:tiyB-max(naB-kB,0)) ; u_vecB(tiuB-dB+kB:-1:tiuB-dB+1-nbB+kB)];
              h1B               = W1B(:,1:inputsB)*phiB + W1B(:,inputsB+1);  
              y1B(H_hiddenB,kB-N1B+1) = pmntanh(h1B(H_hiddenB)); 
              y1B(L_hiddenB,kB-N1B+1) = h1B(L_hiddenB);
              yhat_vecB(tiyhB+kB) = W2B(:,1:hiddenB)*y1B(:,kB-N1B+1) + W2B(:,hiddenB+1);
          end
   
          %>>>>>>>>>>>>>>>>>>>>>>>>    EVALUATE CRITERION   <<<<<<<<<<<<<<<<<<<<<<<<<<  
          duvecB = u_vecB(tiuB:tiuB+NuB-1)-u_vecB(tiuB-1:tiuB+NuB-2);
          evec_newB  = refB(iB+N1B:iB+N2B) - yhat_vecB(tiyhB+N1B:tiyhB+N2B);
          J_newB = evec_newB'*evec_newB + rhoB*(duvecB'*duvecB);

 
          %>>>>>>>>>>>>>>>>>>>>>>>       UPDATE lambda      <<<<<<<<<<<<<<<<<<<<<<<<<<
          LB = -fB'*dJduB +0.5*lambdaB*fB'*fB;
     
          % Menurunkan lambda jika SSE menurun 'sufficiently'
          if (JB - J_newB) > (0.75*LB),
              lambdaB = lambdaB/2;
              % Meningkatkan lambda if SSE has meningkat 'sufficiently'
          elseif (JB-J_newB) <= (0.25*LB),
              lambdaB = 2*lambdaB;
          end
          %>>>>>>>>>>>>>>>>>>>>>   UPDATES FOR NEXT ITERATION   <<<<<<<<<<<<<<<<<<<<<< 
          % Update only if criterion has decreased
          if J_newB < JB,
              evecB = evec_newB;
              duB = up_newB-upB;                
              upB = up_newB;
              JB = J_newB;
              dwB = 1;
              mB = mB + 1;
              % Check if stop condition is satisfied 
              if sqrt(duB'*duB)<delta2B, break, end 
          end

          if lambdaB>1e3, break, end       
      end
      %>>>>>>>>>>>>>>>>>>>>>>>     Memilih nilai minimum terbaik     <<<<<<<<<<<<<<<<<<<<<<<<<
      if trB==1,
          Jmin_oldB = JB;
          upminB = upB;
      else
          if JB<Jmin_oldB,
              upminB = upB;
          end
      end
  end


  %>>>>>>>>>>>>>>>>>>>>>>     MENGHITUNG SINYAL KONTROL     <<<<<<<<<<<<<<<<<<<<<<
  eB = refB(iB) - yB;
 
  % Prediktif Kontroller
  uB=upminB(1);
  
  % Meyakinkan bahwa input sinyal kontrol berada dalam batasan
  if uB>ulim_maxB,
      uB=ulim_maxB;
  elseif uB<ulim_minB
      uB=ulim_minB;
  end
  upminB(1) = uB;

  
  %>>>>>>>>>>>>>>>>>>>       MENYIMPAN DATA DALAM DATA VEKTOR      <<<<<<<<<<<<<<<<<<<
  u_dataB(iB)       = uB;
  y_dataB(iB)       = yB;
  yhat_dataB(iB)    = yhat_vecB(tiyhB);
  t_dataB(iB)       = tB;
  J_dataB(iB)       = JB;
  e_dataB(iB)       = eB;

  %>>>>>>>>>>>>>>>>>>>>>>>>>>        MENG-UPDATE        <<<<<<<<<<<<<<<<<<<<<<<<<
  y_oldB     = [yB;y_oldB(1:end-1)];
  u_oldB     = [uB;u_oldB(1:end-1)];
  u_vecB(uviB) = upminB(upiB);
  u_vecB    = [u_vecB(2:length(u_vecB)) ; upminB(length(upB))];
  y_vecB    = [y_vecB(2:length(y_vecB)) ; yB];
  yhat_vecB(1:length(yhat_vecB)-1) = yhat_vecB(2:length(yhat_vecB));

  eB   = (((max(xbB)-min(xbB)))*(eB+1)/2)+min(xbB);%descaling
  IAEB = IAEB+abs(eB)*deltaB;

  %>>>>>>>>>>>>>>>>>>      Mencatat persentase simulasi terpenuhi      <<<<<<<<<<<<<<<<<
  progress(fighandle,floor(100*iB/samplesB));
end

%------------------------------------------------------------------------------
%------------------        >>>   AKHIR LOOP UTAMA   <<<       -----------------
%------------------------------------------------------------------------------
u_dataA=(((max(Reflu)-min(Reflu)))*(u_dataA+1)/2)+min(Reflu);
y_dataA=(((max(xdA)-min(xdA)))*(y_dataA+1)/2)+min(xdA);%descaling
ref_dataA=(((max(xdA)-min(xdA)))*(ref_dataA+1)/2)+min(xdA);%descaling

u_dataB=(((max(Qreboiler)-min(Qreboiler)))*(u_dataB+1)/2)+min(Qreboiler);
y_dataB=(((max(xbB)-min(xbB)))*(y_dataB+1)/2)+min(xbB);     % descaling
ref_dataB=(((max(xbB)-min(xbB)))*(ref_dataB+1)/2)+min(xbB); % descaling

%>>>>>>>>>>>>>>>>>>>>>>            MENGGAMBAR PLOTS           <<<<<<<<<<<<<<<<<<<<<<<
figure(1)
plot(ref_dataA,'b');
hold on
plot(y_dataA,'g');
grid
title('Produk Atas');
legend('biru=set point xd','hijau=output plant xd');

figure(2)
plot(u_dataA,'r');
grid
legend('biru=Reflux');

figure(3)
plot(ref_dataB,'b');
hold on
plot(y_dataB,'g');
grid;
title('Produk Bottom');
legend('biru=set point xb','hijau=output plant xb');

figure(4)
plot(u_dataB,'r');
grid;
legend('biru=Qreboiler');

c=fix(clock);
fprintf('\n\nTime is ended at %2i.%2i.%2i\n',c(4),c(5),c(6));
