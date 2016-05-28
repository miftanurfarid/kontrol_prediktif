% -------------------------------> npcinit_heat.M <------------------------------
% ----------      Switches       -----------
regtyA      ='npcA';            % Type kontroller (npc, pid, none)
reftyA      ='myrefA';          % Sinyal reference
simul       ='nnet';            % Obyek yang dikendalikan
if exist('simulink')~=5,
  simul     ='matlab';      
end

% ----------   Inisialisasi   -----------
TsA = 0.02;                     % Periode sampling = delta (dalam detik)
samplesA = 200;                 % Jumlah sample dalam simulasi
u_0A = 0;                       % Inisialisasi sinyal kontrol input
y_0A = 0;                       % Inisialisasi output
ulim_minA = -Inf;               % Minimum sinyal kontrol input
ulim_maxA = Inf;                % Maximum sinyal kontrol input

% ----- Neural Network Specification ------
nnfileA  = 'sisoA_n8_h4';       % Nama file model JST berisi NN, NetDef, W1, W2

% ---------- GPC initializations -----------
N1A = 1;                        % Horison prediksi minimum (harus sama dengan time delay)
N2A = 10;                       % Horison prediksi maksimum (>= nb)
NuA = 2;                        % Horison kontrol
rhoA = 0.9;                     % Faktor pengali perubahan sinyalkontrol


% -- Minimization algorithm initialzations -
maxiterA = 5;                   % Iterasi maksimum untuk menghitung u 
deltaA = 1e-4;                  % Norm-of-control-change-stop criterion
initvalA = '[upminA(NuA)]';

% ----------- Reference signal -------------
myrefA =[0.99*ones(100,1)];
% myrefA =[myrefA;0.991*ones(100,1)];
% myrefA =[myrefA;0.989*ones(400,1)];


% ----------      Switches       -----------
regtyB      ='npcB';            % Type kontroller (npc, pid, none)
reftyB      ='myrefB';          % Sinyal reference (siggener/<var. name>)
simul      ='nnet';             % Obyek yang dikendalikan (simulink/matlab/nnet)
if exist('simulink')~=5,
  simul      ='matlab';      
end


% ----------   Initializations   -----------
TsB = 0.02;
samplesB = 200;                 % Jumlah sample dalam simulasi
u_0B = 0;                       % Inisialisasi sinyal kontrol input
y_0B = 0;                       % Inisialisasi output
ulim_minB = -Inf;               % Minimum sinyal kontrol input
ulim_maxB = Inf;                % Maximum sinyal kontrol input

% ----- Neural Network Specification ------
nnfileB  = 'sisoB_n8_h4';       % Nama file model JST berisi NN, NetDef, W1, W2

% ---------- GPC initializations -----------
N1B = 1;                        % Horison prediksi minimal (harus sama dengan time delay)
N2B = 10;                       % Horison prediksi maksimum (>= nb)
NuB = 5;                        % Horison kontrol
rhoB = 0.1;                     % Faktor pengali perubahan sinyal kontrol


% -- Minimization algorithm initialzations -
maxiterB = 5;                   % Iterasi maksimum untuk menghitung u 
deltaB = 1e-4;                  % Norm-of-control-change-stop criterion
initvalB = '[upminB(NuB)]';


% ----------- Reference signal -------------
myrefB =[0.01*ones(100,1)];
myrefB =[myrefB;0.011*ones(400,1)];
% myrefB =[myrefB;0.01*ones(400,1)];
% myrefB =[myrefB;0.009*ones(400,1)];

myrefA=((2/(max(xdA)-min(xdA)))*(myrefA-min(xdA)))-1;		%Scaling
myrefB=((2/(max(xbB)-min(xbB)))*(myrefB-min(xbB)))-1;		%Scaling

