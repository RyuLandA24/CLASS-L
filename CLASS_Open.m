%% README

% This code is a MATLAB version of a part of the
% the Chemistry Land-surface Atmosphere Soil Slab (CLASS) model
% de Arellano et al., (2015). Atmospheric boundary layer: Integrating Air
% Chemistry and Land Interactions. Cambridge University Press.

% The original codes are available at https://classmodel.github.io/  
% Kyoungho Ryu and Guido Salvucci converted the CLASS model to the MATLAB version.
% Department of Earth and Environment at Boston University

% We refer the CLASS textbook to readers for further detilas of the model
% and its origin.

clear
clc


%% Basic setting

pathsave = strcat('/projectnb/moisture/khr/OpenDataCh1/');

cd(pathsave)

% write a file name to save the results.
pathSave1 = append(pathsave, 'CLS_Results.mat');


% Parameters

pa_DZMIN = 300;

qscl=2500;  % Exponential scale for initial humidity profile, used in gammaq

beta = 0.2;  % Entrainment ratio for virtual heat [-]

HO= 100; % Initial ABL Height

pa_cefJq  = 0 ; % No Jump in init q

pa_CHISTAR = 1; % the parameter for moisture being trasponrted from the sub-cloud to cloud layers

% Advection
advtheta = 0;    % advection of heat [K s-1]
advq        = 0;    % advection of moisture [kg kg-1 s-1]


% Constants

julian = [1, 31   ; 32, 59  ; 60, 90  ; 91, 120 ; 121, 151; 152, 181; ...
           182, 212; 213, 243; 244, 273; 274, 304; 305, 334; 335, 365]; 
    
Lv  = 2.5e6;                 % heat of vaporization [J kg-1]
k    = 0.4;                    % Von Karman constant [-]
g    = 9.81;                  % gravity acceleration [m s-2]

         
% Define Constants
         BOLTZMAN = 1.380658e-23;
         AVOGADRO = .602214199e24;
         MD = 28.9644e-3;
         MV = 18.0153e-3;
         RV = (AVOGADRO)*(BOLTZMAN) / (MV);            % gas constant f or water vapor (J/(kg-K))
         RD = (AVOGADRO)*(BOLTZMAN) / (MD);            % gas constant for dry air (J/(kg-K))
         Cp = 7./2*(RD);                                                   % specific heat of air (J/(kg-K))
         CPD=Cp;
         cp=Cp;
         Rd=RD;
         Rv=RV;

         
% Constants for Lifiting Condensation Level, gds
        CPD         = cp;         % Specific heat of dry air [J kg-1 K-1]
        EPSILON  = Rd/Rv; % Unitless

% Time setting
         dt   = 1*60;                % seconds
         t     = 15;                    % for Hours, daytime          
         runtime  = t * 3600;                    % Total run time [s]. 
         tsteps = round(runtime ./ dt) ;  % Seconds
         
        tSave=  tsteps / 5; %  save every 5 minutes 
        incs=round(300/dt); % save every 5 minutes
        tRN_greater_zero = t;

%%  Simulation Inputs 

inputEF = (0:0.1:0.9);

Rn = [75 150 225];
gam = [30 50 70] ./ 10000;
RH = (30:10:90) ./ 100; 
T = [270 280 285 290 295 300];
pr = [80000 90000 100000];

cases = zeros( length(Rn)* length(gam) * length(RH) * length(T) * length(pr), 5);

ics = 0;
for i1 = 1: length(gam)
    for i2 = 1: length(RH)
        for i3 = 1:length(Rn)
            for i4 = 1:length(T)
                for i5= 1: length(pr)
ics = ics +1 ;
                    cases(ics, :) = [gam(i1) RH(i2) Rn(i3) T(i4) pr(i5)];

                end
            end
        end
    end
end

inputGamma = cases(:, 1);
inputRH = cases(:,2);
inputRn = cases(:, 3);
inputT = cases(:, 4);
inputPsfc = cases(:, 5);

    
   
%% Rn interpolation
funtime = (0: 5 :15*60);
fun = sin(pi * funtime / ( tRN_greater_zero*60) ) ; 

inputRn30 = inputRn .* fun;
[sz1, sz2] = size( inputRn30);

xt = funtime;
xitp = (0: 1: tRN_greater_zero * 60); % 15 hrs 

LE = zeros( [ size(inputRn30) , 1, length(inputEF)]);
SH = zeros( [ size(inputRn30) , 1, length(inputEF)]);
    for ief = 1: length(inputEF)
        LE(:, :, :, ief) = ( inputEF(ief) .* inputRn30 .* inputRn ./ ( sum(inputRn30, 2) / ( sz2 -1 ) )  ) ...
            .* 24 ./ tRN_greater_zero ;
        SH(:, :, :, ief) = ( (1 - inputEF(ief) ) .* inputRn30 .* inputRn ./ ( sum(inputRn30, 2) / ( sz2-1 ) ) ) ...
            .* 24 ./ tRN_greater_zero ; 
    end


inputRn_perm = permute( inputRn30, [2, 1]);
inputRnD = interp1(xt, inputRn_perm, xitp, 'linear', 'extrap');
inputRn = permute ( inputRnD, [2, 1] );

LE_perm = permute( LE, [2, 1, 3, 4]);
LED = interp1 (xt, LE_perm, xitp, 'linear', 'extrap');
LE = permute( LED, [2, 1, 3, 4]);

SH_perm = permute( SH, [2, 1, 3, 4]);
SHD = interp1 (xt, SH_perm, xitp, 'linear', 'extrap');
SH = permute( SHD, [2, 1, 3, 4]);


%% Preallocation
[D1, D2, D3, D4] = size(LE);

% Results of diurnal evolutions
CLS_h         = zeros(D1, tSave, D3, D4 );
CLS_TABL   = zeros(D1, tSave, D3, D4 );
CLS_qABL   = zeros(D1, tSave, D3, D4 );
CLS_JumpT = zeros(D1, tSave, D3, D4 );
CLS_Jumpq = zeros(D1, tSave, D3, D4 );

% Variables related to cloud core fraction
CLS_acbulk = zeros(D1, tSave, D3, D4 );
CLS_M         = zeros(D1, tSave, D3, D4 );
CLS_wqM    = zeros(D1, tSave, D3, D4 );
CLS_wqMA  = zeros(D1, tSave, D3, D4 );
CLS_wqMB  = zeros(D1, tSave, D3, D4 );
CLS_qvar_h  = zeros(D1, tSave, D3, D4 );

% Variables at h
CLS_qsat_h  = zeros(D1, tSave, D3, D4 );
CLS_T_h       = zeros(D1, tSave, D3, D4 );
CLS_wTv       = zeros(D1, tSave, D3, D4 );

% Variables related to the ML
CLS_wstar   = zeros(D1, tSave, D3, D4 );
CLS_we        = zeros(D1, tSave, D3, D4 );
CLS_wqe      = zeros(D1, tSave, D3, D4 );
CLS_wTe      = zeros(D1, tSave, D3, D4 );
CLS_dzh      = zeros(D1, tSave, D3, D4 );
CLS_Q0       = zeros(D1, tSave, D3, D4 );
CLS_DQ      = zeros(D1, tSave, D3, D4 );
CLS_BF       = zeros(D1, tSave, D3, D4 );

% LCL
CLS_lclbulk = zeros(D1, tSave, D3, D4 );

% Variables related to the surface
CLS_Rn   = zeros(D1, tSave, D3, D4 );

CLS_wq   = zeros(D1, tSave, D3, D4 );
CLS_wT   = zeros(D1, tSave, D3, D4 );

CLS_wqF   = zeros(D1, tSave, D3, D4 );
CLS_wTF   = zeros(D1, tSave, D3, D4 );


%% MIXED LAYER
% Pre-allocation of initial variables

% Intergrate mixed-layer equations
h      = zeros(D1, D3, D4);
TABL = zeros(D1, D3, D4);
GamT = zeros(D1, D3, D4);
PsfcA = zeros(D1, D3, D4);
RHA = zeros(D1, D3, D4);

% cumulus parameterization
dz_h  = zeros(D1, D3, D4);
ac      = zeros(D1, D3, D4);
M       = zeros( D1,D3, D4);

wqe_lasttime = zeros( size ( dz_h) );
wqM_lasttime = zeros( size( dz_h) );
wTM_lasttime = zeros( size( dz_h) );

wqMA = zeros( size (dz_h) );
wqMB = zeros( size (wqMA) );

htend = zeros( size (wqMA) );
Ttend = zeros( size (wqMA) );
qtend = zeros( size (wqMA) );

JumpTtend = zeros( size (wqMA) );
JumpTtendTOT = zeros( size (wqMA) );

Jumpqtend = zeros( size (wqMA) );
JumpqtendTOT = zeros( size (wqMA) );
        
'Assign initial values' 

h(:, :, :) = HO; % initial ABL height [m]

dz_h(:, :, :) = 150; % Transition layer thickness  from dry to moist convection [m]

for ief = 1 :D4

TABL(:, :, ief)  = inputT;
GamT(:, :, ief) = inputGamma;
PsfcA(:, :, ief) = inputPsfc;
RHA(:, :, ief) = inputRH;

end

rho =  PsfcA  ./ (Rd .* TABL );

% Humidity
qsat_a = calcqsatB( TABL, PsfcA); %[kg/kg]
qABL    =  RHA .* qsat_a; % Initial mixed-layer specific humidity [kg/kg]

% Jumps
JumpT = (GamT .* h .* beta) ./ (1 + 2 * beta);  % initial temperature jump at h= the entrainment zone. the default = 1 [K]

Jumpq =  pa_cefJq .* qABL; % initial specific humidity jump at h [kg/kg]



'Time Evolution'
tind=0; % This index (tind..time index) is used later to save 1/2 hour averages
for i1 = 1 : tsteps
% Surface Turbulent Flux 
SHa = reshape( SH(:, i1, :, :), [ D1, D3, D4]);
LEa = reshape( LE(:, i1, :, :), [ D1, D3, D4]);

% Kinematic Surface Fluxes 
 wT   =  SHa  ./ (rho .* cp);  % Heat flux [Km/s]
 wq   = LEa ./ (rho .* Lv);  % Moisture flux [kg/kgm/s]

% Virtual Temperatures
wTv       = wT  + 0.61 .* TABL .* wq; % [Km/s]
Tv          = TABL  + 0.61 .* TABL .* qABL; % [K]
JumpTv  = ( (TABL + JumpT ) .* (1 + 0.61 .* ( qABL + Jumpq ) ) ) ...
                                - ( TABL .* ( 1+ 0.61 .* qABL ) ); % [K]


%% Lifiting Condensation Level

T_lcl = 2840 ./ (3.5 .* log( TABL ) - log( PsfcA ./ 1000 .* qABL ./ (EPSILON + qABL ) ) - 7.108) + 55;
P_lcl  = PsfcA .* ( T_lcl ./ TABL ).^3.5;
lcl = CPD  .* TABL ./ g .* ( 1 - ( P_lcl ./ PsfcA ).^(Rd ./ CPD)); % [m]
lclbulk = max( lcl, HO); % To prevent imaginary value.
qsat_lcl = calcqsatB( T_lcl , P_lcl );
RH_lclbulk = qABL ./ qsat_lcl;

P_h      = PsfcA - ( rho .* g .* h );   
T_h      = TABL - ( (g/cp) .* h ); 
qsat_h = calcqsatB( T_h, P_h);
RH_h   = qABL ./ qsat_h;

%% Mixed Layer

% Convective velocity scale, w* [m/s].
% Deardorff, 1970, doi:10.1175/1520-0469(1970)027<1211:CVATSF>2.0.CO;2
BF  =  g .* ( wTv ./ Tv); 
wstar = ( BF .* h  ).^(1/3);  %[m/s]
    wstar( wTv <= 0 ) = 10.^(-6);
    wstar( wstar < 1e-6 ) = 1e-6;  % Not divide by zero later which causes NAN and crash
   
% Virtual heat entrainment flux
wTve = -beta .* wTv; % [Km/s]

% The entrainment velocity [m/s]
% Typical values are often in the range of 0.01 m/s to 0.20 m/s.
we = -wTve ./ JumpTv ; % [m/s]
    we( we < 0 ) = 0; % No BL sinking if wtheta <0

% Entrainment fluxes
wTe   = -(we) .* JumpT; % [ Km/s] 
wqe  = -(we) .* Jumpq;  % [kg/kgm/s]  

%% Cumulus Clouds 
% KHR & GDS
            
  % gamma q calcs
            if i1==1
                Q0 = qABL;
                Theta0 = TABL;                
                gammatheta = inputGamma;                
            end
                        
            % note, if q/qo=exp(-z/zo) then dqdz at surface=qo/z
            % this is confirmed in gfdl model with mean scale of 2500
            % and very small spatial standard deviation of about 100 M
            
            if i1==1
                gamsurf = -Q0  ./ qscl;
            end
           
            Gammaq = gamsurf .* exp(-( h-HO ) ./ qscl);
   
            
            % note that even though altering jumpq to account for moisture
            % above h, if wqM linear from value at h to zero at ZD, then
            % dqdt is constant with height in this zone, and thus the gammaq
            % does not actually change...stays simply dependent on z through above
            
            %%%%%%%%%%%%%%%%%%%%
                     
            
            %%%%%%%%%%%%%
            % variance calcs
                    
           DQ =-Jumpq ; % 
     
           wqe_lasttimeH = wqe_lasttime;
           wqM_lasttimeH = wqM_lasttime;
           
      
            qvar_h = ( wqe_lasttimeH + wqM_lasttimeH ) .* ...
                    ( DQ ./ dz_h ) .* ( h ./ wstar  ); % [Unitless]
    
  
  %%%%%%%% ac calcs
            
            
  qvar_hH = qvar_h;

  qABLH = qABL;

  qsat_hH = qsat_h;
  qvar_hH( qvar_hH <=  1e-10 )= 1e-10;       
  sigq= ( qvar_hH.^0.5);

    acbulk = normcdf( qABLH- qsat_hH, 0 , sigq);
            
    M = acbulk .* wstar; % [m/s]

    wqMA = pa_CHISTAR .* M .* sigq; % [kg/kgm/s]
    wqMB = -M .* Jumpq; % [kg/kgm/s]


    htend  = we - M; % [m/s]

    Ttend   =  ( wT - wTe ) ./ h + advtheta ;  % [ K/s]
    qtend   =  ( wq - wqe - wqMA  ) ./ h + advq ; % [kg/kg/s]
    
    JumpTtend  = inputGamma .* ( htend ) - Ttend ;  % [ K/s]
    Jumpqtend  = Gammaq .* ( htend ) - qtend ; % [ kg/kg/s]

    JumpqtendTOT = Jumpqtend ;
    JumpTtendTOT = JumpTtend;


 wqM_lasttime = wqMA;
 wqe_lasttime = wqe;    
 wTe_lasttime = wTe;


 JumpqH = Jumpq;
 JumpTH = JumpT;      

dztend = ( ( lclbulk- h ) - dz_h ) ./ 7200; % [m]

% Omit windtend

%% Save data every 5 minutes  

if mod( (i1-1), incs)== 0 

    tind=tind+1;


% Diurnal Evolutions
CLS_h     (:,tind,:,:)    = h;
CLS_TABL (:,tind,:, :)  = TABL;
CLS_qABL    (:,tind, :, :)  = qABL;
CLS_Jumpq (:,tind, :, :)  = Jumpq ;
CLS_JumpT (:,tind, :, :)  = JumpT ;

% Cloud Mass Flux
CLS_acbulk(:, tind, :, :)      = acbulk ;

CLS_M         (:,tind,:, :)  = M ;

CLS_wqM   (:,tind, :, :)  = wqMA + wqMB ; %wqM ;

CLS_wqMA   (:,tind, :, :)  = wqMA ;
CLS_wqMB   (:,tind, :, :)  = wqMB ;

CLS_qvar_h(:, tind, :, : )  = qvar_h;

% The top of the ML
CLS_qsat_h(:,tind, :, :)  = qsat_h;
CLS_T_h      (:,tind, :, :)  = T_h;
CLS_wTv (:, tind, :, :) = wTv ;

% Mixed Layer
CLS_wstar (:,tind, :, :) =  wstar;
CLS_we(:,tind,:, :)      =  we ;
CLS_wqe(:,tind, :, :)   = wqe  ;
CLS_wTe(:,tind, :, :)   = wTe ;

CLS_dzh(:,tind,:, :) = dz_h ;
CLS_Q0(:,tind,:, :) = Q0 ;
CLS_DQ(:,tind,:, :) = DQ ;
CLS_BF (:,tind,:, :) = BF;

% LCL
CLS_lclbulk (:,tind,:, :) = lclbulk;

% The Surface
for ief = 1 : 10
CLS_Rn(:,tind,:, ief)  = inputRn(:,i1) ;
end
CLS_wq(:,tind,:, :)   =  wq ;
CLS_wT(:,tind,:, :)   =  wT ;

CLS_wqF(:, tind, :, :) = wq .* rho ;
CLS_wTF(:, tind, :, :) = wT .* rho ;

end



%% intergrate mixed layer
h        = real( h + dt .* htend ); % [m]
h        = max( h,  HO); % to prevent h from being a weired value.

TABL  = real( TABL + dt .* Ttend ); % [K]
TABL  = max( TABL , 100); % to prevent TABL from being a weired value.

qABL   = real( qABL + dt .* qtend ); %[ kg/kg ]
qABL   = max ( qABL, 0.000001); % to prevent qABL from being a weired value.

JumpT = real( JumpT + dt .* JumpTtendTOT ); % [K]
Jumpq = real( Jumpq + dt .* JumpqtendTOT ); % [kg/kg]  
 
dz_h    = real( dz_h + dt .* dztend ); % [m]
dz_h    = max( pa_DZMIN, dz_h ); % keep dz_h positive and non-zero values

end % Time evolution


%% Profile
% KHR and GDS

HO =100;
ZS = (0:30:8100);

[sz1, sz2, ~ ,sz4] = size( CLS_h );

slope = +( (CLS_wqe + CLS_wqMA) - CLS_wq ) ./ (CLS_h - HO);

qf1 = zeros( length(ZS), sz1, sz2, sz4 );
qf2 = zeros( size ( qf1) );

LSTORE = 1000000;

for ii = 1: length(ZS)

    ii

    if ZS(ii) > HO

        LINZ = 1 - ( ZS(ii) - CLS_h) ./ LSTORE;

        delta = CLS_h > ZS(ii); % find below h


        qf1(ii, :, :, :) = ( CLS_wq + slope .* (ZS(ii)-HO) ) .* delta ;
        qf2(ii, :, :, :) = (CLS_wqMA + CLS_wqMB) .* LINZ .* (1-delta);

    else

         qf1(ii, :, :, :) = CLS_wq;

    end

end

 profqf = qf1 + qf2 ;
 profqf( profqf < 0 ) = 0;

save(pathSave1, 'CLS_acbulk', 'CLS_h', 'CLS_Jumpq', 'CLS_TABL', 'CLS_qABL', ...
    'CLS_wq',  'CLS_wqe', 'CLS_wqMA', 'CLS_wqMB', 'CLS_wqF', 'CLS_Rn', ...
    'profqf', '-v7.3');




