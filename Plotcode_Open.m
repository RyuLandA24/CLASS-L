% This code produces plots shown the paper, 
% Ryu, K and G. Salvucci (2023), 
% Modifications to the CLASS Boundary Layer Model for Improved Interaction between the Mixed Layer and Clouds
% sumbitted to Journal of Advances inModeling Earth Systems.

% To get the plots correctly, please down the heat map code 
% Jonathan C. Lansey provides available at 
% https://www.mathworks.com/matlabcentral/fileexchange/45325-efficient-2d-histogram-no-toolboxes-needed

clear
clc
% Please change the path accordingly.
cd /projectnb/moisture/khr/OpenDataCh1/
addpath /project/moisture/khr/Functions/
%%
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

EFA = (0:0.1:0.9);
caseA = zeros( 11340, 6);
i = 0 ;
for i1 = 1: length(gam)
    for i2 = 1: length(RH)
        for i3 = 1:length(Rn)
            for i4 = 1:length(T)
                for i5= 1: length(pr)
                    for i6 = 1: 10

                    i = i +1;
                    caseA(i, :) = [gam(i1) RH(i2) Rn(i3) T(i4) pr(i5) EFA(i6)];
                    
                    end
                end
            end
        end
    end
end



%%
% Best Match
load('CLASS-B.mat')

% CLASS
COR = strcat('CLS_Results.mat');    

load(COR, 'CLS_h')
    ChOR = CLS_h;
    clear CLS_h

% CLASS-L
CReg = strcat('CLS_L_Results.mat');    
load(CReg, 'CLS_h')
    ChR = CLS_h;
    clear CLS_h

% LES Humidity Flux
load('LES_wqA.mat')
LwqDM = squeeze( mean( wqA, 3, 'omitnan') );

    % Dimension 
    Lwqf = zeros( 11340, 271, 181);
    CwqBDMA = zeros(11340, 271);
    
    i = 0;
    for ics = 1: 1134
        for ief = 1: 10
            i = i + 1 ;
            CwqBDMA(i, :) = CwqbestA(ics, :, ief);
            Lwqf(i, :, :) = wqA(ics, :, :, ief);
    
        end
    end

    LwqDMA = mean( Lwqf, 3, 'omitnan');

% CLASS
load(COR, 'profqf',  'profwqMA', 'profwqMB', 'profwqSE')
    CwqORDM = squeeze( mean( profqf, 3, 'omitnan'));
    CwqORDM = permute( CwqORDM, [2, 1, 3]);
    CprofOR = permute (profqf, [2, 1, 3, 4]);

    CprofAOR = profwqMA;
    CprofBOR = profwqMB;
    CprofSEOR = profwqSE;

    % Dimension
    CwqfO = zeros(11340, 271, 179);
    i = 0;
    for ics = 1: 1134
        for ief = 1: 10
    
            i = i+1;
            CwqfO(i, :, : ) = CprofOR(ics, :, 2:180, ief);
            
        end
    end


    test = zeros( 11340, 179);
    for i = 1 : 11340
        
    
        for it = 1 : 179

             test1 =  find( isnan ( CwqfO(i, :, it) ) );

            if isempty(test1) == 1
                 test(i, it ) = 1; % No nan
            else 
                test(i, it ) = 0; % Nan;
            end

    
        end 
    
    end

    test2 = squeeze( sum(test, 2) );
    test3 = find ( test2 < 179);

    CwqORDMA = mean( CwqfO, 3, 'omitnan');

        clear profqf

% CLASS-L
load(CReg, 'profqf', 'profwqMA', 'profwqMB', 'profwqSE', 'LSTORE')
    CwqRDM = squeeze( mean ( profqf, 3, 'omitnan') );
    CwqRDM = permute(CwqRDM, [2, 1, 3]);
    CprofR = permute(profqf, [2, 1, 3, 4]);
    
    CwqfR = zeros(11340, 271, 179);
    i = 0;
    for ics = 1: 1134
        for ief = 1: 10
    
            i = i+1;
            CwqfR( i, :, : ) = CprofR( ics, :, 2:180, ief);
            
    
        end
    end
    
    CwqRDMA = mean ( CwqfR, 3, 'omitnan') ;

% Remove nan cases from the CLASS model
LwqDMA(test3, :  ) = [];
CwqBDMA(test3, :) = [];
CwqRDMA(test3, :) = [];
CwqORDMA(test3, :) = [];


% Humidity
load('LES_qA.mat')
% CLASS
    load(COR, 'CLS_qABL', 'CLS_Jumpq');
        CqOR = CLS_qABL;
        CJqOR = CLS_Jumpq;
        clear CLS_qABL CLS_Jumpq
% CLASS-L
    load(CReg, 'CLS_qABL', 'CLS_Jumpq')
        CqR = CLS_qABL;
        CJqR = CLS_Jumpq;


% ac
    load('LES_acA.mat')
        accore = acareaA;
%CLASS
    load(COR, 'CLS_acbulk')
        CacOR = CLS_acbulk;
% CLASS-L
    load(CReg, 'CLS_accore');
        CacR = CLS_accore;

    
%% Settings
lwd = 3;
fszPT = 24; % Poster
fszPA = 8; % Paper

Zflux = (0 : 30 : 8100);
Zvar = (15: 30 : 8100);

Lcolor = "#77AC30";
Ccolor = "#0072BD";
Mcolor = "#D95319";

Pst1 = [10, -5, 15, 12];
Pst2 = [10, -8, 15 18.1]; % 1 or -8
Pst3 = [10 -10 19 23]; % 1 or - 10
Pst3h = [10 -10 19 15];

xticksH = (0: 36: 180);
xticksH(1) =1;


%% Figure1. Motivation
fsz = fszPA;

ics = 733;
ief = 5;

figure('units', 'centimeters', 'Position', Pst2);

t = tiledlayout( 1, 2, 'TileSpacing','compact');
xlabel(t, 'Humidity Flux (10^-^5kg/kgm/s)', 'FontSize', fsz)
ylabel(t, 'Height (km)', 'FontSize', fsz)

    nexttile
        plot(CwqORDM(825, :, 1), Zflux, 'lineWidth', lwd); hold on       
        plot(LwqDM(825, :, 1), Zflux, 'Color', Lcolor, 'lineWidth', lwd); 
            ylim([0 4000])
            yticks((0:1000:4000))
            yticklabels({'0', '1', '2', '3', '4'})    
            xticks([0 1e-5])
            xticklabels({'0', '1'})
            legend('CLASS', 'LES')
            
        title('(a)')
    set(gca, 'FontSize', fsz, 'TitleHorizontalAlignment', 'left')

    nexttile
       plot(CwqORDM(ics, :, ief), Zflux, 'lineWidth', lwd); hold on 
       plot(LwqDM(ics, :, ief), Zflux, 'Color', Lcolor, 'lineWidth', lwd);
            ylim([0 4000])
            yticks((0:1000:4000))
            yticklabels({'0', '1', '2', '3', '4'})    
            xticks([0 5e-5 10e-5 15e-5 20e-5])
            xticklabels({'0', '5', '10', '15', '20'})
            yticklabels({'0', '1', '2', '3', '4', '5', '6', '7', '8'})
            legend('CLASS', 'LES', 'location', 'northwest')
       title('(b)')
    set(gca, 'FontSize', fsz, 'TitleHorizontalAlignment', 'left')


%% Figure 2. Moistening above the mixed layer top via L. 
fsz = 8;

ics = 733;
ief = 5;

it = 60;

ChO = ChOR(ics, it, :, ief); % CLASS
ChM = ChR(ics, it, :, ief); % CLASS-L
LM = LSTORE(ics); % LES

figure('units', 'centimeters', 'Position', Pst3h);
t = tiledlayout( 1, 3 );
ylabel(t, 'Height (km)', 'FontSize', fsz)
xlabel(t, 'Humidity Flux (10^-^5kg/kgm/s)', 'FontSize', fsz)

nexttile
    plot(CprofOR(ics, :, it, ief), Zflux, ...
        'lineWidth', lwd, 'Color', Ccolor); hold on
    plot( CprofSEOR(:, ics, it, ief), Zflux, 'lineWidth', lwd, 'Color', "#EDB120"); hold on
    plot ( CprofAOR(:, ics, it, ief), Zflux, 'lineWidth', lwd, 'Color', "#7E2F8E" ); hold on
    plot( CprofBOR(:, ics, it, ief), Zflux, 'lineWidth', lwd, 'Color', "#27590F"); hold on
    legend('CLASS', 'SFC+ENT', 'Mq', 'M_J_q',...
         'location', 'northwest')
        ylim([0 4000])
        yticks([ChM])
        yticklabels({'h'})
        xticks([0 2e-4 4e-4])
        xticklabels({'0', '20', '40'})
set(gca, 'FontSize', fsz, 'TitleHorizontalAlignment', 'left')

% CLASS-L
nexttile
    plot( CprofR( ics, :, it, ief), Zflux, 'lineWidth', lwd, 'Color', Mcolor); hold on
    plot( profwqSE(:, ics, it, ief), Zflux, 'lineWidth', lwd, 'Color', "#EDB120"); hold on
    plot ( profwqMA(:, ics, it, ief), Zflux, 'lineWidth', lwd, 'Color', "#7E2F8E" ); hold on
    plot( profwqMB(:, ics, it, ief), Zflux, 'lineWidth', lwd, 'Color', "#27590F"); hold on
    ylim([0 4000])
        yticks([ChM ChM+LM])
        yticklabels({'h', 'h + L'})
        xticks([0 2e-4 4e-4])
        xticklabels({'0', '20', '40'})
        legend('CLASS-L', 'SFC+ENT', 'M_q', 'M_J_q', ...
         'location', 'northwest')
set(gca, 'FontSize', fsz, 'TitleHorizontalAlignment', 'left')

% Daily Average Profile
    nexttile
       plot(CwqORDM(ics, :, ief), Zflux, 'Color', Ccolor, 'lineWidth', lwd); hold on % CLASS
       plot(CwqbestA(ics, :, ief), Zflux, '--', 'Color', Mcolor, 'lineWidth', lwd ); hold on % Best
       plot(LwqDM(ics, :, ief), Zflux, 'Color', Lcolor, 'lineWidth', lwd);
          ylim([0 4000])
          xticks([0 1e-4 2e-4])
          xticklabels({'0', '10', '20'})
          yticks([0 1000 2000 3000 4000])
          yticklabels({'0', '1', '2', '3', '4'})
          legend('CLASS', 'Best Fit L', 'LES', 'Location','northwest')
          set(gca, 'FontSize', fsz, 'YAxisLocation', 'right', 'TitleHorizontalAlignment', 'left')





%% Figure 3. The best match L vs The regression-based L

load('activeacidx.mat')

% Indexes for the PBLW
idxPBLW = zeros(1134, 181, 10);

    for ics = 1 : 1134
        for ief =  1 : 10

            for id = 1 : 181
            PBLW = find ( activeacidx(ics, :, id, ief) > 0 );

                if isempty (PBLW) == 1
    
                    idxPBLW(ics, id, ief) = 0;
    
                else
                
                    idxPBLW(ics, id, ief) = 1;
    
                end

            end
        
        
        end
    end

idxPBLW2 = squeeze( sum( idxPBLW, 2) );
idxPBLW2( idxPBLW2 >= 2 ) = 1;

ARn = zeros( length(cases(:,3 )), 10);
ARH = zeros( length(cases(:,3 )), 10);
AGam = zeros( length(cases(:,3 )), 10);
AEF =zeros( length(cases(:,3 )), 10);
AT = zeros ( length(cases(:,3 )), 10);
AP = zeros( size(AT));

for ief = 1: 10
    ARn(:, ief) = cases(:, 3);
    ARH(:, ief) = cases(:, 2);
    AGam(:, ief) = cases(:, 1);
    AT(:, ief) = cases(:, 4);
    AP(:, ief) = cases(:, 5);
end

EF = (0: 0.1 : 0.9);
for ics = 1: length(cases(: , 3))
    AEF(ics, :) = EF;
end

EFgam = AGam;
EFP = AP;
EFRH = ARH;
EFRn = ARn;
EFT = AT;
EFEF = AEF;
EFLbest = Lbest;

EFLbest( idxPBLW2 == 0 ) = 0;
EFLbest( idxPBLW2 == 0 ) = nan;

EFgam(idxPBLW2 ==0 )= nan;
EFP(idxPBLW2 ==0 )= nan;
EFRH(idxPBLW2 ==0 )= nan;
EFRn(idxPBLW2 ==0 )= nan;
EFT(idxPBLW2 ==0 )= nan;
EFEF(idxPBLW2 ==0 )= nan;

X = [log(EFRn(:)), log(EFgam(:)), log(EFRH(:)),  log(EFT(:)) ];
Y = log(EFLbest(:));
    mdl = fitlm ( X, Y);
    cfs = table2array( mdl.Coefficients(:, 1) );

Lmd = exp(cfs(1)) .* EFRn(:).^(cfs(2)) .* EFgam(:).^(cfs(3)) ...
    .* EFRH(:).^(cfs(4)) .*  EFT(:).^(cfs(5));

figure('units', 'centimeters', 'Position', Pst1);
    ndhist(Lmd(:), EFLbest(:), 'normalizex', 'filt', 8); hold on
    line( EFLbest(:), EFLbest(:), 'Color', 'k', 'lineWidth', 3)
    colorbar
    xlabel('The Regression-Based L (km)')
        xticks([0 1000 2000 3000 4000])
        xticklabels({'0', '1', '2', '3', '4'})
    ylabel('The Best Fit L (km)')
        yticks([0 1000 2000 3000 4000])
        yticklabels({'0', '1', '2', '3', '4'})
    set(gca, 'FontSize', 10)
    
    mdl = fitlm( Lmd(:), EFLbest(:))



%% Figure 4. Humidity Jump and Humidity Dynamic throughout the daytime

qstar = zeros( 1134, 181, 10);
qstaridx = zeros( 1134, 181, 10);

for ics = 1 : 1134
    for ief = 1 : 10

        for it = 2 : 181
        qmxZ = qA(ics, :, it, ief) - qA(ics, :, 1, ief)    ;
        
        [qmxZ1 , qmxZidx] = max( qmxZ);
        qstar(ics, it, ief) = qA(ics, qmxZidx, it, ief);
        qstaridx(ics, it, ief) = qmxZidx;

        end

    end
end

qstar(:, 1, : )= qA(:, 4, 1, : );
qstaridx(:, 1, :) = 4;

qzd = qstar(:, :, :) - qstar(:, 1, :);

CqOd = CqOR(:, :, :, :) - CqOR(:, 1, :, :);
CqMd = CqR(:, :, :, :) - CqR(:, 1, :, :);
CqBd = Cqbest(:, :, :, :) - Cqbest(:, 1, :, :);

qzdDM = squeeze( mean( qzd, 2, 'omitnan'));
CqOdDM = squeeze( mean( CqOd, 2, 'omitnan') );
CqMdDM = squeeze( mean( CqMd, 2, 'omitnan') );
CqBdDM = squeeze( mean( CqBd, 2, 'omitnan'));


test4 = zeros( 10,1);
for i =  1: 10
cs1 = find ( cases(:, 1) == caseA(test3(i), 1) );
cs2 = find( cases(:, 2) == caseA(test3(i), 2));
cs3 = find ( cases(:, 3) == caseA(test3(i), 3) );
cs4 = find( cases(:, 4) == caseA(test3(i), 4));
cs5 = find ( cases(:, 5) == caseA(test3(i), 5) );

test4(i) = intersect( intersect( intersect( intersect( cs1, cs2), cs3), cs4), cs5);
end

for i = 1:10
qzdDM(test4(i), 10) = nan;
CqOdDM(test4(i), 10) = nan;
CqMdDM(test4(i), 10) = nan;
CqBdDM(test4(i), 10) = nan;
end


%
ics = 625;
ief = 5;
fsz = fszPA;

scmatch = 10^3; % scale match

figure('units', 'centimeters', 'Position', Pst3);
    tiledlayout(3, 3)

    nexttile(1, [1 3])
        plot(qzd(ics, :, ief) .* scmatch, 'Color', Lcolor, 'lineWidth', lwd) ; hold on
        plot( squeeze( CqOd(ics, :, :, ief)) .* scmatch, 'Color', Ccolor, 'lineWidth', lwd ); hold on
        plot( squeeze( CqBd(ics, :, :, ief)) .* scmatch, 'Color', Mcolor, 'lineStyle', '--', 'lineWidth', lwd); hold on
        plot( squeeze( CqMd(ics, :, :, ief)) .* scmatch, 'Color', Mcolor, 'lineWidth', lwd )
            xlim([1 180])
            xticks( xticksH);
            xticklabels({'0', '3', '6', '9', '12', '15'})
            ylabel('Humidity Difference (g/kg)')
            legend('LES', 'CLASS', 'Best Fit L', 'Regression L', 'Location', 'southwest')
    set(gca, 'FontSize', fsz)

     nexttile(4, [1 3])
         plot( squeeze( CJqOR(ics, :, :, ief)) .* scmatch, 'Color', Ccolor, 'lineWidth', lwd); hold on
         plot( squeeze( CJqbest(ics, :, :, ief)) .* scmatch, 'Color', Mcolor, 'lineStyle', '--', 'lineWidth', lwd); hold on
         plot( squeeze( CJqR(ics, :, :, ief)) .* scmatch, 'Color', Mcolor, 'lineWidth', lwd)
            xlim([1 180])
            xticks( xticksH);
            xticklabels({'0', '3', '6', '9', '12', '15'})
            xlabel('Time [Hours]')            
            yticks([-4e-3 -3e-3 -2e-3 -1e-3 0].* scmatch)
            ylabel('Humidity Jump (g/kg)')
            legend('CLASS', 'Best Fit L', 'Regression L', 'Location', 'southwest')
    set(gca, 'FontSize', fsz)

           
     nexttile
        ndhist(CqOdDM.* scmatch, qzdDM.* scmatch, 'normalizex', 'filt', 8); hold on
        line(qzdDM.* scmatch, qzdDM.* scmatch, 'Color', 'k', 'lineWidth', lwd)
        colorbar
            xlabel({'(g/kg)', 'CLASS' })
            ylabel({'LES', '(g/kg)'})
    set(gca, 'FontSize', fsz)

    nexttile
        ndhist(CqBdDM.* scmatch, qzdDM.* scmatch, 'normalizex', 'filt', 8); hold on
        line(qzdDM.* scmatch, qzdDM.* scmatch, 'Color', 'k', 'lineWidth', lwd)
        colorbar
            xlabel({'(g/kg)', 'Best Fit L'})
            yticklabels([])
    set(gca, 'FontSize', fsz)

    nexttile
        ndhist(CqMdDM.* scmatch, qzdDM.* scmatch, 'normalizex', 'filt', 8); hold on
        line(qzdDM.* scmatch, qzdDM.* scmatch, 'Color', 'k', 'lineWidth', lwd)
        colorbar
            xlabel({'(g/kg)', 'Regression L'})
            yticklabels([])
    set(gca, 'FontSize', fsz)
       




%% Figure 5. Comparison of Humidity Fluxes

iz30 = 101; % 3km
iz20 = 68; % 2km 
iz10 = 35; % 1 km

L3A = squeeze( LwqDMA(:, iz30) );
L2A = squeeze( LwqDMA(:, iz20) );
L1A = squeeze( LwqDMA(:, iz10) );


%
CM3A = squeeze( CwqRDMA(:, iz30) );
CM2A = squeeze( CwqRDMA(:, iz20) );
CM1A = squeeze( CwqRDMA(:, iz10) );

%
CO3A = squeeze( CwqORDMA(:, iz30) );
CO2A = squeeze( CwqORDMA(:, iz20) );
CO1A = squeeze( CwqORDMA(:, iz10) );

%
CB3A = squeeze( CwqBDMA(:, iz30) );
CB2A = squeeze( CwqBDMA(:, iz20) );
CB1A = squeeze( CwqBDMA(:, iz10) );

fsz = fszPA;

figure('units', 'centimeters', 'Position', Pst3);
    t= tiledlayout( 3, 3, 'TileSpacing','compact');
 
    ylabel(t, 'LES (10^-^5kg/kgm/s)', 'FontSize', fsz)
    xlabel(t,  'Humidity Flux (10^-^5kg/kgm/s)', 'FontSize', fsz)

    nexttile
        ndhist(CO3A(:), L3A(:), 'normalizex', 'filt', 8); hold on
        line(L3A(:), L3A(:), 'lineWidth', lwd, 'Color', 'k');
        xlim([0 20e-5])
            xticks([0  10e-5 20e-5])
            xticklabels({'0',  '10',   '20'})
            yticks([ 0 5e-5])
            yticklabels({'0', '5'})
        title('(a)           CLASS')
        set(gca, 'YAxisLocation', 'right', 'FontSize', fsz, 'TitleHorizontalAlignment', 'left')

    nexttile
        ndhist(CB3A(:), L3A(:), 'normalizex', 'filt', 8);
        line(L3A(:), L3A(:), 'lineWidth', lwd, 'Color', 'k');
            xticks([0  10e-5])
            xticklabels({'0',  '10'})
            yticks([ 0 5e-5 ])
            yticklabels({'0', '5' })
            title('(b)         Best Fit L')
            set(gca, 'YAxisLocation', 'right', 'FontSize', fsz, 'TitleHorizontalAlignment', 'left')

       
    nexttile
        ndhist(CM3A(:), L3A(:), 'normalizex', 'filt', 8);
        line(L3A(:), L3A(:), 'lineWidth', lwd, 'Color', 'k');
            xticks([0 5e-5 ])
            xticklabels({'0', '5'})
            yticks([ 0 5e-5 ])
            yticklabels({'0', '5'})
         title('(c)       Regression L')
        set(gca, 'YAxisLocation', 'right', 'FontSize', fsz, 'TitleHorizontalAlignment', 'left')

% 2 km
    nexttile
        ndhist(CO2A(:), L2A(:), 'normalizex', 'filt', 8); hold on
        line(L2A(:), L2A(:), 'lineWidth', lwd, 'Color', 'k');
            xticks([0  10e-5  20e-5])
            xticklabels({'0',  '10',  '20'})
            yticks([0 10e-5])
            yticklabels({'0',  '10'})
         title('(d)')
    set(gca, 'YAxisLocation', 'right', 'FontSize', fsz, 'TitleHorizontalAlignment', 'left')

    nexttile
        ndhist(CB2A(:), L2A(:), 'normalizex', 'filt', 8); hold on
        line(L2A(:), L2A(:), 'lineWidth', lwd, 'Color', 'k');
            xticks([0  10e-5 ])
            xticklabels({'0', '10'})
            yticks([0 10e-5])
            yticklabels({'0',  '10'})
         title('(e)')
    set(gca, 'YAxisLocation', 'right', 'FontSize', fsz, 'TitleHorizontalAlignment', 'left')

    nexttile
        ndhist(CM2A(:), L2A(:), 'normalizex', 'filt', 8);
        line(L2A(:), L2A(:), 'lineWidth', lwd, 'Color', 'k');
        colorbar
        xticks([0  10e-5])
            xticklabels({'0', '10'})
            ylim([0 12e-5])
            yticks([0 10e-5])
            yticklabels({'0',  '10'})
            title('(f)')
    set(gca, 'YAxisLocation', 'right', 'FontSize', fsz, 'TitleHorizontalAlignment', 'left')
    
        
        % 1 km
    nexttile     
        ndhist(CO1A(:), L1A(:), 'normalizex', 'filt', 8); hold on
        line(L1A(:), L1A(:), 'lineWidth', lwd, 'Color', 'k');
        title('(g)')
            xticks([0 10e-5  ])
            xticklabels({'0', '10' })
            yticks([0 10e-5])
            yticklabels({'0',  '10'})
    set(gca, 'YAxisLocation', 'right', 'FontSize', fsz, 'TitleHorizontalAlignment', 'left')

    nexttile   
        ndhist(CB1A(:), L1A(:), 'normalizex', 'filt', 8); hold on
        line(L1A(:), L1A(:), 'lineWidth', lwd, 'Color', 'k');
         title('(h)')
             xticks([0 10e-5 ])
             xticklabels({'0', '10'})
             yticks([0 10e-5])
             yticklabels({'0', '10'})
    set(gca, 'YAxisLocation', 'right', 'FontSize', fsz, 'TitleHorizontalAlignment', 'left')

    nexttile
        ndhist(CM1A(:), L1A(:), 'normalizex', 'filt', 8); hold on
        line(L1A(:), L1A(:), 'lineWidth', lwd, 'Color', 'k');
            xticks([0  10e-5])
            xticklabels({'0',  '10'})
            yticks([0  10e-5 ])
            yticklabels({'0', '10'})
            title('(i)')
    set(gca, 'YAxisLocation', 'right', 'FontSize', fsz, 'TitleHorizontalAlignment', 'left')

    

%% Figure 6. Comparison of Cloud Core Fraction
        
Cvar = Chbest;

Lacmx = zeros(1134, 180, 10);

%
for ics = 1 : 1134
     for ief = 1 : 10


       for it = 1 : 180
                   
        Lcheck = find(  acareaA(ics, :, it ,ief) > 0 );
    

        if isempty(Lcheck) == 1
        Lacmx(ics, it, ief) = 0;
        else
        % active clouds exist at this time    
        Lcd = find ( diff( Lcheck ) > 1);   

        % Consider the gap between heights
            if isempty(Lcd) == 0
                Lcheck = Lcheck(1: Lcd);
            end
            
        Lcd1 = Lcheck(1);

        [~, Chidx ] = min( abs( Cvar(ics, it, :, ief) - Zvar) );
    

        if Lcd1 < Chidx == 1 & Chidx < Lcheck(end) == 1
        % ML top is between active clouds in the LES
       Lacmx1 = max( acareaA( ics, Lcd1 : Chidx, it, ief), [], 2)    ;
       Lacmx2 = max( acareaA( ics, Chidx : Lcheck(end), it, ief), [], 2)    ;

        Lacmx(ics, it, ief) = Lacmx1;

        end
 
       if Chidx == Lcd1 == 1
       Lacmx(ics, it, ief) = acareaA( ics, Lcd1, it, ief)    ;
       end

       if Chidx > Lcheck(end) == 1
       Lacmx(ics, it, ief) = max( ...
                    acareaA(ics, Lcheck, it, ief), [], 2);
        end

        if Chidx < Lcd1
            Lacmx(ics, it, ief) = 0;
        end


        end


       end

      
     end
end

    LacBmx = Lacmx;

% 
LBA = zeros( 11340, 180);
CBA = zeros( size( LBA ));
COA = zeros( size( LBA));
COR = zeros( size( LBA));

i = 0;
for ics = 1: 1134
    for ief  = 1: 10
i = i + 1;

    LBA(i, :) = LacBmx(ics, :, ief);
    CBA(i, :) = Cacbest(ics, :, :, ief);
    COA(i, :) = CacOR(ics, :, :, ief);
    COR(i, :) = CacR(ics, :, :, ief);

    end
end

LBA(test3, :) = [];
CBA(test3, :) = [];
COA(test3, :) = [];
COR(test3, :) = [];

LacmxA = mean( LBA, 1, 'omitnan');
CacBA = mean ( CBA, 1, 'omitnan');
CacOA = mean( COA, 1, 'omitnan');
CacRA = mean( COR, 1, 'omitnan');

%

ics = 164;
ief = 8;


f= figure('units', 'centimeters', 'Position', Pst2);
   t =  tiledlayout(2, 1, 'TileSpacing','compact');
    ylabel(t, 'Cloud Core Fraction', 'FontSize', fsz)
    xlabel(t, 'Time [Hour]', 'FontSize', fsz)
    Lsm1 = find ( LacBmx(ics, :, ief) > 0 );
    Lsm2 = smooth( LacBmx(ics, Lsm1, ief), 12);
    Lsm = LacBmx(ics, :, ief);
    Lsm(Lsm1) = Lsm2;

    COsm = smooth( CacOR(ics, :, : , ief), 12);
    CBsm = smooth( Cacbest(ics, :, :, ief), 12);
    CRsm = smooth( CacR(ics, :, :, ief), 12);

    COAsm  = smooth( CacOA, 12);
    CBAsm = smooth( CacBA, 12);
    CRAsm = smooth( CacRA, 12);
    LAsm = smooth( LacmxA, 12);

nexttile
    plot( Lsm, 'Color', Lcolor, 'lineWidth', lwd ); hold on
    plot( COsm, 'Color', Ccolor, 'lineWidth', lwd); hold on
    plot( CBsm, 'Color', Mcolor, 'lineStyle', '--', 'lineWidth', lwd); hold on
    plot( CRsm, 'Color', Mcolor, 'lineWidth', lwd)
         xlim([1 180])
         xticks( xticksH);
         xticklabels({'0', '3', '6', '9', '12', '15'})
         ylim([0 0.2])
         yticks([0 0.1 0.2])
        legend('LES', 'CLASS', 'Best Fit L', 'Regression L', 'location', 'northwest', 'FontSize', fsz)
set(gca, 'FontSize', fsz)

nexttile
    plot(LAsm, 'Color', Lcolor, 'lineWidth', lwd); hold on
    plot(COAsm, 'Color', Ccolor, 'lineWidth', lwd); hold on
    plot(CBAsm, 'Color', Mcolor, 'lineStyle', '--', 'lineWidth', lwd); hold on
    plot(CRAsm, 'Color', Mcolor, 'lineWidth', lwd); hold on
        xlim([1 180])
        xticks( xticksH);
        xticklabels({'0', '3', '6', '9', '12', '15'})
        ylim([0 0.1])
        yticks([0 0.1])
        legend('LES', 'CLASS', 'Best Fit L', 'Regression L', 'location', 'northwest')
set(gca, 'FontSize', fsz)

