% -------------------------------------------------------------------------
% WeltraumFahrstuhl.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet verschiedene Parameter eines Space Elevators auf Erde 
% und Mars. 
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

%%
% Parameter
RE = 6370000;                       % Radius Erde m
ME = 5.97*10^24;                    % Masse Erde kg
rhoE = 5510 ;                       % Dichte Erde kg/m^3
RM =  3389000;                      % Radius Mars m
MM = 6.417*10^23;                   % Masse Mars kg
rhoM = 3344;                        % Dichte Mars kg/m^3
G = 6.671*10^-11;                   % Gravitationskonstante
omegaE = 2*pi/86125;                % Rotationsfrequenz Erde
omegaM = 2*pi/88620;                % Rotationsfrequenz Mars
gE = G*ME/RE^2;

%Materialien [Eisen, Kevlar, Dyneema, Karbonfasern]
sigma  = [100,4,3,5]*10^9;          % Spannnungsfestigkeit in N/m^2=Pa
dichte = [1300,1440,900,7900];      % Dichte in kg / m^3
material(1) = "Karbon";
material(2) = "Kevlar";
material(3) = "Dyneema";
material(4) = "Stahl";
titel(1) = "Erde";
titel(2) = "Mars";

%%
% Berechnungen Erde index = 1 Mars index = 2
figure();

for index=1:2
    if index == 2
        ME = MM;
        RE = RM;
        omegaE = omegaM;          
        gE = G*ME/RE^2;
    end
    fprintf('\n ----------------------------------------------------- \n');

    if index ==1
        fprintf('\n Berechnungen Erde \n');
    else
        fprintf('\n Berechnungen Mars \n'); 
    end        
    %Planetoationärer Radius
    RG     = (G*ME/omegaE^2)^(1/3);
    RGinkm = RG/1000;
    fprintf('\n Planetostationärer Radius');
    fprintf('\n %8.0f km ', RGinkm);

    %Seillänge
    LSeil = RE/2*(sqrt(1+8*RG^3/RE^3)-3);
    Rmax  = LSeil + RE;
    LSeilinkm = LSeil/1000;
    fprintf('\n Seillänge');
    fprintf('\n %8.0f km', LSeilinkm);
    fprintf('\n Abstand vom Zentrum');
    fprintf('\n %8.0f km', Rmax/1000);
    fprintf('\n');

    %Maximale Seilspannung in N/m^2
    sigma_max = G* ME*dichte*(1/RE-3/2/RG+RE^2/2/RG^3);
    %in N/mm^2
    sigma_max_GPa = sigma_max/10^9;
    fprintf('\n Maximal auftretende Seilspannung');
    for pindex= 1:4 
        fprintf('\n %s: \t  %4.0f GPa', material(pindex), sigma_max_GPa(pindex));
    end
    fprintf('\n');

    %Maximaler Durchmesser in m, 5 cm Dicke am Boden
    du = 0.005;
    Au = pi*du^2/4;
    A_max = Au*exp(dichte*gE*RE^2./sigma*(1/RE+RE^2/2/RG^3-3/2/RG));
    d_max = 2*sqrt(A_max/pi);
    %in m
    fprintf('\n Taper Ratio und Max Durchmesser @ RG bei Seildicke %5.3f cm am Seilanfang/ende', du*100);
    for pindex= 1:4 
        fprintf('\n %s: \t %8.2e  \t %2.1e cm', material(pindex), (d_max(pindex)/du)^2, d_max(pindex)*100);
    end
    fprintf('\n');
    r = linspace(RE,Rmax,1000);
    for k=1:4
        A(k,:) = Au*exp((dichte(k)*gE*RE^2/sigma(k))*(1/RE+RE^2/2/RG^3-1./r-r.^2/2/RG^3));
    end


    %Ab hier nur bei Erde noch fuer Karbonfasern, 
    %da alles andere keinen Sinn macht!!
    %Bei Mars für alle Materialien!
       
    %Berechnung Nutzlastgewinn 
    W1 = -G*ME*(1/RG-1/RE);
    W2 = omegaE^2*(RG)^2/2 + W1;
    Gewinn = (W1-W2)/W2;
    fprintf('\n Energiegewinn bei Nutzlasttransport auf RG');
    fprintf('\n %8.2f Prozent ', abs(Gewinn*100));

    %Berechnung Seildurchmesser über Abstand vom Erdmittelpunkt
    subplot(2,2,(index-1)+1);
    if index==1
        semilogy(r/10^6, 200*sqrt(A(1,:)/pi),'LineWidth',2);
    else
        for plotindex =1:3 
            semilogy(r/10^6, 200*sqrt(A(plotindex,:)/pi),'LineWidth',2);
            hold on;
        end
    end
    hold on;
    xlabel('Abstand r in 1000 km','FontSize',12);
    ylabel('Dicke Seil in cm','FontSize',12);
    ymax = 0.008 + (index-1)*0.024;
    ymin = 0.004;
    ylim([100*ymin, 100*ymax]);
    xlim([0,1.1*Rmax/10^6]);
    x=linspace(0,100,101);
    rx=0*x + RG;
    yx=(ymax-ymin)*x/length(x)+ymin;
    plot(rx/10^6,100*yx,'Color','black','LineStyle', '--','LineWidth',1); 
    if index==1 
        legend(material(1),'location','southeast');
    else
        legend(material(1),material(2),material(3),'location','east');
    end
    legend box off
    grid on;
    title(titel(index));
    set(gca,'Fontsize', 16);


    %Berechnung Gegenmasse M_c im Abstand DeltaR vom geostationären Orbit
    DeltaR = linspace(RE,Rmax-RG,1000);
    
%     Mc = Au*sigma(1)./(RG+DeltaR)./((omegaE^2-(G*ME)./(RG+DeltaR).^3))/1000;
    
    MS1 = RE^2*dichte(1)*gE/(2*sigma(1)*RG^3)*((2*RG^3+RE^3)/RE...
          -(2*RG^3+(RG+DeltaR).^3)./(RG+DeltaR));
    MS2= Au*sigma(1)/gE.*exp(MS1);
    Mc = MS2./(((RE^2*(RG+DeltaR))./RG.^3).*(1-(RG./(RG+DeltaR)).^3))/1000;

    %Berechnung Masse Seil
    fun = @(x) exp((dichte(1)*gE*RE^2/sigma(1))...
          .*(1/RE+RE^2/2/RG^3-1./x-x.^2/2/RG^3));
    for k=1:1000
       MasseSeil(k) = dichte(1)*Au*integral(fun,RE,RG+DeltaR(k))/1000;
%         MasseSeil(k) = dichte(1)*Au*(RG+DeltaR(k)-RE)/1000;
    end
    subplot(2,2,(index-1)+3);
    semilogy((RG+DeltaR)/10^6, MasseSeil,'LineWidth',2);
    hold on;
    semilogy((RG+DeltaR)/10^6, Mc,'LineWidth',2,'Color',Colors(2,:));
    semilogy((RG+DeltaR)/10^6, Mc+MasseSeil,'LineWidth',2,'Color',Colors(8,:));
    xlim([0, 1.1*Rmax/10^6]);
    ymin=MasseSeil(1);
    ymax=Mc(1);
    ymin=5*10^2;
    ymax=5*10^4;
    ylim([ymin, ymax]);
    x =linspace(1,100,100);
    rx=0*x + RG;
    yx=(ymax-ymin)*(x-1)/length(x)+ymin;
    plot(rx/10^6,yx,'Color','black','LineStyle', '--','LineWidth',1); 
    xlabel('Abstand r in 1000 km','FontSize',12);
    ylabel('Masse in 1000 kg','FontSize',12);
    legend(strcat(material(1),": Seilmasse"),strcat("Gegengewicht M_c"),...
    "Gesamtmasse",'location','northeast');
    legend box off
    grid on
    set(gca,'Fontsize', 16);
end
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------



