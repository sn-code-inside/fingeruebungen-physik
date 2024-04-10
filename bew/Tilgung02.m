% -------------------------------------------------------------------------
% Tilgung02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Schwingungsdämpfung an Hochhäusern
% 
% Programm berechnet Frequenzantwort Schwingungsdämpfung an Hochhäusern
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];


%% Numerik
% Parameter

gamma   =   0.05;          % Masseverhältnis
mT      =   01e3           % Masse Tilger in kg
M       =   mT/gamma;      % Masse Gerät in kg
kM      =   (2*pi*10)^2*M; % Federkonstante Gerät
kT      =   kM*gamma;      % Federkonstante Gegenmasse  
                           % (abgestimmt auf Eigenfrequenz)
etaM    =  0  ;            % Dämpfungskoeffizient
etaT    =  1e4;            % Dämpfungskoeffizient

Ax      = 0.30;            % Anfangsamplitude externe Anregung;
omegaM  = sqrt(kM/M);
omegaT  = sqrt(kT/mT);

fprintf('\n ');
fprintf('\n ');
fprintf('\n Masseverhältnis :            gamma= %8.2f ', gamma);
fprintf('\n eff. Masse Gebäude :         M    = %8.2f kg', M);
fprintf('\n Masse Tilger :               mT   = %8.2f kg', mT);
fprintf('\n Federkonstante Gebäude :     kM   = %8.2f kg/s² (N/m)', kM);
fprintf('\n Federkonstante Tilger :      kT   = %8.2f kg/s² (N/m)', kT);
fprintf('\n Dämpfung Gebäude :           etaM = %8.2f kg/s  (N/(m/s)', etaM);
fprintf('\n Dämpfung Tilger :            etaT = %8.2f kg/s  (N/(m/s)', etaT);
fprintf('\n Eigenfrequenz Gebäude :      fM   = %8.2f 1/s  ', omegaM/2/pi);
fprintf('\n Eigenfrequenz Tilger :       fT   = %8.2f 1/s  ', omegaT/2/pi);
fprintf('\n Anregungsamplitude  :        Ax   = %8.2f m ', Ax);
fprintf('\n ');





%% Berechnungen Frequenzantwort 
%  Variation etaT, Massenverhältnis gamma=0.05

omega = linspace(0.5*omegaM,1.5*omegaM,1000);
etaTv = [0, etaT, 1e20];
i = sqrt(-1);
for k=1:length(etaTv)
    A(:,k) = -omega(:).^2 + i*omega(:)*etaTv(k)/M + omegaM^2 + gamma*omegaT^2;
    B(:,k) = - i*omega(:)*etaTv(k)/M  - gamma*omegaT^2;
    C(:,k) = - i*omega(:)*etaTv(k)/mT - omegaT^2;
    D(:,k) = -omega(:).^2 + i*omega(:)*etaTv(k)/mT + omegaT^2;
    
    NN(:,k)  = A(:,k).*D(:,k)-B(:,k).*C(:,k);
    R(:,k)   = (-omega(:).^2 + i*omega(:)*etaTv(k)/mT + omegaT^2)./NN(:,k);
    R(:,k)   = R(:,k)/Ax;
    parastr1(k,:)  = string(strcat('\eta_T = ',num2str(etaTv(k),3)));
end

% Frequenzantwort 
figure();
for k= 1:length(etaTv)
     p1(k) = semilogy(omega/omegaM,abs(R(:,k)),'Color',Colors(k,:),...
            'LineWidth',2);
     hold on
end
grid on
ymin = 10*min(abs(R(:,1)));
ymax = 10*max(abs(R(:,2)));
axis ([omega(1)/omegaM omega(end)/omegaM, ymin ymax  ]);
ylabel('R(\omega)','FontSize',14)
xlabel('\omega/\omega_M','FontSize',14)
h=title('Schwingungstilgung ');
set(h,'FontSize',12,'FontWeight','normal'); 
legend(p1,parastr1, 'location','southeast',...
        'NumColumns',1);
legend box off
set(gca,'FontSize',16);


%% Berechnungen Frequenzantwort Anpassung Den-Hartog
% Optimum kT, etaT, Massenverhältnis gamma=0.02

omega = linspace(0.7*omegaM,1.3*omegaM,1000);
gamma = 0.05;
etaM  = 0;

fopt      = sqrt(1-0.5*gamma)/(1+gamma)
omegaTopt = omegaM*fopt
kTopt     = omegaTopt^2*mT;
xiopt     = sqrt(gamma*(3-sqrt(0.5*gamma))/(8*(1+gamma)*(1-0.5*gamma)))
etaTopt   = 2*xiopt*omegaTopt*mT;

fprintf('\n ');
fprintf('\n ');
fprintf('\n Masseverhältnis :            gamma= %8.2f ', gamma);
fprintf('\n eff. Masse Gebäude :         M    = %8.2f kg', M);
fprintf('\n Masse Tilger :               mT   = %8.2f kg', mT);
fprintf('\n Federkonstante Gebäude :     kM   = %8.2f kg/s² (N/m)', kM);
fprintf('\n Federkonstante Tilger :      kT   = %8.2f kg/s² (N/m)', kT);
fprintf('\n Federkonstante opt Tilger :  kTopt= %8.2f kg/s² (N/m)', kTopt);
fprintf('\n Dämpfung Gebäude :           etaM = %8.2f kg/s  (N/(m/s)', etaM);
fprintf('\n Dämpfung Tilger :            eta  = %8.2f kg/s  (N/(m/s)', etaT);
fprintf('\n Dämpfung Tilger opt :     etaTopt = %8.2f kg/s  (N/(m/s)', etaTopt);
fprintf('\n Eigenfrequenz Gebäude :      fM   = %8.2f 1/s  ', omegaM/2/pi);
fprintf('\n Eigenfrequenz Tilger :       fT   = %8.2f 1/s  ', omegaT/2/pi);
fprintf('\n Eigenfrequenz opt. Tilger :  fTopt= %8.2f 1/s  ', omegaTopt/2/pi);
fprintf('\n Anregungsamplitude  :        Ax   = %8.2f m ', Ax);
fprintf('\n ');

etaTv    = [0, etaT, 1e20, etaTopt];
omegaTv  = [omegaT, omegaT, omegaT, omegaTopt];
kTv      = [kTopt, kTopt, kTopt, kTopt];
for k=1:length(etaTv)
    A(:,k) = -omega(:).^2 + i*omega(:)*etaTv(k)/M + omegaM^2 + gamma*omegaTv(k)^2;
    B(:,k) = - i*omega(:)*etaTv(k)/M  - gamma*omegaTv(k)^2;
    C(:,k) = - i*omega(:)*etaTv(k)/mT - omegaTv(k)^2;
    D(:,k) = -omega(:).^2 + i*omega(:)*etaTv(k)/mT + omegaTv(k)^2;
    
    NN(:,k)  = A(:,k).*D(:,k)-B(:,k).*C(:,k);
    R(:,k)   = (-omega(:).^2 + i*omega(:)*etaTv(k)/mT + omegaTv(k)^2)./NN(:,k);
    R(:,k)   = R(:,k)/Ax;
    parastr2(k,:)  = string(strcat('\eta_T = ',num2str(etaTv(k),'%6.2e')));
end


% Frequenzantwort 
figure();
for k= 1:length(etaTv)
     p2(k) = semilogy(omega/omegaM,abs(R(:,k)),'Color',Colors(k,:),...
            'LineWidth',2);
     hold on
end
grid on
axis ([omega(1)/omegaM omega(end)/omegaM, ymin ymax]);
ylabel('R(\omega)','FontSize',14)
xlabel('\omega/\omega_M','FontSize',14)
h=title('Schwingungstilgung ');
set(h,'FontSize',12,'FontWeight','normal'); 
legend(p2,parastr2, 'location','southeast',...
        'NumColumns',2);
legend box off
set(gca,'FontSize',16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

