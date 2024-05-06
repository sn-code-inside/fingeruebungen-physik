% -------------------------------------------------------------------------
% Tilgung01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
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
mT      =   10000;         % Masse Tilger in kg
M       =   mT/gamma;      % Masse Gerät in kg
kM      =  75e7;           % Federkonstante Gerät
kT      =  kM*gamma;       % Federkonstante Gegenmasse  
                           % (abgestimmt auf Eigenfrequenz)
etaM    =  3000;           % Dämpfungskoeffizient
etaT    =  3000;           % Dämpfungskoeffizient
tmax    =   3;

Ax      = 0.30;            % Anfangsamplitude externe Anregung;
omegaM  = sqrt(kM/M);
omegaT  = sqrt(kT/mT);
Omega   = omegaT;

fprintf('\n ');
fprintf('\n ');
fprintf('\n eff. Masse Gebäude :         M    = %8.2f kg', M);
fprintf('\n Masse Tilger :               mT   = %8.2f kg', mT);
fprintf('\n Federkonstante Gebäude :     kM   = %8.2f kg/s² (N/m)', kM);
fprintf('\n Federkonstante Tilger :      kT   = %8.2f kg/s² (N/m)', kT);
fprintf('\n Dämpfung Gebäude :           etaM = %8.2f kg/s  (N/(m/s)', etaM);
fprintf('\n Dämpfung Tilger :            etaT = %8.2f kg/s  (N/(m/s)', etaT);
fprintf('\n Eigenfrequenz Gebäude :      fM   = %8.2f 1/s  ', omegaM/2/pi);
fprintf('\n Eigenfrequenz Tilger :       fT   = %8.2f 1/s  ', omegaT/2/pi);
fprintf('\n Anregungsfrequenz  :         fext = %8.2f 1/s  ', Omega/2/pi);
fprintf('\n Anregungsamplitude  :        Ax   = %8.2f m ', Ax);
fprintf('\n ');


%% Berechnungen Frequenzantwort Variation Massenverhältnis 
%  etaM = 3000
%  etaT = 0;

omega = linspace(0.8*omegaM,1.25*omegaM,1000);
i   = sqrt(-1);
mTv = [0.2*mT, 0.4*mT 0.8*mT, 1*mT, 1.2*mT];
kTv = kM*mTv/M;

for k=1:length(mTv)
    NN(:,k)  = (-M*omega(:).^2 + i*omega(:)*etaM + kTv(k) + kM).*...
               (-mTv(k)*omega(:).^2 + kTv(k)) -...
               (kTv(k)).^2;
    R(:,k)   = kM*(-mTv(k)*omega(:).^2 + kTv(k))./NN(:,k);
    parastr(k,:)  = string(strcat('\gamma = ',num2str(mTv(k)/M,'%5.2f')));
end


% Frequenzantwort 
figure();
for k= 1:length(mTv)
     p(k) = semilogy(omega/omegaM,abs(R(:,k)),'Color',Colors(k,:),...
            'LineWidth',2);
     hold on
end
grid on
ymin = 1e-3;
ymax = 1.5*max(abs(R(:,5)));

axis ([omega(1)/omegaM omega(end)/omegaM,ymin,ymax]);
ylabel('R(\omega)','FontSize',14)
xlabel('\omega in 1/s','FontSize',14)
h=title('Schwingungstilgung ');
set(h,'FontSize',12,'FontWeight','normal'); 
legend(p,parastr, 'location','southeast',...
        'NumColumns',1);
legend box off
set(gca,'FontSize',16);

%% Berechnungen Frequenzantwort Variation etaT 
% Massenverhältnis gamma=0.05
% etaM = 0

omega = linspace(0.8*omegaM,1.25*omegaM,1000);
etaTv = [0.1*etaT, 1*etaT, 10*etaT];
etaM  = 0;

for k=1:length(etaTv)
    NN(:,k)  = (-M*omega(:).^2 + i*omega(:)*etaM + kT + kM).*...
               (-mT*omega(:).^2 + i*omega(:)*etaTv(k)+kT) -...
               (-i*omega(:)*etaTv(k) + kT).^2;
    R(:,k)   = kM*(-mT*omega(:).^2 + i*omega(:)*etaTv(k) + kT)./NN(:,k);
    parastr2(k,:)  = string(strcat('\eta_T = ',num2str(etaTv(k),'%6.1e')));
end


% Frequenzantwort 
figure();
for k= 1:length(etaTv)
     p2(k) = semilogy(omega/omegaM,abs(R(:,k)),'Color',Colors(k,:),...
            'LineWidth',2);
     hold on
end
grid on
axis ([omega(1)/omegaM omega(end)/omegaM,ymin,ymax]);
ylabel('R(\omega)','FontSize',14)
xlabel('\omega/\omega_M','FontSize',14)
h=title('Schwingungstilgung ');
set(h,'FontSize',12,'FontWeight','normal'); 
legend(p2,parastr2, 'location','southeast',...
        'NumColumns',1);
legend box off
set(gca,'FontSize',16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

