% -------------------------------------------------------------------------
% Tilgung03.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Schwingungsdämpfung
% 
% Programm berechnet Frequenzantwort Schwingungsdämpfung für Taipei 101 
% und den Berliner Fernsehturm
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];



%% Numerik
% Parameter Taipei 101
g       =   9.81;           % g in m/s^2            
gamma   =   0.05;           % Masseverhältnis
mT      =   662e3;          % Masse Pendel in kg
LP      =   42;             % Laenge Pendel in mg
M       =   mT/gamma;       % effektive Masse Gebäude in kg
TG      =   12.5;           % Schwingungsperiode Gebäude
omegaM  =   2*pi/TG;  
kM      =   sqrt(omegaM)*M; % effektive Federkonstante Gebäude
omegaT  =   sqrt(g/LP);
etaM    =   1e3;            % Dämpfungskoeffizient Gebäude
etaT    =   0;              % Dämpfungskoeffizient Pendel/Tilger

fprintf('\n ');
fprintf('\n Taipei 101');
fprintf('\n ');
fprintf('\n Masseverhältnis :            gamma= %8.2f ', gamma);
fprintf('\n eff. Masse Gebäude :         M    = %8.2e kg', M);
fprintf('\n Masse Tilger :               mT   = %8.2e kg', mT);
fprintf('\n Dämpfung Gebäude :           etaM = %8.2f kg/s  (N/(m/s)', etaM);
fprintf('\n Eigenfrequenz Gebäude :      fM   = %8.2f 1/s  ', omegaM/2/pi);
fprintf('\n Eigenperiode Gebäude :       TP   = %8.2f 1/s  ', TG);
fprintf('\n Eigenfrequenz Tilger :       fT   = %8.2f 1/s  ', omegaT/2/pi);
fprintf('\n Eigenfrequenz Tilger :     omegaT = %8.2f 1/s  ', omegaT);

% Berechnungen Frequenzantwort Taipei 101
%  Variation etaT

omega = linspace(0.75*omegaT,1.25*omegaT,1000);
etaTv = [1e3, 5e3, 2e4, 1e2];

i = sqrt(-1);
for k=1:length(etaTv)-1
    A(:,k) = -omega(:).^2  + i*omega(:)*etaTv(k)*gamma/mT/(1+gamma) + omegaM^2;
    B(:,k) = -omega(:).^2/LP + i*omega(:)*etaTv(k)/mT/LP;
    C(:,k) = -gamma/(1+gamma)*LP*omega(:).^2 + i*omega(:)*etaTv(k)*LP*gamma/mT/(1+gamma);
    D(:,k) = -omega(:).^2 + i*omega(:)*etaTv(k)/mT + omegaT^2;
    
    NN(:,k)  = A(:,k).*D(:,k)-B(:,k).*C(:,k);
    R(:,k)   = D(:,k)./NN(:,k);
    parastr1(k,:)  = string(strcat('\eta_T = ',num2str(etaTv(k),3)));
end

% Frequenzantwort Zero
R(:,k+1)   = 1./(-omega(:).^2+omegaM^2+i*omega(:)*etaM/M);
parastr1(k+1,:)  = string(strcat('m_T = ','0'));

% Frequenzantwort 
figure();
for k= 1:length(etaTv)-1
     p1(k) = semilogy(omega/omegaT,abs(R(:,k)),'Color',Colors(k,:),...
            'LineWidth',2);
     hold on
end
p1(k+1) = semilogy(omega/omegaT,abs(R(:,k+1)),'Color',Colors(k+1,:),...
            'LineWidth',1,'LineStyle',Style(2));
fprintf('\n Dämpfung Tilger :            etaT = %8.2f kg/s  (N/(m/s)', etaTv(2));
fprintf('\n Anregungsamplitude m. Tilg.: R Ax = %8.2f m ',  min(abs(R(:,2))));
fprintf('\n Beschleunigung m. Tilg. :    aT   = %8.2f m/s^2 ', min(abs(R(:,2)))*omegaT^2);
fprintf('\n ');

grid on
ymin = min(abs(R(:,1)));
ymax = 10*max(abs(R(:,1)));
axis ([omega(1)/omegaT omega(end)/omegaT, ymin ymax]);
ylabel('R(\omega)','FontSize',14)
xlabel('\omega/\omega_T','FontSize',14)
h=title('Taipei 101 ');
set(h,'FontSize',12,'FontWeight','normal'); 
legend(p1,parastr1, 'location','southeast',...
        'NumColumns',1);
legend box off
set(gca,'FontSize',16);

%%
% Parameter Berliner Fernsehturm

g       =   9.81;           % g in m/s^2            
gamma   =   0.05;           % Masseverhältnis
mT      =   1.5e3;          % Masse Pendel in kg
LP      =   5.0;            % Laenge Pendel in mg
M       =   mT/gamma;       % effektive Masse Gebäude in kg
TG      =   5.5;            % Schwingungsperiode Gebäude
omegaM  =   2*pi/TG;  
kM      =   omegaM^2*M; % effektive Federkonstante Gebäude
omegaT  =   sqrt(g/LP);
etaM    =   01e2;           % Dämpfungskoeffizient Gebäude
etaT    =   0;              % Dämpfungskoeffizient Pendel/Tilger

fprintf('\n ');
fprintf('\n Berliner Fernsehturm');
fprintf('\n ');
fprintf('\n Masseverhältnis :            gamma= %8.2f ', gamma);
fprintf('\n eff. Masse Gebäude :         M    = %8.2e kg', M);
fprintf('\n Masse Tilger :               mT   = %8.2e kg', mT);
fprintf('\n Dämpfung Gebäude :           etaM = %8.2f kg/s  (N/(m/s)', etaM);
fprintf('\n Eigenfrequenz Gebäude :      fM   = %8.2f 1/s  ', omegaM/2/pi);
fprintf('\n Eigenperiode Gebäude :       TP   = %8.2f 1/s  ', TG);
fprintf('\n Eigenfrequenz Tilger :       fT   = %8.2f 1/s  ', omegaT/2/pi);
% fprintf('\n Eigenfrequenz Tilger :     omegaT = %8.2f 1/s  ', omegaT);

% Berechnungen Frequenzantwort Berliner Fernsehturm
%  Variation etaT

omega = linspace(0.75*omegaT,1.2*omegaT,1000);
etaTv = [10, 25, 150, 10];

i = sqrt(-1);
for k=1:length(etaTv)-1
    A(:,k) = -omega(:).^2  + i*omega(:)*etaTv(k)*gamma/mT/(1+gamma) + omegaM^2;
    B(:,k) = -omega(:).^2/LP + i*omega(:)*etaTv(k)/mT/LP;
    C(:,k) = -gamma/(1+gamma)*LP*omega(:).^2 + i*omega(:)*etaTv(k)*LP*gamma/mT/(1+gamma);
    D(:,k) = -omega(:).^2 + i*omega(:)*etaTv(k)/mT + omegaT^2;   
    NN(:,k)  = A(:,k).*D(:,k)-B(:,k).*C(:,k);
    R(:,k)   = D(:,k)./NN(:,k);
    parastr1(k,:)  = string(strcat('\eta_T = ',num2str(etaTv(k),3)));
end

R(:,k+1)   = 1./(-omega(:).^2+omegaM^2+i*omega(:)*etaM/M);
parastr1(k+1,:)  = string(strcat('m_T = ','0'));

% Frequenzantwort 
figure();
for k= 1:length(etaTv)-1
     p1(k) = semilogy(omega/omegaT,abs(R(:,k)),'Color',Colors(k,:),...
            'LineWidth',2);
     hold on
end
p1(k+1) = semilogy(omega/omegaT,abs(R(:,k+1)),'Color',Colors(k+1,:),...
            'LineWidth',1,'LineStyle',Style(2));
fprintf('\n Dämpfung Tilger :            etaT = %8.2f kg/s  (N/(m/s)', etaTv(2));
fprintf('\n Anregungsamplitude m. Tilg.: R Ax = %8.2f m ',min(abs(R(:,2))));
fprintf('\n Beschleunigung m. Tilg. :    aT   = %8.2f m/s^2 ', min(abs(R(:,2)))*omegaT^2);
fprintf('\n ');

grid on
ymin = min(abs(R(:,1)));
ymax = 10*max(abs(R(:,1)));
axis ([omega(1)/omegaT omega(end)/omegaT, ymin ymax]);
ylabel('R(\omega)','FontSize',14)
xlabel('\omega/\omega_T','FontSize',14)
h=title('Berliner Fernsehturm ');
set(h,'FontSize',12,'FontWeight','normal'); 
legend(p1,parastr1, 'location','southeast',...
        'NumColumns',1);
legend box off
set(gca,'FontSize',16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

