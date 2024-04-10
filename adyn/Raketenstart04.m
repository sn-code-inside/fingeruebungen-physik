% -------------------------------------------------------------------------
% Raketenstart04.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Physikalische FingerÃ¼bungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
%
% -------------------------------------------------------------------------
%
% Programm berechnet Nutzlast-Antriebsvermögen-Relation für gleiche 
% Struktur- und Treibstoffmassen
%
% -------------------------------------------------------------------------

%%
clc
clear 
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors=GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

%%
% Initialisierung
% alle Daten in kg, m und s 

% Daten Ariane 1 
sigmaA  = [0.0655	0.0659	0.0872	0.0334];
vGA     = [2400	2585	3998	2708];       % m/s
m0A     = [202510	49240	12969	3600];   % kg
m0i     = m0A;
mLA0     = 2793;
mTi = [140000	33028	8238	685];        % kg
% Nutzlastbereich 
mLA      = linspace(0,100000,20001);

% Parameter für gleiche Struktur- und Treibstoffmassenverhältnisse
MStern = m0i./(m0i-mTi);
% Wir benutzen deshalb Werte für MStern zwischen 1.5 und 3.0
MStern = [1.5 2.0 2.5 3.0];
% vGN    = 3000; % m/s
vGN    = [2.5 3.0 3.5 4.0]*1000;
NEnd   = 10;


%%
% Berechnung Antriebsvermoegen und Nutzlast Ariane1 
Delta_vgesA = zeros(length(sigmaA),length(mLA));
for k=1:length(sigmaA)
   m0k           = m0A(k)  - mLA0 + mLA; %Korrektur 
   if k == length(sigmaA)
       m0kplus1   = mLA;
   else
       m0kplus1   = m0A(k+1)- mLA0 + mLA; %Korrektur
   end
   mu            = m0kplus1./m0k;
   Delta_v(k,:)  = vGA(k)*log(1./(sigmaA(k)+mu(:))); 
end
Delta_vgesA(1,:) =  Delta_v(1,:);
for k=2:length(sigmaA)
   Delta_vgesA(k,:) = Delta_vgesA(k-1,:) + Delta_v(k,:);   
end

% Berechnung Antriebsvermoegen und Nutzlast Ariane1 für gleiche 
% Struktur- und Treibstoffmassenverhältnisse (theoretische Betrachtung)
for k= 1:length(MStern)
    for N= 1:NEnd
        x(N)         = N;
        muN(k,N)     = MStern(3)^(-N);
        DeltavN(k,N) = vGN(k)*N*log(MStern(3));
    end
end

%%
% Graphische Ausgabe

figure('Name','Ariane1 und Vergleich Theorie');
subplot(1,2,1)  % Plot Nutzlastverhältnis über Antriebsvermögen
for k=1:length(sigmaA)
    semilogy(Delta_vgesA(k,:)/1000, mLA/m0A(1) ,'Color', Colors(k,:),'LineWidth',2);
    hold on
    lgdstrA(k,:) = sprintf('Stufe %s',num2str(k,2));
end
grid on;
xlim([0, 25]);
ylim([1e-3, 1]);
ylabel('Nutzlastverhältnis','FontSize',14)
xlabel('Antriebsvermögen Delta v in km/s','FontSize',14)
h2 = legend(lgdstrA,'location','northeast'); 
set(h2,'FontSize',14,'FontWeight','normal'); 
legend box off;
h3 =title('Ariane 1 ');
set(h3,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);

subplot(1,2,2)  % Plot Nutzlastverhältnis über Antriebsvermögen
for k=1:length(vGN)
    semilogy(DeltavN(k,:)/1000, muN(k,:), 'Color', Colors(k,:),'LineWidth',2);
    hold on
    lgdstrV(k,:) = sprintf('v_G %s km/s',num2str(vGN(k)/1000,'%3.1f'));
end
grid on;
xlim([0, 25]);
ylim([1e-3, 1]);
ylabel('Nutzlastverhältnis','FontSize',14)
xlabel('Antriebsvermögen Delta v in km/s','FontSize',14)
text(2,2e-3,strcat('M* = ',num2str(MStern(3),3)),'FontSize',14)
h2 = legend(lgdstrV,'location','northeast'); 
set(h2,'FontSize',14,'FontWeight','normal'); 
legend box off;
h3 =title('Gleiche Massenverhältnisse ');
set(h3,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);

for k= 1:length(MStern)
     for N= 1:NEnd
        x(N)         = N;
        muN(k,N)     = MStern(k)^(-N);
        DeltavN(k,N) = vGN(3)*N*log(MStern(k));
    end
end


% Berechnung Antriebsvermoegen und Nutzlast Ariane1 für gleiche 
% Struktur- und Treibstoffmassenverhältnisse (theoretische Betrachtung)
figure('Name','Gleiche Struktur- und Treibstoffmassenverhältnisse');
subplot(1,2,1)  % Plot Nutzlastverhältnis über Stufenzahl
for k=1:length(MStern)
    semilogy(x(:), muN(k,:),'Color', Colors(k,:),'LineWidth',2);
    hold on
    lgdstrV(k,:) = sprintf(' M* = %s   ', num2str(MStern(k),'%5.1f'));
end
grid on;
xlim([0, NEnd]);
ylim([1e-3, 1]);
ylabel('Nutzlastverhältnis','FontSize',14)
ylabel('Nutzlastverhältnis','FontSize',14)
xlabel('Stufenzahl','FontSize',14)
h2 = legend(lgdstrV,'location','northeast'); 
set(h2,'FontSize',14,'FontWeight','normal'); 
legend box off;
text(1,1.5e-3,strcat('v_G = ',num2str(vGN(3)/1000,'%3.1f km/s')),'FontSize',14)
h3 =title('Nutzlastverhältnis');
set(h3,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);

subplot(1,2,2)  % Plot Delta v über Stufenzahl
for k=1:length(MStern)
    plot(x(:), DeltavN(k,:)/1000,'Color', Colors(k,:),'LineWidth',2);
    hold on
end
grid on;
xlim([0, NEnd]);
ylabel('Antriebsvermögen Delta v in km/s','FontSize',14)
xlabel('Stufenzahl','FontSize',14)
h2 = legend(lgdstrV,'location','northwest'); 
set(h2,'FontSize',14,'FontWeight','normal'); 
legend box off;
text(3,3,strcat('v_G = ',num2str(vGN(3)/1000,'%3.1f km/s')),'FontSize',14)
h3 =title('Antriebsvermoegen');
set(h3,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);

