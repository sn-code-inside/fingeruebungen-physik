% -------------------------------------------------------------------------
% Raketenstart03.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Programm berechnet Nutzlast-Antriebsvermögen-Relation für Ariane 1 und
% Saturn V Rakete
%
%--------------------------------------------------------------------------


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
vGA     = [2400	2585	3998	2708];
m0A     = [202510	49240	12969	3600];
mLA0     = 2793;
% Daten Saturn V 
sigmaS  = [0.05	0.06 0.06];
vGS     = [2761	4561 4549];
m0S     = [2820950	675950	217250];
mLS0    = 100000;

% Nutzlastbereich 
mLA      = linspace(0,100000,20001);
mLS      = linspace(0,1000000,20001);

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

Delta_vgesS = zeros(length(sigmaS),length(mLS));
for k=1:length(sigmaS)
   m0k           = m0S(k)  - mLS0 + mLS; %Korrektur 
   if k == length(sigmaS)
       m0kplus1  = mLS;
   else
       m0kplus1  = m0S(k+1)- mLS0 + mLS; %Korrektur
   end
   mu            = m0kplus1./m0k;
   Delta_v(k,:)  = vGS(k)*log(1./(sigmaS(k)+mu(:))); 
end
Delta_vgesS(1,:) =  Delta_v(1,:);
for k=2:length(sigmaA)
   Delta_vgesS(k,:) = Delta_vgesS(k-1,:) + Delta_v(k,:);   
end



%%
% Graphische Ausgabe

figure('Name','Nutzlastverhältnis');
subplot(1,2,1)  % Plot Nutzlastverhältnis über Antriebsvermögen
for k=1:length(sigmaA)
    p(k)=semilogy(Delta_vgesA(k,:)/1000, mLA/m0A(1) ,'Color', Colors(k,:),'LineWidth',2);
    hold on
    lgdstrA(k,:) = sprintf('Stufe %s',num2str(k,2));
end
line 
grid on;
xlim([0, 25]);
ylim([1e-3, 1]);
line([7.73 7.73],[0.001 1],'Color', Colors(8,:),...
    'LineStyle',Style(3),'LineWidth',1);
ylabel('Nutzlastverhältnis','FontSize',14)
xlabel('Antriebsvermögen Delta v in km/s','FontSize',14)
h2 = legend(p, lgdstrA,'location','northeast'); 
set(h2,'FontSize',14,'FontWeight','normal'); 
ht = title('Ariane 1'); 
set(ht,'FontSize',14,'FontWeight','normal'); 
legend box off;
set(gca,'FontSize',16);


subplot(1,2,2)  % Plot Nutzlastverhältnis über Antriebsvermögen
for k=1:length(sigmaS)
    p2(k)=semilogy(Delta_vgesS(k,:)/1000, mLS/m0S(1) ,'Color', Colors(k,:),...
    'LineStyle',Style(1),'LineWidth',2);
    hold on
    lgdstrS(k,:) = sprintf('Stufe %s',num2str(k,2));
end
grid on;
xlim([0, 25]);
ylim([1e-3, 1]);
line([7.73 7.73],[0.001 1],'Color', Colors(8,:),...
    'LineStyle',Style(3),'LineWidth',1);
ylabel('Nutzlastverhältnis','FontSize',14)
xlabel('Antriebsvermögen Delta v in km/s','FontSize',14)
h2 = legend(p2,lgdstrS,'location','northeast'); 
set(h2,'FontSize',14,'FontWeight','normal'); 
legend box off;
ht = title('Saturn V');
set(ht,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);
