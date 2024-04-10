%%-----------------------------------------------------------------------%%
%
% Raketenoptimierung.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
%
%--------------------------------------------------------------------------
%
% Programm optimiert Parameter einer Mehrstufenrakete mit der 
% Methode der Lagrange-Multiplikatoren
% 
% Beispiel Saturn V - Flug zum Mond
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
% Parameter
sigmaS  = [0.0765	0.114 0.111];           % Strukturkoeffizienten
vGS     = [2980	4130 4130];                 % Gasgeschwindiglkeiten in m/s
m0S     = [2902000	658000	165000];        % Untermassen in kg
mLS     = 47000;                            % Nutzlast in kg
muS     = m0S/m0S(1);                       % µ_i
muLS    = mLS/m0S(1);

% Optimierungsbereich Nutzlastverhältnis
muL      = linspace(0.0,0.06,120);         

vG1 =vGS(1); 
vG2 =vGS(2);
vG3 =vGS(3);
sigma1 = sigmaS(1);
sigma2 = sigmaS(2);
sigma3 = sigmaS(3);

%%
% Methode Lagrange-Multiplikatoren 

% Wir geben ein Geschwindigkeitsbudget vor und suchen für das gegebene 
% Antriebsvermögen die zugehörige Nutzlast m_L bzw. das Nutzlastverhältnis
% mu_L.
% Gegeben sind nur die Strukturkoeefizienten pro Stufe und die 
% Ausströmgeschwindigkeiten der Gse pro Stufe.
% Wir wollen vorgeben, dass wir mindestens 12.0 km/s brauchen, um zum Mond
% zu kommen, besser etwas mehr (Die NASA ging von 12.4 km/s aus).

vmin  = 12000;
vmax  = 14000;
DVNASA= 12400;
DV = linspace(vmin,vmax,201);
for k=1:length(DV)
    DVact = DV(k);
    myfun = @(gamma,vG1,vG2,vG3,sigma1,sigma2,sigma3, DVact) - DVact + ...
            vG1*log((vG1-gamma)/sigma1/vG1)+...
            vG2*log((vG2-gamma)/sigma2/vG2)+...
            vG3*log((vG3-gamma)/sigma3/vG3);
    fun = @(gamma) myfun(gamma,vG1,vG2,vG3,sigma1,sigma2,sigma3, DVact); 
    % Achtung Einschränkung des Bereiches, so dass Nullstellen reell
    % bleiben. Siehe Bemerkung im Text.
    Gamma(k) = fzero(fun,[2000,2760]);
end

for k=1:length(Gamma)
    prod(1,k)=Gamma(k)*sigmaS(1)./((1-sigmaS(1)).*(vGS(1)-Gamma(k)));
    for j=2:length(vGS)
        prod(j,k) = prod(j-1,k).*Gamma(k)*sigmaS(j)./...
            ((1-sigmaS(j)).*(vGS(j)-Gamma(k)));
    end
    muLL(k) = prod(length(vGS),k);
    if round(DV(k)) == DVNASA
        muMax =  muLL(k);
        muStr = strcat('max \mu_L =',num2str(muMax,'%5.4f'));
    end

end


%%

% Graphische Darstellung

figure('Name','Drei-Stufen Optimierung')
semilogy(DV/1000,muLL,'color',Colors(2,:),'Linewidth',2,...
                     'LineStyle', Style(1));
x = [DVNASA DVNASA]/1000;
y = [min(muLL) max(muLL)];
p1 = line(x,y,'color',Colors(2,:),'Linewidth',2,'LineStyle', Style(3));
x  = [vmin DVNASA]/1000;
y  = [0.0162 0.0162];
p2 = line(x,y,'color',Colors(4,:),'Linewidth',2,'LineStyle', Style(3));
x  = [vmin DVNASA]/1000;
y  = [muMax muMax];
p2 = line(x,y,'color',Colors(3,:),'Linewidth',2,'LineStyle', Style(3));
ylabel('max \mu_L','FontSize',14)
xlabel('\Delta\it v','FontSize',14)
grid on
legend('optimiert Saturn V','Budget 12.4 km/s',...
       '\mu_L=0.0162 (NASA)',muStr,...
       'location','northeast','numColumns',1)
legend box off
set(gca,'FontSize',16);



