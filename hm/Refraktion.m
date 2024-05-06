% -------------------------------------------------------------------------
% Refraktion.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Refraktionsberechnung für das Kugelschalen-Schichtenmodell (3.202)
% und Vergleicht mit Normalrefraktion.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Konstanten der Erdatmosphaere 
K0 = 0.0003161;      % dimensionslos 
n0 = 1+K0;           % Brechzahl auf Erdoberfläche
% Konstanten der Erde 
R_Erde = 6366;  % [R_Erde] = km 

% Wir bestimmen uns zunächst H aus Datenwerten:
% z.B. aus https://github.com/sky-s/standard-atmosphere
% oder https://web.archive.org/web/20170516220116/
% oder http://home.anadolu.edu.tr/~mcavcar/common/ISAweb.pdf
rhoex=[10000 9711 9428 9151 8881 8617 8359 8106 7860 7620 7385 ...
       7156 6932 6713 6500 6292];
rhoex =rhoex/10000;
hn = linspace(0,15,16); % in 1000 feet
hn = hn/3.28084; % in 1000 meter

% Parameter [H] = km 
H0=-hn(length(hn))/log(rhoex(length(rhoex)));%Abgeleiteter Wert ISA Tabelle 
H = H0;
h =  linspace(0,50,100);  % scheinbare Zenitdistanz als Winkel im Gradmass 
rho = exp(-h/H);
nref=1+K0*rho;

%% 
% Dichte und Brechzahl von der Höhe

figure()
subplot(1,3,1);
semilogy(hn, rhoex, h, rho);
xlabel('Höhe in km');
ylabel('Dichte/Dichte_0');
legend('ISA Tabelle', 'Fit an '+ string(num2str(H,'%4.2f km')));
legend box off;
xlim([0 10]);
grid on;

subplot(1,3,2);
plot(h, nref);
xlabel('Höhe in km');
ylabel('Brechungsindex');
legend('Brechzahl mit {\kappa} = ' + string(num2str(K0,'%6.5f ')));
legend box off;
xlim([0 10]);
ylim([1 (1+1.1*K0)]);
grid on;

% Graphische Überprüfung des Integranden auf Pole und Nullstellen

NSt = 5000;                 %Stützstellen für Integration
NP  = 10;                   %Stützstellen für Zenitwinkel
z =  linspace(80,90,NP);    %scheinbare Zenitdistanz als Winkel im Gradmass 
xe = 2;
xw = linspace(0,xe,NSt);
for m=1:NP
    z1 = z(m); 
    C1(m)  = sind(z1);
    for n=1:NSt
      xi = xw(n); 
      RefIntegrand(m,n) = exp(-xi)./(1+K0.*exp(-xi))...
        ./sqrt((((1+K0.*exp(-xi))/(1+K0)).*(1+H.*xi/R_Erde)).^2-C1(m)^2);
    end
end
subplot(1,3,3);
for m=1:NP 
    semilogy(xw,RefIntegrand(m,:));
    hold on;
end
xlabel('x');
ylabel('Integrand');
legend(string(num2str(90-z(:),'Höhe %4.1f °')));
legend box off;
xlim([0 xe]);
grid on;

figure()
NP  = 100;                     %Stützstellen für Zenitwinkel
z =  linspace(0,10,NP);        %scheinbare Zenitdistanz im Gradmass 
R2=1.0./tand(z+7.31./(z+4.4)); %Refraktionsformel (3.77)
RefInt0 = zeros(NP);

for k=1:2
    if k==1 
        H = H0; % Werte aus ISA Tabellen
        K = K0;      
    else
        H = 8.59735;  %Benutzte Werte für H und K von Newton
        K = 0.000263462;
    end
    RefInt0(:) = K*sqrt(pi*R_Erde/2/H); %Näherungswert für h' = 0 
                                        %aus Formel (3.202)
    %Integrand
    fun = @(x,C,H,K,R_E) exp(-x)./(1+K.*exp(-x))./...
            sqrt((((1+K.*exp(-x))/(1+K)).*(1+H.*x/R_E) ).^2 - C^2 );
    %Integration Formel (3.202)
    for m=1:length(z)
        Csin = sind(90-z(m));
        RefInt(m) = Csin*K*integral(@(x)fun(x,Csin,H,K,R_Erde),0,100);
    end
    plot(z,60*rad2deg(RefInt),'Color',Colors(k+1,:),...
         'LineWidth',2,'LineStyle','-');
    hold on
end
plot(z,R2,'Color',Colors(4,:),'LineWidth',2,'LineStyle','-.');
plot(z,60*rad2deg(RefInt0),'Color',Colors(8,:),'LineWidth',2,...
     'LineStyle',':');
xlabel("Scheinbare Höhe \it h' \rm in °");
ylabel('Refraktion R in Bogenminuten');  
lgtitle1= ' Werte aus ISA Tabellen \it H \rm = ' + ...
        string(num2str(H0,'%5.2f'))+ ...
        ' {\kappa} \rm = ' + string(num2str(K0,'%8.6f'));
lgtitle2=  " Newton's Werte (1694) " + ' \it H \rm = ' + ...
        string(num2str(H,'%5.2f'))+...
        ' {\kappa} \rm = ' + string(num2str(K,'%8.6f'));
legend(lgtitle1, lgtitle2, " Formel (3.77)",...
        " Näherungswert für Formel (3.202) für \it h' \rm = 0 ",...
        'location','northeast');
legend box off;
xlim([0 10]);
ylim([0 40]);
grid on;
set(gca,'FontSize',16);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
