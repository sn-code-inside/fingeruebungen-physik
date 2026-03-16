% -------------------------------------------------------------------------
% Raumfahrtschrott.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Abschätzung der Verweildauer von Weltraummüll (exponentielle Atmosphäre)
%
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('IncludeFolder')
addpath('Data')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];


% Konstanten
Re   = 6371e3;                 % Erdradius [m]
mu   = 3.986004418e14;         % GM [m^3/s^2]

% Atmosphärenparameter (grob!)
rho_ref = 4e-12;               % Dichte bei href [kg/m^3] (z.B. ~400 km)
href    = 400e3;               % Referenzhöhe [m]
H       = 60e3;                % Skalenhöhe [m]

% Objektparameter (typischer Satellit / Trümmer)
CD      = 2.2;
A_over_m = 0.01;               % A/m [m^2/kg], z.B. 10 m^2 / 1000 kg

% Wiedereintrittsschwelle
h_reentry = 150e3;             % [m]

% Höhensampling für Lebensdauer-Kurve
h0_vec = linspace(400e3, 2000e3,50);  % Start-Höhen [m]
Tyears  = zeros(size(h0_vec));

wb1 = waitbar(0,'Bitte warten, Rechnung läuft');
for k = 1:numel(h0_vec)
    waitbar(k/length(h0_vec),wb1);
    h0 = h0_vec(k);
    a0 = Re + h0;
    Tyears(k) = lifetime_from_h(a0, h_reentry, mu, Re, ...
                                rho_ref, href, H, CD, A_over_m) ...
                / (365.25*24*3600); % in Jahre
end    
close(wb1)
figure;
plot(h0_vec/1e3, Tyears, 'LineWidth', 1.5);
set(gca,'YScale','log');
grid on;
xlabel('Bahnhöhe h_0 [km]');
ylabel('Lebensdauer [Jahre] (sehr grobes Modell)');
ht=title('Grobe Abschätzung der Verweildauer von Weltraummüll');
set(ht,'FontSize',12, 'FontWeight', 'normal');
set(gca,'FontSize',12, 'FontWeight', 'normal');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T = lifetime_from_h(a0, h_reentry, mu, Re, ...
                             rho_ref, href, H, CD, A_over_m)
    % Numerische Integration von a(t) bis h <= h_reentry
    a       = a0;
    t       = 0;
    dt_max  = 24*3600;       % max. Zeitschritt: 1 Tag
    dt_min  = 60;            % min. Zeitschritt: 1 Minute
    
    while true
        h = a - Re;
        if h <= h_reentry
            break;
        end
        
        % Atmosphärendichte (Exponentialsmodell)
        rho = rho_ref * exp(-(h - href)/H);
        
        % Bahngeschwindigkeit
        v = sqrt(mu / a);
        
        % da/dt (siehe Herleitung im Lösungsvorschlag)
        adot = - sqrt(mu) * CD * A_over_m * rho * sqrt(a);
        
        % adaptiver Zeitschritt (sehr grob)
        % z.B. so wählen, dass Änderung von a im Schritt klein bleibt
        dt = max(dt_min, min(dt_max, 0.001 * a / abs(adot)));
        
        % Euler-Schritt
        a = a + adot * dt;
        t = t + dt;
        
        % Sicherheitsabbruch, falls etwas schiefgeht
        if t > 1e12  % ~300 Jahre
            break;
        end
    end
    
    T = t;  % [s]
end
