% -------------------------------------------------------------------------
% LambertSolver2.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Physikalische FingerÃ¼bungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Routine zur Berechnung von a
% für das Lambert Problem durch ein Bisektionsverfahren mit 
% Angabe der einzelnen Iterationen.
%
% -------------------------------------------------------------------------
% Inputs:         Description                     Range/Unit
%   TOF         - TOF                                  d
%   r1,r2       - Abstände P1, P2                      AE
%   gamma       - Winkel zwische Vektoren r1 und r2    rad
%   mu          - Gravitationsparameter                AE^3/d^3
% 
% OutPuts:
%   aiter       - Vektor Iteration für a               AE
%   alpha_n     - Iterationswert alpha                 rad
%   beta_n      - Iterationswert beta                  rad
%   tau_n       - Iterationswert tau                   d
%   tm          - Zeit für MinimumEnergietransferbahn  d
%   tpar        - Zeit für Parabelbahn                 d
%   nend        - Zahl der Iterationsschritte (<99)         
%
% -------------------------------------------------------------------------

% Lambert-Problem Iterationsverfahren Nr.2 
function [aiter,amin,amax,alpha_n,beta_n,tau_n,tm,tpar,nend] = ...
          LambertSolver2(r1, r2, gamma, TOF, mu)
    c     = sqrt(r1^2+r2^2-2*r1*r2*cos(gamma));
    s     = (c+r1+r2)/2;
    if gamma < pi
        fac_beta  = 1;
    else
        fac_beta = -1;
    end
    amin0    = s/2;
    alpham  = wrapTo2Pi(2*asin(sqrt(s/2/amin0)));
    alphamd = rad2deg(alpham);
    betam   = fac_beta*2*asin(sqrt((s-c)/2/amin0));
    betamd  = rad2deg(betam);
    tm      = sqrt(amin0^3/mu)*(alpham-betam-(sin(alpham)-sin(betam)));
    tpar    = (sqrt(2)/3)*sqrt(s^3/mu)*(1-fac_beta*((s-c)/s)^(3/2));
    if TOF < tm
        fac_alpha = 1;
        sum_alpha = 0;
    else
        fac_alpha = -1;
        sum_alpha = 2*pi;
    end
    amax = zeros(1,100);
    amin = zeros(1,100);
    amax(1)    =  2*amin0;
    amin(1)    =  1*amin0;
    aiter   = zeros(1,100);
    beta_n  = zeros(1,100);
    alpha_n = zeros(1,100);
    tau_n   = zeros(1,100);
    aiter(1)= (amin(1)+amax(1))/2;
    n       = 1;
    nend    = 100;
    while abs(1-tau_n(n)/TOF) > 0.001 && n <100
        alpha_n(n)  = wrapTo2Pi(sum_alpha+fac_alpha*2*asin(sqrt(s/2/aiter(n))));
        beta_n(n)   = fac_beta*2*asin(sqrt((s-c)/2/aiter(n)));
        tau_n(n+1)  = sqrt(aiter(n)^3/mu)*...
                   (alpha_n(n)-beta_n(n)-(sin(alpha_n(n))-sin(beta_n(n))));
        if tau_n(n+1) - TOF  < 0 
            amax(n+1) = amax(n);
            amin(n+1) = aiter(n);
        else
            amax(n+1) = aiter(n);
            amin(n+1) = amin(n);
        end
        aiter(n+1)= (amin(n+1)+amax(n+1))/2;
        n = n+1;
    end
    nend  = n-1; 
    % Falls keine Konvergenz, umgekehrte Richtung
    if nend >98
        amax = zeros(1,100);
        amin = zeros(1,100);
        amax(1)    =  1*amin0;
        amin(1)    =  2*amin0;
        aiter   = zeros(1,100);
        beta_n  = zeros(1,100);
        alpha_n = zeros(1,100);
        tau_n   = zeros(1,100);
        aiter(1)= (amin(1)+amax(1))/2;
        n       = 1;
        nend    = 100;
        while abs(1-tau_n(n)/TOF) > 0.001 && n <100
            alpha_n(n)  = wrapTo2Pi(sum_alpha+fac_alpha*2*asin(sqrt(s/2/aiter(n))));
            beta_n(n)   = fac_beta*2*asin(sqrt((s-c)/2/aiter(n)));
            tau_n(n+1)  = sqrt(aiter(n)^3/mu)*...
                       (alpha_n(n)-beta_n(n)-(sin(alpha_n(n))-sin(beta_n(n))));
            if tau_n(n+1) - TOF  < 0 
                amax(n+1) = amax(n);
                amin(n+1) = aiter(n);
            else
                amax(n+1) = aiter(n);
                amin(n+1) = amin(n);
            end
            aiter(n+1)= (amin(n+1)+amax(n+1))/2;
            n = n+1;
        end
        nend  = n-1;
    end
end