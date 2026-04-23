clc; clear; close all;

%% =========================
% DATOS DE ENTRADA (ORURO)
% =========================
LAT = -17.99;              % grados (Sur negativo)
LONG = -67.13;             % grados
n = 172;                   % día del año (ej: 21 junio)
gamma = deg2rad(25);       % inclinación panel
alpha = 0;                 % orientación (0 = norte)
RF = 0.17;                 % reflectividad suelo

% Datos medidos (ejemplo típico)
t = 15/60; % horas (intervalo 15 min)
IHN = rand(1,96)*800; % radiación instantánea simulada W/m2

%% =========================
% GEOMETRÍA SOLAR
% =========================

% Declinación (ecuación 1)
DEC = 45.23 * sind((360/365.25)*(n-81));

% Vector de horas solares
H = linspace(0,24,96);

% Ángulo horario (ecuación 2)
w = 15*(H - 12); % grados

% Altura solar (ecuación 5)
ALT = asind(cosd(LAT).*cosd(DEC).*cosd(w) + sind(LAT).*sind(DEC));

% Evitar valores negativos (noche)
ALT(ALT<0) = 0;

%% =========================
% RADIACIÓN DIARIA
% =========================

% (10) HHN
HHN = sum(IHN * t); % Wh/m2-día

% (11) Radiación extraterrestre HO
WHE = acosd(-tand(LAT)*tand(DEC));

Io = 1367; % W/m2 constante solar

HO = (24/pi)*Io*( ...
    cosd(LAT)*cosd(DEC)*sind(WHE) + ...
    (pi/180)*WHE*sind(LAT)*sind(DEC));

% (12) Índice de nubosidad
K = HHN / HO;

% (13) Fracción difusa (Liu & Jordan)
FD = 1.3903 - 4.0273*K + 5.5315*K^2 - 3.108*K^3;

% (14) Difusa
DFN = FD * HHN;

% (15) Directa
DRN = HHN - DFN;

%% =========================
% SUPERFICIE INCLINADA
% =========================

% Aproximación RB (simplificada física)
RB = cos(gamma); % aproximación (ver nota abajo)

% (17)
HDR = DRN * RB;

% (18)
HDF = DFN * (1 + cos(gamma))/2;

% (19)
HRF = DFN * RF * (1 - cos(gamma))/2;

% (20)
HIT = HDR + HDF + HRF;

%% =========================
% RADIACIÓN INSTANTÁNEA
% =========================

% (21) factor rd
WHE_rad = deg2rad(WHE);
w_rad = deg2rad(w);

rd = (pi/24) * (cos(w_rad) - cos(WHE_rad)) ./ ...
    (sin(WHE_rad) - WHE_rad*cos(WHE_rad));

% Evitar errores numéricos
rd(isnan(rd)) = 0;

% Difusa instantánea
IDFH = rd * DFN;

% (22) directa instantánea
IDRH = IHN - IDFH;

%% =========================
% SUPERFICIE INCLINADA INSTANTÁNEA
% =========================

% (23)
rb = cos(gamma); % simplificación

IDR = IDRH .* rb;

% (24)
IDF = IDFH .* (1 + cos(gamma))/2;

% (25)
IRF = IHN .* RF .* (1 - cos(gamma))/2;

% Total
IT = IDR + IDF + IRF;

%% =========================
% GRÁFICAS
% =========================

figure;
plot(H, IHN, 'k', 'LineWidth',1.5); hold on;
plot(H, IDRH, '--', 'LineWidth',1.5);
plot(H, IDFH, ':', 'LineWidth',1.5);
xlabel('Hora');
ylabel('W/m^2');
legend('Global','Directa','Difusa');
title('Radiación solar instantánea');

figure;
plot(H, IT, 'LineWidth',1.5);
xlabel('Hora');
ylabel('W/m^2');
title('Radiación en superficie inclinada');

%% =========================
% RESULTADOS
% =========================
fprintf('HHN = %.2f Wh/m2-dia\n', HHN);
fprintf('DRN = %.2f\n', DRN);
fprintf('DFN = %.2f\n', DFN);
fprintf('HIT = %.2f\n', HIT);
