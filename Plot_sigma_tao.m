clear; clc; close all;

%% 1. DATI GEOMETRICI E DI CARICO
F = 1500;           % Forza su ciascun corrente inferiore [N]
F_tot = 2 * F;      % Taglio verticale totale all'estremità libera [N]
a = 0.3;            % Lato della sezione quadrata [m]
A = 0.2;            % Area singolo corrente [m^2]
L = 1.5;            % Lunghezza del cassone [m]

I_xx = 4 * A * (a/2)^2;   % Momento d'inerzia della sezione idealizzata [m^4]

%% 2. SOLUZIONE ANALITICA IDEALE (SFORZO NORMALE N)
z_id = linspace(0, L, 200);
M_x = F_tot * (L - z_id); % Momento flettente lungo z [N*m]

% Sforzo normale assiale analitico N = sigma * A [N]
% Con carico rivolto verso l'alto:
% - Corrente superiore C3 compresso (segno negativo)
% - Corrente inferiore C4 teso (segno positivo)
N_C3_id = -(M_x / I_xx) * (a/2) * A;   % Da -7500 N (incastro) a 0 N (estremità)
N_C4_id = +(M_x / I_xx) * (a/2) * A;   % Da +7500 N (incastro) a 0 N (estremità)

%% 3. LETTURA DEI FILE .RPT DA PATRAN
[z_C3, N_C3_fem] = leggi_rpt_patran('C3_sforzo_normale.rpt');
[z_C4, N_C4_fem] = leggi_rpt_patran('C4_sforzo_normale.rpt');

%% 4. GRAFICO 1: SFORZO NORMALE ASSIALE SUI CORRENTI (SFONDO BIANCO)
fig1 = figure('Name', 'Sforzo Normale Correnti C3 e C4', 'Color', 'w', 'Position', [100 100 850 520]);
if isprop(fig1, 'Theme'), fig1.Theme = 'light'; end % Disattiva il tema scuro per questa figura

ax1 = axes(fig1, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.7 0.7 0.7]);
hold(ax1, 'on'); grid(ax1, 'on'); box(ax1, 'on');

% Curve analitiche (linee tratteggiate)
plot(ax1, z_id, N_C3_id, '--b', 'LineWidth', 2, 'DisplayName', 'Analitico - Corrente Superiore (C3)');
plot(ax1, z_id, N_C4_id, '--r', 'LineWidth', 2, 'DisplayName', 'Analitico - Corrente Inferiore (C4)');

% Curve FEM (linee continue con marcatori)
plot(ax1, z_C3, N_C3_fem, '-ob', 'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', 'FEM - Corrente Superiore (C3)');
plot(ax1, z_C4, N_C4_fem, '-sr', 'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', 'FEM - Corrente Inferiore (C4)');

xlabel(ax1, 'Coordinata longitudinale z [m]', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
ylabel(ax1, 'Sforzo normale N [N]', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
title(ax1, 'Confronto Sforzo Normale sui Correnti (C3 e C4): FEM vs Analitico', 'FontSize', 13, 'Color', 'k');
lgd1 = legend(ax1, 'Location', 'best', 'FontSize', 10);
set(lgd1, 'Color', 'w', 'TextColor', 'k', 'EdgeColor', 'k');
xlim(ax1, [0 L]);

% Esporta con sfondo forzatamente bianco
exportgraphics(fig1, 'plot_sforzo_normale_correnti.png', 'Resolution', 300, 'BackgroundColor', 'white');

%% FUNZIONE LOCALE DI LETTURA FILE .RPT
function [z_sorted, val_sorted] = leggi_rpt_patran(nome_file)
fid = fopen(nome_file, 'r');
if fid == -1
    error('File "%s" non trovato nella cartella corrente di MATLAB (%s).', nome_file, pwd);
end

dati = [];
while ~feof(fid)
    riga = strtrim(fgetl(fid));
    % Legge solo le righe che iniziano con un numero (ignora le intestazioni di Patran)
    if ~isempty(riga) && ~isnan(str2double(strtok(riga)))
        numeri = sscanf(riga, '%f')';
        if length(numeri) >= 3
            % Penultimo numero = Z Location, Ultimo numero = Valore Risultato
            dati = [dati; numeri(end-1), numeri(end)]; %#ok<AGROW>
        elseif length(numeri) == 2
            % Caso con solo [Elem_ID, Valore]
            dati = [dati; numeri(1), numeri(2)]; %#ok<AGROW>
        end
    end
end
fclose(fid);

if isempty(dati)
    error('Nessun dato numerico letto da "%s". Verifica il contenuto del file.', nome_file);
end

% Se la prima colonna è un ID elemento (> 2.0), ricostruisce z da 0 a 1.5 m
if max(dati(:, 1)) > 2.0
    val_sorted = dati(:, 2);
    z_sorted = linspace(0.015, 1.485, length(val_sorted))';
else
    [z_sorted, idx] = sort(dati(:, 1));
    val_sorted = dati(idx, 2);
end
end

%% 5. ANALISI TENSIONI TANGENZIALI TAU_XY ALL'INCASTRO (PANNELLO 2-3 E 3-4)
s = 2e-3;                             % Spessore pannelli [m]
coord_id = linspace(-a/2, a/2, 200);  % Coordinata trasversale da -0.15 m a +0.15 m

% Lettura dei file .rpt da Patran (all'incastro z2):
[x_alto, tau_alto_fem]     = leggi_rpt_patran('pannello_alto_tau_xy.rpt');
[y_destra, tau_destra_fem] = leggi_rpt_patran('pannello_destra_tau_xy.rpt');

% Valori analitici ideali della tensione tangenziale tau_xy [Pa]
% (Il segno del pannello 3-4 viene allineato al verso locale del pannello in Patran)
segno_34 = sign(mean(tau_destra_fem));
if segno_34 == 0, segno_34 = -1; end

tau_23_id = zeros(size(coord_id));                                    % 0 Pa sul pannello 2-3
tau_34_id = segno_34 * (F_tot / (2 * a * s)) * ones(size(coord_id));      % +-5e6 Pa sul pannello 3-4

%% GRAFICO 2A: ANDAMENTO TAU_XY NEL PANNELLO 2-3 (ALTO, ORIZZONTALE)
fig2a = figure('Name', 'Andamento tau_xy pannello 2-3', 'Color', 'w', 'Position', [150 150 850 500]);
if isprop(fig2a, 'Theme'), fig2a.Theme = 'light'; end

ax2a = axes(fig2a, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.75 0.75 0.75]);
hold(ax2a, 'on'); grid(ax2a, 'on'); box(ax2a, 'on');

plot(ax2a, x_alto, tau_alto_fem, '-', 'Color', [0.15 0.25 0.35], 'LineWidth', 1.8, 'DisplayName', 'FEM');
plot(ax2a, coord_id, tau_23_id, '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.8, 'DisplayName', 'Metodo trave a semiguscio');

xlabel(ax2a, 'Asse x [m]', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
ylabel(ax2a, '\tau_{xy} [MPa]', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
title(ax2a, 'Andamento \tau_{xy} nel pannello 2-3', 'FontSize', 13, 'Color', 'k');
xlim(ax2a, [-a/2, a/2]);

lgd2a = legend(ax2a, 'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 10);
set(lgd2a, 'Color', 'w', 'TextColor', 'k', 'EdgeColor', 'none');

exportgraphics(fig2a, 'plot_tau_pannello_2_3.png', 'Resolution', 300, 'BackgroundColor', 'white');

%% GRAFICO 2B: ANDAMENTO TAU_XY NEL PANNELLO 3-4 (DESTRA, VERTICALE)
fig2b = figure('Name', 'Andamento tau_xy pannello 3-4', 'Color', 'w', 'Position', [200 200 850 500]);
if isprop(fig2b, 'Theme'), fig2b.Theme = 'light'; end

ax2b = axes(fig2b, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.75 0.75 0.75]);
hold(ax2b, 'on'); grid(ax2b, 'on'); box(ax2b, 'on');

plot(ax2b, y_destra, tau_destra_fem, '-', 'Color', [0.15 0.25 0.35], 'LineWidth', 1.8, 'DisplayName', 'FEM');
plot(ax2b, coord_id, tau_34_id, '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.8, 'DisplayName', 'Metodo trave a semiguscio');

xlabel(ax2b, 'Asse y [m]', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
ylabel(ax2b, '\tau_{xy} [MPa]', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
title(ax2b, 'Andamento \tau_{xy} nel pannello 3-4', 'FontSize', 13, 'Color', 'k');
xlim(ax2b, [-a/2, a/2]);

lgd2b = legend(ax2b, 'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 10);
set(lgd2b, 'Color', 'w', 'TextColor', 'k', 'EdgeColor', 'none');

exportgraphics(fig2b, 'plot_tau_pannello_3_4.png', 'Resolution', 300, 'BackgroundColor', 'white');