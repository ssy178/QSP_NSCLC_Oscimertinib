%% Osimertinib Integrated PBPK Model — GI (Module 1) + Systemic (Module 2)
% Module 1: 3-compartment GI dissolution/absorption model (9 state variables)
% Module 2: Whole-body PBPK PS-limited distribution model (15 state variables)
% Total state variables: 24
%
% Integration core: A_gw(t) is dynamically computed by Module 1 ODE (replaces gamma approximation in Module 2)
%   - Module 2 portal vein ODE uses C_u_gw = f_u_gw * A_gw / V_gw
%   - A_gw = y(5) (gut wall compartment state variable in Module 1)
%
% File structure:
%   Osimertinib_PBPK_integrated_run.m  ← this file (main run script)
%   osimertinib_combined_odes.m        ← ODE function file (24 equations)

clear; clc; close all;


%% 1. Simulation Time and Initial Conditions
% 336h (14 days) = ~7×t½ (t½ ~48h) → sufficient to capture terminal phase
tspan = 0:0.5:336;   % 0 ~ 336 hours, 0.5-hour intervals [h]
Dose  = 160;         % Dose: 80 mg ≈ 160 μmol (MW 499.6)

% Initial state variables (24):
%   y(1–9):  Module 1 GI compartments
%     y(1) D_sto, y(2) A_sto, y(3) D_si, y(4) A_si, y(5) A_gw, y(6) A_col
%     y(7) A_portal_cum, y(8) A_met_cum, y(9) A_fecal_cum  ← cumulative tracking
%   y(10–24): Module 2 systemic compartments (concentration unit: μM)
%     y(10) Cart, y(11) Cven, y(12) Cli,    y(13) Cpv,    y(14) Cki
%     y(15) Clun, y(16) Cbr,
%     y(17) Cmu_v,  y(18) Cfat_v, y(19) Cre_v
%     y(20) Ctuv,   y(21) Ctuev
%     y(22) Cmu_t,  y(23) Cfat_t, y(24) Cre_t
y0 = zeros(24, 1);
y0(1) = Dose;   % D_sto: entire dose present as solid in stomach immediately after administration

%% 2. Parameter Setup
p = struct();

% =========================================================
% --- Module 1: GI Absorption Parameters ---
% =========================================================

% Rate constants
p.k_ge       = 0.45;  % Gastric emptying first-order rate constant (t1/2~1.5h; fitted to clinical Tmax~3–6h) [1/h]
p.k_diss_sto = 1.0;   % In vivo dissolution rate constant in stomach (acidic environment; reduced from 5.0 to suppress excessive early absorption) [1/h]
p.k_diss_si  = 0.15;  % In vivo dissolution rate constant in small intestine (near-neutral pH; reduced from 0.4 to delay Tmax) [1/h]
p.k_si       = 0.35;  % Small intestine transit rate constant (MSITT ≈ 3h; reduced from 0.7) [1/h]
p.k_fecal    = 0.03;  % Fecal excretion rate constant (fasted; colonic transit ~30–40 h) [1/h]
p.k_deg_sto  = 0.0;   % Gastric drug degradation rate constant (set to 0; osimertinib is acid-stable) [1/h]

% Physicochemical and physiological parameters
p.pKa1   = 4.44;   % Osimertinib first ionization constant (aniline nitrogen, weakly basic) [dimensionless]
p.pKa2   = 9.20;   % Osimertinib second ionization constant (aliphatic amine; experimental estimate range 8.5–9.5) [dimensionless]
p.pH_si  = 6.0;    % Small intestine (jejunum) luminal pH (fasted: 6.0–7.0) [dimensionless]
p.V_sto  = 0.250;  % Gastric fluid volume immediately after dosing (fasted ~35 mL + 200–240 mL water) [L]
p.V_si   = 0.116;  % Small intestinal water volume during absorption (CAT/ACAT model SIWV; range 0.08–0.17) [L]
p.V_gw   = 0.5;    % Gut wall (enterocyte) compartment volume (physiological range 0.4–0.6 L) [L]
p.f_u_gw = 0.01;   % Unbound drug fraction in gut wall (reflecting high lipophilicity and protein binding >99%) [dimensionless]
p.SAFE   = 45;     % Surface area expansion factor (folds/villi/microvilli amplification; fitted to F_abs ~80%) [dimensionless]

% Biopharmaceutical absorption parameters
p.CL_trans = 0.65;  % Transcellular permeability clearance (CL_trans = Peff × SA_gut; upper limit for F_abs calibration) [L/h]

% Michaelis-Menten parameters (gut wall)
p.V_max_gw  = 0.01; % Maximum rate of intestinal CYP3A4-mediated metabolism (Fg ≥ 0.95, near zero) [μmol/h]
p.K_m_gw    = 34.0; % Michaelis constant for CYP3A4 affinity to osimertinib (range 26–42) [μM]
p.V_max_Pgp = 0.5;  % Maximum rate of P-gp (ABCB1)-mediated intestinal efflux (saturable at clinical concentrations) [μmol/h]
p.K_m_Pgp   = 5.0;  % Michaelis constant for P-gp affinity to osimertinib (lipophilic TKI range 0.5–10) [μM]

% =========================================================
% --- Module 2: Systemic PBPK Parameters ---
% =========================================================

% Blood flows [L/h]
Q_co  = 312;             % Cardiac output [L/h]
Q_li  = 97;              % Total hepatic blood flow [L/h]
Q_ki  = 59;              % Renal blood flow [L/h]
Q_br  = 45;              % Cerebral blood flow [L/h]
Q_mu  = 47;              % Muscle blood flow [L/h]
Q_fat = 19;              % Adipose blood flow [L/h]
Q_tu  = 3.075;           % Tumour blood flow [L/h]

p.Q_lu_n = Q_co - Q_tu;          % Normal lung blood flow [L/h]
p.Q_HA   = 22;                   % Hepatic artery blood flow [L/h]
p.Q_pv   = Q_li - p.Q_HA;        % Portal vein blood flow [L/h] (~75 L/h); applied to both GI and systemic modules
p.Q_li   = Q_li;
p.Q_ki   = Q_ki;
p.Q_br   = Q_br;
p.Q_mu   = Q_mu;
p.Q_fat  = Q_fat;
p.Q_tu   = Q_tu;
p.Q_co   = Q_co;
p.Q_re   = Q_co - (Q_li + Q_ki + Q_br + Q_mu + Q_fat + Q_tu); % [L/h]

% Volumes [L]
p.V_art   = 1.4;         % Arterial blood [L]
p.V_ven   = 3.0;         % Venous blood [L]
p.V_li    = 1.5;         % Liver tissue [L]
p.V_pv    = 0.075;       % Portal vein blood volume [L]
p.V_ki    = 0.28;        % Kidney tissue [L]
p.V_lu_n  = 0.53;        % Effective normal lung tissue volume, not total anatomical lung volume [L]
p.V_br    = 1.45;        % Brain tissue [L]
p.V_tu_v  = 0.005;       % Tumour vascular space [L]
p.V_tu_ev = 0.095;       % Tumour extravascular space [L]
V_total   = 70;          % Total body volume [L]

% Total volumes for muscle, fat, rest-of-body
p.V_mu    = 29;          % Muscle total tissue [L]
p.V_fat   = 15;          % Adipose total tissue [L]
p.V_re    = V_total - (p.V_art + p.V_ven + p.V_li + p.V_pv + p.V_ki + ...
                        p.V_lu_n + p.V_br + p.V_mu + p.V_fat + p.V_tu_v + p.V_tu_ev);

% PS-limited subcompartment split (vascular fraction fv = 0.04)
% Rationale: muscle/fat/rest have PS << Q → permeability-limited
%   τ_eq = V_t * Kp / PS ≈ 32h (>> Tmax ~6h) → slow tissue filling
%   This brings Cmax_blood close to clinical values and prolongs terminal t½
fv = 0.04;                           % vascular fraction [-]
p.V_mu_v  = fv * p.V_mu;            % Muscle vascular [L]  = 1.16 L
p.V_mu_t  = (1 - fv) * p.V_mu;      % Muscle tissue   [L]  = 27.84 L
p.V_fat_v = fv * p.V_fat;           % Adipose vascular [L] = 0.60 L
p.V_fat_t = (1 - fv) * p.V_fat;     % Adipose tissue  [L]  = 14.40 L
p.V_re_v  = fv * p.V_re;            % Rest vascular   [L]
p.V_re_t  = (1 - fv) * p.V_re;      % Rest tissue     [L]

% Partition coefficients (Kp) and unbound fractions (fu)
p.R_bp    = 0.917;       % Blood-to-plasma ratio [-]
p.Kp_lu_n = 8.53;        % Lung partition [-]
p.Kp_li   = 11;          % Liver partition [-]
p.Kp_ki   = 4.5;         % Kidney partition [-]
p.Kp_br   = 0.8;         % Brain partition [-]
% Kp_mu/fat/re: maintain Vd_ss ~840 L (target 640–918 L) while reinforcing terminal tail
% Previous: Kp_mu=10, Kp_fat=20, Kp_re=6 → τ_eq ~33h, t½ too short
% Updated: increased tissue storage → τ_eq ~84–90h → supports terminal t½ ~48h
p.Kp_mu   = 12.0;        % Muscle partition   [-]  (previously 10.0; increased τ_eq)
p.Kp_fat  = 25.0;        % Adipose partition  [-]  (previously 20.0; better reflects high logP ~4.3)
p.Kp_re   = 7.0;         % Rest-of-body partition [-] (previously 6.0)

p.fu_p  = 0.036;         % Fraction unbound in plasma [-]
p.fu_li = 0.0055;        % Fraction unbound in liver [-]
p.fu_br = 0.008;         % Fraction unbound in brain [-]
p.fu_ki = 0.036;         % Fraction unbound in kidney [-]
p.fu_t  = 0.0255;        % Fraction unbound in tumour [-]

% Clearance and permeability [L/h]
% CL_h ~9.9 L/h (clinical CL/F=14.2, F=0.698): back-calculated from well-stirred model
p.CL_int_total = 2000;   % Total hepatic intrinsic clearance [L/h]
p.CL_int_sec   = 10;     % Renal secretory clearance [L/h]
p.GFR          = 7.2;    % Glomerular filtration rate [L/h]
p.k_efflux_BBB = 1.05;   % BBB efflux rate [1/h]
p.PS_tu        = 0.2505; % Tumour PS [L/h]
p.k_int        = 5.0;    % Tumour internalization rate [1/h]

% PS-limited parameters: PS << Q → slow tissue equilibration
% Previous: τ_eq ~33h → tissue filled and emptied too quickly, insufficient support for terminal tail
% Updated: reduced PS to extend τ_eq to ~84–90h, supporting terminal t½ ~48h
% τ_eq(muscle) = V_mu_t * Kp_mu / PS_mu = 27.84*12/4.0 ≈ 83.5 h
% τ_eq(fat)    = V_fat_t * Kp_fat / PS_fat = 14.4*25/4.0 ≈ 90.0 h
% τ_eq(rest)   = V_re_t * Kp_re / PS_re ≈ 16.96*7/1.5 ≈ 79.1 h
p.PS_mu  = 4.0;          % Muscle permeability-surface area [L/h]  (Q_mu=47 → PS/Q=0.085, deeply limited)
p.PS_fat = 4.0;          % Adipose PS [L/h]                        (Q_fat=19 → PS/Q=0.21, limited)
p.PS_re  = 1.5;          % Rest-of-body PS [L/h]                   (Q_re≈42 → PS/Q=0.036, deeply limited)

%% 3. Run ODE Solver
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);  % Strict numerical tolerances (minimize mass balance error)
[t, y] = ode15s(@(t, y) osimertinib_combined_odes(t, y, p), tspan, y0, opts);

%% 4. Post-processing: Module 1 — Mass Balance and Absorption Fraction
A_portal_cum = y(:,7);  % Cumulative portal absorption (before hepatic first-pass) [μmol]
A_met_cum    = y(:,8);  % Cumulative gut wall CYP3A4 metabolism [μmol]
A_fecal_cum  = y(:,9);  % Cumulative fecal excretion [μmol]

% F_abs: intestinal absorption fraction reaching portal vein (before hepatic first-pass)
F_abs_pct = A_portal_cum / Dose * 100;  % [%]

% Mass balance verification: sum of all compartments must equal Dose at all times
GI_total   = y(:,1) + y(:,2) + y(:,3) + y(:,4) + y(:,5) + y(:,6);
mass_total = GI_total + A_portal_cum + A_met_cum + A_fecal_cum;
mass_error_pct = (mass_total - Dose) / Dose * 100;  % [%]

% Console summary output (GI)
% F_abs: intestinal absorption fraction reaching portal vein (clinical reference: portal absorption ~78–80%)
% F_oral (systemic bioavailability) is reported in Section 5 PK Summary after F_h correction
fprintf('=== GI Absorption Summary at t = %.0f h ===\n', t(end));
fprintf('  Portal absorbed  (F_abs):   %6.2f μmol  (%5.1f%% of dose)\n', ...
        A_portal_cum(end), F_abs_pct(end));
fprintf('  Gut wall metabolism:        %6.2f μmol  (%5.1f%% of dose)\n', ...
        A_met_cum(end), A_met_cum(end)/Dose*100);
fprintf('  Fecal excretion:            %6.2f μmol  (%5.1f%% of dose)\n', ...
        A_fecal_cum(end), A_fecal_cum(end)/Dose*100);
fprintf('  Remaining in GI tract:      %6.2f μmol  (%5.1f%% of dose)\n', ...
        GI_total(end), GI_total(end)/Dose*100);
fprintf('  Mass balance error (max):   %+.5f%%\n', max(abs(mass_error_pct)));

%% 5. Post-processing: Module 2 — PK Metric Calculation and Console Output

% --- Arterial blood PK ---
Cart_vec = y(:,10);
[Cmax_art, idx_tmax] = max(Cart_vec);
tmax_art     = t(idx_tmax);
AUC_art      = trapz(t, Cart_vec);
Cmax_plasma  = Cmax_art / p.R_bp;

% Terminal t½: log-linear fit over t ≥ 72h window
%   Rationale: fully excludes distribution phase (Tmax~3–6h) and early decline
%   Previous approach (last 40% of concentrations above 0.1% Cmax) could misidentify
%   distribution phase as terminal phase, yielding unrealistic values (e.g., t½ = 0.8h)
t_term_start = 72;   % terminal window start time [h] (after distribution phase)
idx_term = find(t >= t_term_start & Cart_vec > 1e-6 * Cmax_art);
if numel(idx_term) > 5
    pfit   = polyfit(t(idx_term), log(max(Cart_vec(idx_term), 1e-15)), 1);
    k_el   = -pfit(1);
    t_half = log(2) / max(k_el, 1e-6);
else
    t_half = NaN;
end

% Kp_uu_br (steady-state analytic)
Kpuu_br = (p.Q_br * p.R_bp * p.fu_br) / ...
          (p.fu_p * (p.Q_br + p.k_efflux_BBB * p.fu_br * p.V_br));

% Hepatic CL (well-stirred model)
CL_h = p.Q_li * (p.CL_int_total * p.fu_li) / ...
       (p.Q_li + p.CL_int_total * p.fu_li);

% Vd_ss: blood volumes + tissue volumes × Kp
% PS-limited: V_mu_v + V_mu_t*Kp_mu = V_mu*(fv + (1-fv)*Kp_mu) ≈ V_mu*Kp_mu
Vd_ss = p.V_art + p.V_ven + p.V_pv + ...
        p.V_lu_n * p.Kp_lu_n + p.V_li * p.Kp_li + ...
        p.V_ki   * p.Kp_ki   + p.V_br * p.Kp_br + ...
        p.V_mu_v + p.V_mu_t  * p.Kp_mu + ...
        p.V_fat_v + p.V_fat_t * p.Kp_fat + ...
        p.V_re_v  + p.V_re_t  * p.Kp_re  + ...
        p.V_tu_v  + p.V_tu_ev;

% F_oral (systemic fraction after hepatic first-pass) — well-stirred analytical calculation
%   F_oral = F_abs (intestinal absorption fraction) × F_h (hepatic first-pass survival fraction)
%   F_h = 1 - E_h,  E_h = CL_h / Q_li
%   Previous approach (∫Q_li·C_li/Kp_li dt / Dose) overestimates F_oral by including
%   recirculating drug, so replaced with the analytical formula
F_h   = 1 - CL_h / p.Q_li;                   % Hepatic first-pass survival fraction [-]
F_sim = (A_portal_cum(end) / Dose) * F_h;    % F_abs × F_h [-]

fprintf('\n========= Simulation PK Summary (PS-limited) =========\n');
fprintf('  Dose                         : %6.1f μmol  (80 mg)\n', Dose);
fprintf('  Blood  Cmax  (C_art)         : %6.3f μM\n',  Cmax_art);
fprintf('  Plasma Cmax  (C_art/R_bp)    : %6.3f μM   [target ~0.5–1.07 μM]\n', Cmax_plasma);
fprintf('  Blood  Tmax  (C_art)         : %6.1f h    [target ~3–6 h]\n',  tmax_art);
fprintf('  AUC_0-%.0fh (blood)          : %6.1f μM·h\n', t(end), AUC_art);
fprintf('  Apparent terminal t½         : %6.1f h    [target ~48 h]\n',   t_half);
fprintf('  Hepatic CL (well-stirred)    : %6.2f L/h  [target ~9.9 L/h]\n', CL_h);
fprintf('  Vd_ss (calculated)           : %6.0f L    [target ~640–918 L]\n', Vd_ss);
fprintf('  F_oral = F_abs×F_h           : %6.3f      [target ~0.698; F_h=%.3f]\n', F_sim, F_h);
fprintf('  Kp_uu_br (analytic)          : %6.3f      [target ~0.21]\n',   Kpuu_br);
fprintf('  PS_mu / Q_mu                 : %6.3f      (permeability-limited if < 0.3)\n', p.PS_mu/p.Q_mu);
fprintf('  τ_eq,muscle                  : %6.1f h    (V_mu_t*Kp_mu/PS_mu)\n', p.V_mu_t*p.Kp_mu/p.PS_mu);
fprintf('  τ_eq,fat                     : %6.1f h    (V_fat_t*Kp_fat/PS_fat)\n', p.V_fat_t*p.Kp_fat/p.PS_fat);
fprintf('  τ_eq,rest                    : %6.1f h    (V_re_t*Kp_re/PS_re)\n',  p.V_re_t*p.Kp_re/p.PS_re);
fprintf('  Terminal fit window          : t ≥ %.0f h\n', t_term_start);
fprintf('======================================================\n\n');

%% 6. Visualization — Figure 1: GI Absorption Dynamics (3×2 layout)

% Force dark background (MATLAB dark theme)
set(groot, 'defaultFigureColor',  [0.15 0.15 0.15], ...  % Figure background
           'defaultAxesColor',    [0.15 0.15 0.15], ...  % Axes background
           'defaultAxesXColor',   'w', ...                % X-axis ticks/labels
           'defaultAxesYColor',   'w', ...                % Y-axis ticks/labels
           'defaultTextColor',    'w');                   % Title/legend text

figure('Position', [80, 80, 1100, 900]);

% --- (1) GI lumen: stomach + small intestine (displayed 0–10h only) ---
subplot(3,2,[1 2]);
plot(t, y(:,1), '--', 'Color', [1.00 0.85 0.00], 'LineWidth', 2.0); hold on;  % D_sto: golden yellow
plot(t, y(:,3), '--', 'Color', [0.00 0.90 0.90], 'LineWidth', 2.0);           % D_si:  cyan
plot(t, y(:,2), 'r',  'LineWidth', 2);                                         % A_sto: red
plot(t, y(:,4), 'b',  'LineWidth', 2);                                         % A_si:  blue
xlim([0 10]);
title('GI Lumen — Stomach & Small Intestine  (0–10 h)');
xlabel('Time (h)'); ylabel('Amount (\mu mol)');
legend('D_{sto} (Solid)', 'D_{si} (Solid)', 'A_{sto} (Dissolved)', 'A_{si} (Dissolved)', ...
       'Location', 'northeast');
grid on;

% --- (2) Gut wall (A_gw) — dynamic ODE result from Module 1 (not gamma approximation) ---
subplot(3,2,3);
plot(t, y(:,5), 'm', 'LineWidth', 2);
title('Gut Wall (A_{gw})  — Dynamic ODE (Module 1)');
xlabel('Time (h)'); ylabel('Amount (\mu mol)');
xlim([0 48]);
grid on;

% --- (3) Colon (A_col) ---
subplot(3,2,4);
plot(t, y(:,6), 'g', 'LineWidth', 2);
title('Colon (A_{col})');
xlabel('Time (h)'); ylabel('Amount (\mu mol)');
xlim([0 48]);
grid on;

% --- (4) Cumulative portal absorption + F_abs (dual y-axis) ---
subplot(3,2,5);
yyaxis left;
plot(t, A_portal_cum, 'b', 'LineWidth', 2);
ylabel('Cumulative Amount (\mu mol)');
ylim([0 Dose]);
yyaxis right;
plot(t, F_abs_pct, 'r--', 'LineWidth', 1.5);
yline(80, 'k:', 'LineWidth', 1);
text(28, 82, 'F=80% (overall, clinical ref.)', 'FontSize', 8);
ylabel('F_{abs} (%)');
ylim([0 100]);
title('Cumulative Portal Absorption & F_{abs}');
xlabel('Time (h)');
xlim([0 48]);
grid on;

% --- (5) Mass balance check (stacked area) ---
subplot(3,2,6);
area(t, [GI_total/Dose*100, ...
         A_portal_cum/Dose*100, ...
         A_met_cum/Dose*100, ...
         A_fecal_cum/Dose*100], ...
    'FaceAlpha', 0.65);
yline(100, 'k--', 'LineWidth', 1.5);
legend('GI Tract', 'Portal (Absorbed)', 'Gut Wall Met.', 'Fecal', ...
       'Location', 'northeast', 'FontSize', 8);
title(sprintf('Mass Balance  (max error: %.4f%%)', max(abs(mass_error_pct))));
xlabel('Time (h)'); ylabel('% of Dose');
xlim([0 48]);
ylim([0 115]);
grid on;

sgtitle('Osimertinib Integrated PBPK — GI Absorption (Module 1, Mass Balance Verified)', ...
        'FontSize', 14, 'FontWeight', 'bold');

%% 7. Visualization — Figure 2: A_gw Absorption Dynamics Diagnostics (3×1 layout)

Agw_ode    = y(:,5);                                    % Actual A_gw [μmol] (Module 1 ODE)
R_abs_ode  = p.Q_pv * p.f_u_gw / p.V_gw * Agw_ode;    % Portal absorption rate [μmol/h]
Cum_abs    = cumtrapz(t, R_abs_ode);                    % Cumulative portal absorption [μmol]

figure('Position', [100, 100, 900, 700]);

subplot(3,1,1);
plot(t, Agw_ode, 'Color', [1.00 0.75 0.30], 'LineWidth', 2.5);
[Agw_max, idx_peak] = max(Agw_ode);
xline(t(idx_peak), 'w--', 'LineWidth', 1.2);
text(t(idx_peak)+0.4, Agw_max*0.92, ...
     sprintf('T_{peak}=%.1f h\nA_{gw,max}=%.1f \\mumol', t(idx_peak), Agw_max), ...
     'Color','w', 'FontSize', 9);
ylabel('A_{gw} (\mumol)', 'FontSize', 10);
title('Gut Wall Drug Amount — A_{gw}(t)  [Module 1 ODE, Dynamic]', 'FontSize', 10, 'FontWeight', 'bold');
xlim([0 48]); grid on;
set(gca, 'GridColor',[0.4 0.4 0.4], 'GridAlpha', 0.5);

subplot(3,1,2);
plot(t, R_abs_ode, 'Color', [0.30 0.85 1.00], 'LineWidth', 2.5);
[R_max, idx_r] = max(R_abs_ode);
xline(t(idx_r), 'w--', 'LineWidth', 1.2);
text(t(idx_r)+0.4, R_max*0.92, ...
     sprintf('T_{peak}=%.1f h\nR_{max}=%.1f \\mumol/h', t(idx_r), R_max), ...
     'Color','w', 'FontSize', 9);
ylabel('R_{abs} (\mumol/h)', 'FontSize', 10);
title('Absorption Flux into Portal Vein — R_{abs}(t) = Q_{pv}\cdotf_{u,gw}\cdotA_{gw}/V_{gw}', ...
      'FontSize', 10, 'FontWeight', 'bold');
xlim([0 48]); grid on;
set(gca, 'GridColor',[0.4 0.4 0.4], 'GridAlpha', 0.5);

subplot(3,1,3);
yyaxis left
plot(t, Cum_abs, 'Color', [0.40 1.00 0.40], 'LineWidth', 2.5);
ylabel('Cumulative Absorption (\mumol)', 'FontSize', 10);
ylim([0 Dose]);
yyaxis right
plot(t, Cum_abs/Dose*100, 'Color', [0.40 1.00 0.40], 'LineWidth', 1.0);
yline(78, 'w--', 'LineWidth', 1.2);
text(30, 80, sprintf('F_{abs}\\approx%.0f%%  (total input)', Cum_abs(end)/Dose*100), ...
     'Color','w','FontSize',9);
ylabel('F_{abs} (%)', 'FontSize', 10);
ylim([0 100]);
xlabel('Time (h)', 'FontSize', 10);
title('Cumulative Portal Absorption', 'FontSize', 10, 'FontWeight', 'bold');
xlim([0 48]); grid on;
set(gca, 'GridColor',[0.4 0.4 0.4], 'GridAlpha', 0.5);

sgtitle(sprintf('A_{gw}(t) from Module 1 ODE  [peak at t=%.1f h,  F_{abs}\\approx%.0f%%]', ...
                t(idx_peak), Cum_abs(end)/Dose*100), ...
        'FontSize', 11, 'FontWeight', 'bold');

%% 8. Visualization — Figure 3: Systemic PBPK (5×3 layout, 15 compartments)
% State variable labels and colors (corresponding to y(10)–y(24))
labels = { 'Arterial Blood (C_{art})',       'Venous Blood (C_{ven})',              ...
           'Liver (C_{li})',                 'Portal Vein (C_{pv})',                ...
           'Kidney (C_{ki})',                'Normal Lung (C_{lu})',                ...
           'Brain (C_{br})',                 'Muscle — Vascular (C_{mu,v})',        ...
           'Adipose — Vascular (C_{fat,v})', 'Rest-of-Body — Vascular (C_{re,v})', ...
           'Tumour Vascular (C_{tu,v})',     'Tumour EV (C_{tu,ev})',               ...
           'Muscle — Tissue (C_{mu,t})',     'Adipose — Tissue (C_{fat,t})',        ...
           'Rest-of-Body — Tissue (C_{re,t})' };

clrs = { [1.00 0.25 0.25],  [0.30 0.55 1.00], ...
         [1.00 0.60 0.00],  [1.00 0.85 0.00], ...
         [0.00 0.90 0.90],  [0.40 1.00 0.40], ...
         [1.00 0.00 1.00],  [1.00 0.50 0.70], ...
         [1.00 0.75 0.45],  [0.70 0.70 0.70], ...
         [0.20 1.00 0.20],  [0.80 0.40 1.00], ...
         [0.85 0.35 0.55],  [0.95 0.65 0.30], ...
         [0.55 0.55 0.55]  };

ref_Cmax_blood = 1.07 * p.R_bp;   % plasma 1.07 μM × R_bp → blood
ref_Tmax       = 6;               % h
ref_thalf      = 48;              % h

figure('Position', [60, 30, 1500, 1000]);

for i = 1:15
    subplot(5, 3, i);
    plot(t, y(:, 9+i), 'Color', clrs{i}, 'LineWidth', 2.0); hold on;

    if i == 1
        yline(ref_Cmax_blood, 'w:', 'LineWidth', 1.2);
        xline(ref_Tmax,       'w:', 'LineWidth', 1.0);
        text(ref_Tmax + 5, ref_Cmax_blood * 1.08, ...
             sprintf('Ref C_{max}=%.2f\\muM', ref_Cmax_blood), ...
             'Color','w','FontSize',7);
    end
    if i == 1 || i == 2
        [Cm, ~] = max(y(:, 9+i));
        yline(Cm * 0.5, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);
        xline(ref_thalf, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);
    end

    title(labels{i}, 'FontSize', 8, 'FontWeight', 'bold');
    xlabel('Time (h)', 'FontSize', 7);
    ylabel('Conc. (\muM)', 'FontSize', 7);
    xlim([0 tspan(end)]);
    grid on;
    set(gca, 'GridColor', [0.4 0.4 0.4], 'GridAlpha', 0.5);
end

sgtitle(sprintf(['Osimertinib Whole-Body PBPK — PS-Limited Distribution (%.0f h)\n' ...
                 'C_{max,plasma}=%.3f\\muM  T_{max}=%.0fh  t_{1/2}=%.0fh  ' ...
                 'Vd_{ss}=%.0fL  K_{p,uu,br}=%.2f'], ...
                 t(end), Cmax_plasma, tmax_art, t_half, Vd_ss, Kpuu_br), ...
        'FontSize', 10, 'FontWeight', 'bold');
