function dydt = osimertinib_combined_odes(~, y, p)
% osimertinib_combined_odes — Integrated PBPK ODE function
% Module 1 (GI, 9 variables) + Module 2 (Systemic PBPK, 15 variables) = 24 state variables total
%
% Integration core:
%   A_gw = y(5) (dynamically computed by Module 1) → directly linked to Module 2 portal vein ODE
%   C_u_gw = f_u_gw * A_gw / V_gw: fully replaces the gamma approximation in Module 2
%
% State variable array:
%   y(1–9):   Module 1 — GI compartments
%     y(1) D_sto        : Solid drug in stomach [μmol]
%     y(2) A_sto        : Dissolved drug in stomach [μmol]
%     y(3) D_si         : Solid drug in small intestine [μmol]
%     y(4) A_si         : Dissolved drug in small intestine [μmol]
%     y(5) A_gw         : Drug in gut wall (enterocyte) [μmol]  ← connection point to Module 2
%     y(6) A_col        : Drug in colon [μmol]
%     y(7) A_portal_cum : Cumulative portal vein absorption [μmol]
%     y(8) A_met_cum    : Cumulative gut wall metabolism [μmol]
%     y(9) A_fecal_cum  : Cumulative fecal excretion [μmol]
%
%   y(10–24): Module 2 — Systemic compartments (concentration unit: μM)
%     y(10) Cart   : Arterial blood
%     y(11) Cven   : Venous blood
%     y(12) Cli    : Liver
%     y(13) Cpv    : Portal vein
%     y(14) Cki    : Kidney
%     y(15) Clun   : Normal lung
%     y(16) Cbr    : Brain
%     y(17) Cmu_v  : Muscle — vascular (PS-limited)
%     y(18) Cfat_v : Adipose — vascular (PS-limited)
%     y(19) Cre_v  : Rest-of-body — vascular (PS-limited)
%     y(20) Ctuv   : Tumour vascular
%     y(21) Ctuev  : Tumour extravascular
%     y(22) Cmu_t  : Muscle tissue (PS-limited)
%     y(23) Cfat_t : Adipose tissue (PS-limited)
%     y(24) Cre_t  : Rest-of-body tissue (PS-limited)

dydt = zeros(24, 1);

%% ===== Module 1: GI Compartments (y(1)–y(9)) =====
D_sto = y(1); A_sto = y(2);
D_si  = y(3); A_si  = y(4);
A_gw  = y(5); A_col = y(6);
% y(7–9): cumulative tracking variables (A_portal_cum, A_met_cum, A_fecal_cum)

% --- Algebraic equations ---
% Unionized fraction at small intestine pH (diprotic weak base)
f_un_ph = 1 / (1 + 10^(p.pKa1 - p.pH_si) + 10^(p.pKa1 + p.pKa2 - 2*p.pH_si));

% Concentration calculation in small intestine and gut wall (prevent negative values)
C_si   = max(0, A_si / p.V_si);
C_gw   = max(0, A_gw / p.V_gw);
C_u_gw = p.f_u_gw * C_gw;   % Unbound drug concentration in gut wall [μM] — shared with Module 2

% Absorption rate: CL_trans × SAFE × unionized fraction × small intestine concentration
v_abs        = p.CL_trans * p.SAFE * f_un_ph * C_si;
v_met_gw     = (p.V_max_gw  * C_u_gw) / (p.K_m_gw  + C_u_gw);
v_efflux_Pgp = (p.V_max_Pgp * C_u_gw) / (p.K_m_Pgp + C_u_gw);
v_portal     = p.Q_pv * C_u_gw;   % Gut wall → portal vein transfer rate [μmol/h]

% --- ODE: 6 main compartments ---
dydt(1) = -p.k_diss_sto * D_sto - p.k_ge * D_sto;
dydt(2) =  p.k_diss_sto * D_sto - p.k_ge * A_sto - p.k_deg_sto * A_sto;
dydt(3) =  p.k_ge * D_sto - (p.k_diss_si + p.k_si) * D_si;
dydt(4) =  p.k_ge * A_sto + p.k_diss_si * D_si - v_abs - p.k_si * A_si + v_efflux_Pgp;
dydt(5) =  v_abs - v_met_gw - v_portal - v_efflux_Pgp;
dydt(6) =  p.k_si * (D_si + A_si) - p.k_fecal * A_col;

% --- ODE: 3 cumulative tracking variables ---
dydt(7) = v_portal;            % Cumulative portal vein absorption
dydt(8) = v_met_gw;            % Cumulative gut wall metabolism
dydt(9) = p.k_fecal * A_col;  % Cumulative fecal excretion

%% ===== Module 2: Systemic PBPK Compartments (y(10)–y(24)) =====
Cart   = y(10);  Cven   = y(11);  Cli    = y(12);  Cpv    = y(13);
Cki    = y(14);  Clun   = y(15);  Cbr    = y(16);
Cmu_v  = y(17);  Cfat_v = y(18);  Cre_v  = y(19);
Ctuv   = y(20);  Ctuev  = y(21);
Cmu_t  = y(22);  Cfat_t = y(23);  Cre_t  = y(24);

% Unbound concentration calculation
Cubr = p.fu_br * (Cbr / p.Kp_br);   % Unbound brain conc. [μM]
Culi = p.fu_li * Cli;                % Unbound liver conc. [μM]
Cuki = p.fu_ki * (Cki / p.Kp_ki);   % Unbound kidney conc. [μM]

% PS transport driving forces (unbound plasma conc. vs unbound tissue conc.)
% J = PS * (fu_p * C_vascular/R_bp  -  C_tissue/Kp)
J_mu  = p.PS_mu  * (p.fu_p * Cmu_v  / p.R_bp - Cmu_t  / p.Kp_mu);
J_fat = p.PS_fat * (p.fu_p * Cfat_v / p.R_bp - Cfat_t / p.Kp_fat);
J_re  = p.PS_re  * (p.fu_p * Cre_v  / p.R_bp - Cre_t  / p.Kp_re);

% 10. Arterial Blood
dydt(10) = (p.Q_lu_n * (Clun / p.Kp_lu_n) + p.Q_tu * Ctuv - p.Q_co * Cart) / p.V_art;

% 11. Venous Blood
%   PS-limited compartments return blood at vascular concentration (Cmu_v, not Cmu/Kp)
dydt(11) = (p.Q_li  * (Cli  / p.Kp_li)  + ...
            p.Q_ki  * (Cki  / p.Kp_ki)  + ...
            p.Q_br  * (Cbr  / p.Kp_br)  + p.k_efflux_BBB * Cubr * p.V_br + ...
            p.Q_mu  * Cmu_v  + ...
            p.Q_fat * Cfat_v + ...
            p.Q_re  * Cre_v  - ...
            p.Q_lu_n * Cven) / p.V_ven;

% 12. Liver
dydt(12) = (p.Q_HA * Cart + p.Q_pv * Cpv - p.Q_li * (Cli / p.Kp_li) - ...
            p.CL_int_total * Culi) / p.V_li;

% 13. Portal Vein
%   Drug inflow: mesenteric blood flow (Q_pv*Cart) + gut wall unbound drug (Q_pv*C_u_gw)
%   C_u_gw = f_u_gw * A_gw / V_gw: A_gw = y(5) is dynamically computed by Module 1
%   (fully replaces the old gamma approximation: 27.7*(t/1.5)^2*exp(-t/1.5))
dydt(13) = (p.Q_pv * Cart + p.Q_pv * C_u_gw - p.Q_pv * Cpv) / p.V_pv;

% 14. Kidney
dydt(14) = (p.Q_ki * Cart - p.Q_ki * (Cki / p.Kp_ki) - ...
            p.fu_p * p.GFR * Cart / p.R_bp - p.CL_int_sec * Cuki) / p.V_ki;

% 15. Normal Lung
dydt(15) = (p.Q_lu_n * Cven - p.Q_lu_n * (Clun / p.Kp_lu_n)) / p.V_lu_n;

% 16. Brain
dydt(16) = (p.Q_br * Cart - p.Q_br * (Cbr / p.Kp_br) - ...
            p.k_efflux_BBB * Cubr * p.V_br) / p.V_br;

% 17. Muscle — Vascular (PS-limited)
%   In: Q_mu*Cart;  Out: Q_mu*Cmu_v + PS transport to tissue
dydt(17) = (p.Q_mu * Cart - p.Q_mu * Cmu_v - J_mu) / p.V_mu_v;

% 18. Adipose — Vascular (PS-limited)
dydt(18) = (p.Q_fat * Cart - p.Q_fat * Cfat_v - J_fat) / p.V_fat_v;

% 19. Rest-of-Body — Vascular (PS-limited)
dydt(19) = (p.Q_re * Cart - p.Q_re * Cre_v - J_re) / p.V_re_v;

% 20. Tumour Vascular
dydt(20) = (p.Q_tu * (Cart - Ctuv) - ...
            p.PS_tu * (p.fu_p * Ctuv / p.R_bp - p.fu_t * Ctuev)) / p.V_tu_v;

% 21. Tumour Extravascular
dydt(21) = (p.PS_tu * (p.fu_p * Ctuv / p.R_bp - p.fu_t * Ctuev) - ...
            p.k_int * Ctuev * p.V_tu_ev) / p.V_tu_ev;

% 22. Muscle — Tissue (PS-limited)
dydt(22) = J_mu / p.V_mu_t;

% 23. Adipose — Tissue (PS-limited)
dydt(23) = J_fat / p.V_fat_t;

% 24. Rest-of-Body — Tissue (PS-limited)
dydt(24) = J_re / p.V_re_t;

end
