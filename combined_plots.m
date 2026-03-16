% Simple VCRS vs Ejector-based VCRS Performance Analysis 
clearvars; clc;
refrigerant = 'R1234yf'; 
T_evap_array = [275.15, 277.15, 279.15, 281.15, 283.15]; 
T_cond = 308.15; 
Q_target = 3517; 
eta_n = 0.85; 
eta_isen = 0.85; 
n = length(T_evap_array);
P_cond = py.CoolProp.CoolProp.PropsSI('P','T',T_cond,'Q',0,refrigerant);

COP_ej = zeros(1, n); COP_vcrs = zeros(1, n); 
s1_vcrs_arr = zeros(1,n); s2_vcrs_arr = zeros(1,n); s3_vcrs_arr = zeros(1,n); s4_vcrs_arr = zeros(1,n);
s1_ej_arr = zeros(1,n); s2_ej_arr = zeros(1,n); s3_ej_arr = zeros(1,n); s4_ej_arr = zeros(1,n);
s6_ej_arr = zeros(1,n); s7_ej_arr = zeros(1,n); s8_ej_arr = zeros(1,n); s9_ej_arr = zeros(1,n); s10_ej_arr = zeros(1,n);

for i = 1:n
    T_evap = T_evap_array(i); P_evap = py.CoolProp.CoolProp.PropsSI('P','T',T_evap,'Q',1,refrigerant);
    h1_vcrs = py.CoolProp.CoolProp.PropsSI('H','P',P_evap,'Q',1,refrigerant);
    s1_vcrs = py.CoolProp.CoolProp.PropsSI('S','P',P_evap,'Q',1,refrigerant);
    h2s_vcrs = py.CoolProp.CoolProp.PropsSI('H','P',P_cond,'S',s1_vcrs,refrigerant);
    h2_vcrs = h1_vcrs + (h2s_vcrs - h1_vcrs) / eta_isen;
    s2_vcrs = py.CoolProp.CoolProp.PropsSI('S','P',P_cond,'H',h2_vcrs,refrigerant);
    h3 = py.CoolProp.CoolProp.PropsSI('H','P',P_cond,'Q',0,refrigerant);
    s3 = py.CoolProp.CoolProp.PropsSI('S','P',P_cond,'Q',0,refrigerant);
    h4_vcrs = h3; s4_vcrs = py.CoolProp.CoolProp.PropsSI('S','P',P_evap,'H',h4_vcrs,refrigerant);
    
    Q_spec_vcrs = h1_vcrs - h4_vcrs; W_spec_vcrs = h2_vcrs - h1_vcrs; COP_vcrs(i) = Q_spec_vcrs / W_spec_vcrs;
    s1_vcrs_arr(i)=s1_vcrs; s2_vcrs_arr(i)=s2_vcrs; s3_vcrs_arr(i)=s3; s4_vcrs_arr(i)=s4_vcrs;
    
    h4s = py.CoolProp.CoolProp.PropsSI('H','P',P_evap,'S',s3,refrigerant);
    h4 = h3 - eta_n * (h3 - h4s); s4_ej = py.CoolProp.CoolProp.PropsSI('S','P',P_evap,'H',h4,refrigerant);
    u4 = sqrt(2 * (h3 - h4)); h10 = py.CoolProp.CoolProp.PropsSI('H','T',T_evap,'Q',1,refrigerant);
    s10 = py.CoolProp.CoolProp.PropsSI('S','T',T_evap,'Q',1,refrigerant);
    
    mu_low = 0.01; mu_high = 1.0;
    for iter = 1:100
        mu = (mu_low + mu_high) / 2; u6 = u4 / (1 + mu); 
        h6 = ((h4 + 0.5*u4^2) + mu*h10) / (1 + mu) - 0.5*u6^2;
        s6 = py.CoolProp.CoolProp.PropsSI('S','P',P_evap,'H',h6,refrigerant);
        h7 = (h3 + mu*h10) / (1 + mu); h7s = h6 + eta_n * (h7 - h6); 
        P_c = py.CoolProp.CoolProp.PropsSI('P','S',s6,'H',h7s,refrigerant);
        try x7 = py.CoolProp.CoolProp.PropsSI('Q','P',P_c,'H',h7,refrigerant); if x7 < 0, x7 = 0; elseif x7 > 1, x7 = 1; end
        catch H_sat_vap = py.CoolProp.CoolProp.PropsSI('H','P',P_c,'Q',1,refrigerant); x7 = (h7 >= H_sat_vap); end
        residual = (1 + mu) * x7 - 1;
        if residual > 0, mu_high = mu; else, mu_low = mu; end; if abs(residual) < 1e-4, break; end
    end
    s7 = py.CoolProp.CoolProp.PropsSI('S','P',P_c,'H',h7,refrigerant);
    h8 = py.CoolProp.CoolProp.PropsSI('H', 'P', P_c, 'Q', 0, refrigerant); s8 = py.CoolProp.CoolProp.PropsSI('S', 'P', P_c, 'Q', 0, refrigerant);
    h9 = h8; s9 = py.CoolProp.CoolProp.PropsSI('S', 'P', P_evap, 'H', h9, refrigerant);
    Q_spec_ej = h10 - h9; m_comp_ej = (Q_target / Q_spec_ej) / mu; 
    
    h1_ej = py.CoolProp.CoolProp.PropsSI('H', 'P', P_c, 'Q', 1, refrigerant); s1_ej = py.CoolProp.CoolProp.PropsSI('S', 'P', P_c, 'Q', 1, refrigerant);
    h2s_ej = py.CoolProp.CoolProp.PropsSI('H', 'P', P_cond, 'S', s1_ej, refrigerant);
    h2_ej = h1_ej + (h2s_ej - h1_ej) / eta_isen; s2_ej = py.CoolProp.CoolProp.PropsSI('S', 'P', P_cond, 'H', h2_ej, refrigerant);
    COP_ej(i) = Q_target / (m_comp_ej * (h2_ej - h1_ej));
    
    s1_ej_arr(i)=s1_ej; s2_ej_arr(i)=s2_ej; s3_ej_arr(i)=s3; s4_ej_arr(i)=s4_ej;
    s6_ej_arr(i)=s6; s7_ej_arr(i)=s7; s8_ej_arr(i)=s8; s9_ej_arr(i)=s9; s10_ej_arr(i)=s10;
end

%% Plotting
c_bg = 'w'; c_text = 'k'; line_w = 1.5;

% 1. COP vs T_evap
fig1 = figure('Color', c_bg, 'Name', 'COP vs T_evap'); ax1 = axes('Parent', fig1); hold(ax1, 'on');
plot(ax1, T_evap_array, COP_vcrs, '-o', 'LineWidth', line_w, 'DisplayName', 'Simple VCRS');
plot(ax1, T_evap_array, COP_ej, '--s', 'LineWidth', line_w, 'DisplayName', 'Ejector-based VCRS');
set(ax1, 'Color', c_bg, 'XColor', c_text, 'YColor', c_text, 'GridColor', c_text, 'GridAlpha', 0.15);
title(ax1, 'COP vs Evaporator Temperature', 'Color', c_text, 'FontWeight', 'bold');
xlabel(ax1, 'Evaporator Temperature (Kelvin)', 'Color', c_text, 'FontWeight', 'bold'); ylabel(ax1, 'COP', 'Color', c_text, 'FontWeight', 'bold');
grid(ax1, 'on'); legend(ax1, 'Location', 'best', 'TextColor', c_text, 'Color', 'w', 'EdgeColor', c_text); hold(ax1, 'off');

% 2. EER vs T_evap
fig2 = figure('Color', c_bg, 'Name', 'EER vs T_evap'); ax2 = axes('Parent', fig2); hold(ax2, 'on');
plot(ax2, T_evap_array, EER_vcrs, '-o', 'Color', [0.9290 0.6940 0.1250], 'LineWidth', line_w, 'DisplayName', 'Simple VCRS EER');
plot(ax2, T_evap_array, EER_ej, '--s', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', line_w, 'DisplayName', 'Ejector VCRS EER');
set(ax2, 'Color', c_bg, 'XColor', c_text, 'YColor', c_text, 'GridColor', c_text, 'GridAlpha', 0.15);
title(ax2, 'EER vs Evaporator Temperature', 'Color', c_text, 'FontWeight', 'bold');
xlabel(ax2, 'Evaporator Temperature (Kelvin)', 'Color', c_text, 'FontWeight', 'bold'); ylabel(ax2, 'EER', 'Color', c_text, 'FontWeight', 'bold');
grid(ax2, 'on'); legend(ax2, 'Location', 'best', 'TextColor', c_text, 'Color', 'w', 'EdgeColor', c_text); hold(ax2, 'off');

% 3. Compressor Work vs T_evap
fig3 = figure('Color', c_bg, 'Name', 'W_comp vs T_evap'); ax3 = axes('Parent', fig3); hold(ax3, 'on');
plot(ax3, T_evap_array, W_vcrs, '-o', 'LineWidth', line_w, 'DisplayName', 'Simple VCRS W_{comp}');
plot(ax3, T_evap_array, W_ej, '--s', 'LineWidth', line_w, 'DisplayName', 'Ejector VCRS W_{comp}');
set(ax3, 'Color', c_bg, 'XColor', c_text, 'YColor', c_text, 'GridColor', c_text, 'GridAlpha', 0.15);
title(ax3, 'Compressor Work vs Evaporator Temperature', 'Color', c_text, 'FontWeight', 'bold');
xlabel(ax3, 'Evaporator Temperature (K)', 'Color', c_text, 'FontWeight', 'bold'); ylabel(ax3, 'Compressor Work (kW)', 'Color', c_text, 'FontWeight', 'bold');
grid(ax3, 'on'); legend(ax3, 'Location', 'best', 'TextColor', c_text, 'Color', 'w', 'EdgeColor', c_text); hold(ax3, 'off');

% 4. Mass Flow Rates vs T_evap
fig4 = figure('Color', c_bg, 'Name', 'Mass Flows vs T_evap'); ax4 = axes('Parent', fig4); hold(ax4, 'on');
plot(ax4, T_evap_array, m_vcrs, '-o', 'LineWidth', line_w, 'DisplayName', 'VCRS m_{comp}');
plot(ax4, T_evap_array, m_total_ej_array, '-*', 'LineWidth', line_w, 'Color', [0.4660 0.6740 0.1880], 'DisplayName', 'Ejector m_{total}');
set(ax4, 'Color', c_bg, 'XColor', c_text, 'YColor', c_text, 'GridColor', c_text, 'GridAlpha', 0.15);
title(ax4, 'Mass Flow Rates vs Evaporator Temperature', 'Color', c_text, 'FontWeight', 'bold');
xlabel(ax4, 'Evaporator Temperature (Kelvin)', 'Color', c_text, 'FontWeight', 'bold'); ylabel(ax4, 'Mass Flow Rate (kg/s)', 'Color', c_text, 'FontWeight', 'bold');
grid(ax4, 'on'); legend(ax4, 'Location', 'best', 'TextColor', c_text, 'Color', 'w', 'EdgeColor', c_text); hold(ax4, 'off');

% 5. Entropy vs T_evap
figure('Color', 'w', 'Name', 'Entropy vs T_evap'); hold on;
plot(T_evap_array, s1_ej_arr, '-^', 'LineWidth', 1.5, 'DisplayName', 'Ejector Comp Inlet (s1)');
plot(T_evap_array, s2_ej_arr, '-v', 'LineWidth', 1.5, 'DisplayName', 'Ejector Comp Outlet (s2)');
plot(T_evap_array, s1_vcrs_arr, '--^', 'LineWidth', 1.5, 'DisplayName', 'VCRS Comp Inlet (s1)');
plot(T_evap_array, s2_vcrs_arr, '--v', 'LineWidth', 1.5, 'DisplayName', 'VCRS Comp Outlet (s2)');
title('Compressor Entropies vs Evaporator Temp');
xlabel('Evaporator Temp (K)'); ylabel('Entropy (J/kg-K)'); grid on; legend('Location', 'best'); hold off;

% 6. Pressure Lift Ratio (PLR) vs T_evap
fig5 = figure('Color', c_bg, 'Name', 'PLR vs T_evap'); ax5 = axes('Parent', fig5); hold(ax5, 'on');
plot(ax5, T_evap_array, PLR_array, '-^', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', line_w, 'DisplayName', 'Calculated PLR');
set(ax5, 'Color', c_bg, 'XColor', c_text, 'YColor', c_text, 'GridColor', c_text, 'GridAlpha', 0.15);
title(ax5, 'Calculated Pressure Lift Ratio vs Evaporator Temperature', 'Color', c_text, 'FontWeight', 'bold');
xlabel(ax5, 'Evaporator Temperature (Kelvin)', 'Color', c_text, 'FontWeight', 'bold'); ylabel(ax5, 'Pressure Lift Ratio (PLR)', 'Color', c_text, 'FontWeight', 'bold');
grid(ax5, 'on'); legend(ax5, 'Location', 'best', 'TextColor', c_text, 'Color', 'w', 'EdgeColor', c_text); hold(ax5, 'off');

% 7. Entrainment Ratio (omega) vs T_evap
fig6 = figure('Color', c_bg, 'Name', 'Entrainment Ratio vs T_evap'); ax6 = axes('Parent', fig6); hold(ax6, 'on');
plot(ax6, T_evap_array, omega_array, '-v', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', line_w, 'DisplayName', 'Calculated \mu');
set(ax6, 'Color', c_bg, 'XColor', c_text, 'YColor', c_text, 'GridColor', c_text, 'GridAlpha', 0.15);
title(ax6, 'Calculated Entrainment Ratio vs Evaporator Temperature', 'Color', c_text, 'FontWeight', 'bold');
xlabel(ax6, 'Evaporator Temperature (Kelvin)', 'Color', c_text, 'FontWeight', 'bold'); ylabel(ax6, 'Entrainment Ratio', 'Color', c_text, 'FontWeight', 'bold');
grid(ax6, 'on'); legend(ax6, 'Location', 'best', 'TextColor', c_text, 'Color', 'w', 'EdgeColor', c_text); hold(ax6, 'off');