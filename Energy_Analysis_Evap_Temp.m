clc; clear; close all;

refrigerants = {'R1234yf', 'R1234ze(E)'}; 
T_evap_array = [275.15, 277.15, 279.15, 281.15, 283.15];
T_cond = 308.15;
Q_target = 3517;
eta_n = 0.85;
eta_isen = 0.85;

n = length(T_evap_array);
results_store = struct();

for r = 1:length(refrigerants)
    ref = refrigerants{r};
    ref_field = matlab.lang.makeValidName(ref);

    P_cond = py.CoolProp.CoolProp.PropsSI('P', 'T', T_cond, 'Q', 0, ref);
    
    COP_vcrs = zeros(1, n); W_vcrs = zeros(1, n); EER_vcrs = zeros(1, n); m_vcrs = zeros(1, n);
    COP_ej = zeros(1, n); W_ej = zeros(1, n); EER_ej = zeros(1, n); m_total_ej = zeros(1, n);
    P_evap_arr = zeros(1, n); 
    
    % --- NEW: Initialize arrays for Condenser Heat Rejection ---
    Q_cond_vcrs = zeros(1, n); 
    Q_cond_ej = zeros(1, n);
    
    for i = 1:n
        T_evap = T_evap_array(i);

        P_evap = py.CoolProp.CoolProp.PropsSI('P', 'T', T_evap, 'Q', 1, ref);
        P_evap_arr(i) = P_evap; 
        
        h1_vcrs = py.CoolProp.CoolProp.PropsSI('H', 'P', P_evap, 'Q', 1, ref);
        s1_vcrs = py.CoolProp.CoolProp.PropsSI('S', 'P', P_evap, 'Q', 1, ref);
        h2s_vcrs = py.CoolProp.CoolProp.PropsSI('H', 'P', P_cond, 'S', s1_vcrs, ref);
        h2_vcrs = h1_vcrs + (h2s_vcrs - h1_vcrs) / eta_isen;
        
        h3 = py.CoolProp.CoolProp.PropsSI('H', 'P', P_cond, 'Q', 0, ref);
        s3 = py.CoolProp.CoolProp.PropsSI('S', 'P', P_cond, 'Q', 0, ref);
        h4_vcrs = h3;
        
        Q_spec_vcrs = h1_vcrs - h4_vcrs;
        W_spec_vcrs = h2_vcrs - h1_vcrs;
        
        COP_vcrs(i) = Q_spec_vcrs / W_spec_vcrs;
        EER_vcrs(i) = COP_vcrs(i) * 3.41214;
        m_vcrs(i) = Q_target / Q_spec_vcrs;
        W_vcrs(i) = m_vcrs(i) * W_spec_vcrs / 1000;
        
        % --- NEW: Calculate VCRS Condenser Heat Rejection (kW) ---
        Q_cond_vcrs(i) = m_vcrs(i) * (h2_vcrs - h3) / 1000; 
        
        h4s = py.CoolProp.CoolProp.PropsSI('H', 'P', P_evap, 'S', s3, ref);
        h4 = h3 - eta_n * (h3 - h4s);
        u4 = sqrt(2 * (h3 - h4));
        h10 = py.CoolProp.CoolProp.PropsSI('H', 'T', T_evap, 'Q', 1, ref);
        
        mu_low = 0.01; 
        mu_high = 1.0;
        
        for iter = 1:100
            mu = (mu_low + mu_high) / 2;
            u6 = u4 / (1 + mu);
            h6 = ((h4 + 0.5 * u4^2) + mu * h10) / (1 + mu) - 0.5 * u6^2;
            s6 = py.CoolProp.CoolProp.PropsSI('S', 'P', P_evap, 'H', h6, ref);
            h7 = (h3 + mu * h10) / (1 + mu);
            h7s = h6 + eta_n * (h7 - h6);
            P_c = py.CoolProp.CoolProp.PropsSI('P', 'S', s6, 'H', h7s, ref);
            
            try
                x7 = py.CoolProp.CoolProp.PropsSI('Q', 'P', P_c, 'H', h7, ref);
                x7 = max(0.0, min(1.0, x7));
            catch
                H_sat_vap = py.CoolProp.CoolProp.PropsSI('H', 'P', P_c, 'Q', 1, ref);
                if h7 >= H_sat_vap
                    x7 = 1.0;
                else
                    x7 = 0.0;
                end
            end
                
            residual = (1 + mu) * x7 - 1;
            if residual > 0
                mu_high = mu;
            else
                mu_low = mu;
            end
            if abs(residual) < 1e-4
                break;
            end
        end
        
        h8 = py.CoolProp.CoolProp.PropsSI('H', 'P', P_c, 'Q', 0, ref);
        h9 = h8;
        
        Q_spec_ej = h10 - h9;
        m_comp_ej = (Q_target / Q_spec_ej) / mu;
        
        h1_ej = py.CoolProp.CoolProp.PropsSI('H', 'P', P_c, 'Q', 1, ref);
        s1_ej = py.CoolProp.CoolProp.PropsSI('S', 'P', P_c, 'Q', 1, ref);
        h2s_ej = py.CoolProp.CoolProp.PropsSI('H', 'P', P_cond, 'S', s1_ej, ref);
        h2_ej = h1_ej + (h2s_ej - h1_ej) / eta_isen;
        
        COP_ej(i) = Q_target / (m_comp_ej * (h2_ej - h1_ej));
        EER_ej(i) = COP_ej(i) * 3.41214;
        W_ej(i) = m_comp_ej * (h2_ej - h1_ej) / 1000; 
        m_total_ej(i) = m_comp_ej * (1 + mu);
        
        % --- NEW: Calculate Ejector Condenser Heat Rejection (kW) ---
        Q_cond_ej(i) = m_comp_ej * (h2_ej - h3) / 1000; 
    end
    
    results_store.(ref_field).COP_vcrs = COP_vcrs;
    results_store.(ref_field).P_evap = P_evap_arr; 
    results_store.(ref_field).COP_ej = COP_ej;
    results_store.(ref_field).EER_vcrs = EER_vcrs;
    results_store.(ref_field).EER_ej = EER_ej;
    results_store.(ref_field).W_vcrs = W_vcrs;
    results_store.(ref_field).W_ej = W_ej;
    results_store.(ref_field).m_vcrs = m_vcrs;
    results_store.(ref_field).m_total_ej = m_total_ej;
    
    % --- NEW: Store Condenser values ---
    results_store.(ref_field).Q_cond_vcrs = Q_cond_vcrs;
    results_store.(ref_field).Q_cond_ej = Q_cond_ej;
end

fprintf('                     EJECTOR VCRS vs SIMPLE VCRS PERFORMANCE COMPARISON                      \n');

for r = 1:length(refrigerants)
    ref = refrigerants{r};
    ref_field = matlab.lang.makeValidName(ref); 
    
    fprintf('\n--- %s Results ---\n', ref);
    
    % --- NEW: Expanded header for the new columns ---
    fprintf('%-13s | %-15s | %-11s | %-11s | %-11s | %-11s | %-12s | %-12s\n', ...
            'Evap Temp (K)', 'Evap Pres (kPa)', 'Ejector COP', 'COP Imp (%)', 'VCRS W (kW)', 'Ej W (kW)', 'VCRS Qc (kW)', 'Ej Qc (kW)');
    fprintf('%s\n', repmat('-', 1, 115));
    
    for i = 1:n
        cop_vcrs_val = results_store.(ref_field).COP_vcrs(i);
        cop_ej_val   = results_store.(ref_field).COP_ej(i);
        p_evap_val   = results_store.(ref_field).P_evap(i);
        w_vcrs_val   = results_store.(ref_field).W_vcrs(i);
        w_ej_val     = results_store.(ref_field).W_ej(i);
        
        % --- NEW: Extract Condenser values ---
        q_cond_vcrs_val = results_store.(ref_field).Q_cond_vcrs(i);
        q_cond_ej_val   = results_store.(ref_field).Q_cond_ej(i);
        
        cop_imp = ((cop_ej_val - cop_vcrs_val) / cop_vcrs_val) * 100;
        
        fprintf('%-13.2f | %-15.2f | %-11.3f | %-11.2f | %-11.3f | %-11.3f | %-12.3f | %-12.3f\n', ...
            T_evap_array(i), ...
            (p_evap_val / 1000), ... 
            cop_ej_val, ...
            cop_imp, ...
            w_vcrs_val, ...
            w_ej_val, ...
            q_cond_vcrs_val, ...
            q_cond_ej_val);
    end
end
ref1 = matlab.lang.makeValidName(refrigerants{1}); 
ref2 = matlab.lang.makeValidName(refrigerants{2}); 
fixed_param = sprintf('T_{cond} = %.2f K', T_cond);

c_ref1_vcrs = '#0066FF'; 
c_ref1_ej   = '#FF0000'; 
c_ref2_vcrs = '#00AA00'; 
c_ref2_ej   = '#800080'; 

create_combined_plot(T_evap_array, ...
    results_store.(ref1).COP_vcrs, results_store.(ref1).COP_ej, ...
    results_store.(ref2).COP_vcrs, results_store.(ref2).COP_ej, ...
    'Evaporator Temperature (K)', 'COP', ...
    'R1234yf Simple VCRS', 'R1234yf Ejector VCRS', ...
    'R1234ze(E) Simple VCRS', 'R1234ze(E) Ejector VCRS', ...
    fixed_param, c_ref1_vcrs, c_ref1_ej, c_ref2_vcrs, c_ref2_ej, 'Combined_COP.png');

create_combined_plot(T_evap_array, ...
    results_store.(ref1).EER_vcrs, results_store.(ref1).EER_ej, ...
    results_store.(ref2).EER_vcrs, results_store.(ref2).EER_ej, ...
    'Evaporator Temperature (K)', 'EER', ...
    'R1234yf Simple VCRS', 'R1234yf Ejector VCRS', ...
    'R1234ze(E) Simple VCRS', 'R1234ze(E) Ejector VCRS', ...
    fixed_param, c_ref1_vcrs, c_ref1_ej, c_ref2_vcrs, c_ref2_ej, 'Combined_EER.png');

create_combined_plot(T_evap_array, ...
    results_store.(ref1).W_vcrs, results_store.(ref1).W_ej, ...
    results_store.(ref2).W_vcrs, results_store.(ref2).W_ej, ...
    'Evaporator Temperature (K)', 'Compressor Work (kW)', ...
    'R1234yf Simple VCRS', 'R1234yf Ejector VCRS', ...
    'R1234ze(E) Simple VCRS', 'R1234ze(E) Ejector VCRS', ...
    fixed_param, c_ref1_vcrs, c_ref1_ej, c_ref2_vcrs, c_ref2_ej, 'Combined_Work.png');

create_combined_plot(T_evap_array, ...
    results_store.(ref1).m_vcrs, results_store.(ref1).m_total_ej, ...
    results_store.(ref2).m_vcrs, results_store.(ref2).m_total_ej, ...
    'Evaporator Temperature (K)', 'Mass Flow Rate (kg/s)', ...
    'R1234yf VCRS m_{comp}', 'R1234yf Ejector m_{total}', ...
    'R1234ze(E) VCRS m_{comp}', 'R1234ze(E) Ejector m_{total}', ...
    fixed_param, c_ref1_vcrs, c_ref1_ej, c_ref2_vcrs, c_ref2_ej, 'Combined_MassFlow.png');

function create_combined_plot(x, y1, y2, y3, y4, x_lbl, y_lbl, leg1, leg2, leg3, leg4, fixed_text, c1, c2, c3, c4, filename)
    fig = figure('Position', [100, 100, 800, 600], 'Color', 'w', 'Visible', 'off'); 
    hold on;
    
    ax = gca;
    ax.Color = 'w';
    ax.XColor = 'k'; 
    ax.YColor = 'k'; 
    
    grid on;
    ax.GridColor = '#D3D3D3';
    ax.GridLineStyle = '-';
    ax.GridAlpha = 1.0; 
    
    ax.FontName = 'Times New Roman';
    ax.FontSize = 14;
    ax.LineWidth = 1.2;

    plot(x, y1, '-o', 'Color', c1, 'DisplayName', leg1, 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', c1);
    plot(x, y2, '-s', 'Color', c2, 'DisplayName', leg2, 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', c2);
    plot(x, y3, '-^', 'Color', c3, 'DisplayName', leg3, 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', c3);
    plot(x, y4, '-d', 'Color', c4, 'DisplayName', leg4, 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', c4);

    xlabel(x_lbl, 'FontWeight', 'bold', 'FontSize', 16, 'Color', 'k');
    ylabel(y_lbl, 'FontWeight', 'bold', 'FontSize', 16, 'Color', 'k');
    
    x_padding = (max(x) - min(x)) * 0.05;
    xlim([min(x) - x_padding, max(x) + x_padding]);

    y_min = min([min(y1), min(y2), min(y3), min(y4)]);
    y_max = max([max(y1), max(y2), max(y3), max(y4)]);
    y_range = y_max - y_min;
    if y_range == 0
        y_range = 1; 
    end
    ylim([y_min - y_range * 0.05, y_max + y_range * 0.35]);

    leg = legend('Location', 'best', 'EdgeColor', 'k', 'NumColumns', 2);
    t = title(leg, fixed_text, 'FontWeight', 'normal', 'FontSize', 13);
    t.Color = 'k'; 
    leg.FontSize = 11;
    leg.LineWidth = 1.0;
    leg.Color = 'w';      
    leg.TextColor = 'k';  

    hold off;
    
    exportgraphics(fig, filename, 'Resolution', 600, 'BackgroundColor', 'w');
    close(fig); 
end