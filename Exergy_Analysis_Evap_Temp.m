clc; clear; close all;

refrigerants = {'R1234yf', 'R1234ze(E)'}; 
T_evap_array = [275.15, 277.15, 279.15, 281.15, 283.15];
T_cond = 308.15;
T0 = 298.15; 
Q_target = 3517;
eta_n = 0.85;
eta_isen = 0.85;

n = length(T_evap_array);
results_store = struct();

for r = 1:length(refrigerants)
    ref = refrigerants{r};
    ref_field = matlab.lang.makeValidName(ref);

    P_cond = py.CoolProp.CoolProp.PropsSI('P', 'T', T_cond, 'Q', 0, ref);
    
    eta_ex_vcrs = zeros(1, n); EDR_vcrs = zeros(1, n);
    eta_ex_ej = zeros(1, n); EDR_ej = zeros(1, n);
    
    for i = 1:n
        T_evap = T_evap_array(i);
        P_evap = py.CoolProp.CoolProp.PropsSI('P', 'T', T_evap, 'Q', 1, ref);
        
        h1_vcrs = py.CoolProp.CoolProp.PropsSI('H', 'P', P_evap, 'Q', 1, ref);
        s1_vcrs = py.CoolProp.CoolProp.PropsSI('S', 'P', P_evap, 'Q', 1, ref);
        h2s_vcrs = py.CoolProp.CoolProp.PropsSI('H', 'P', P_cond, 'S', s1_vcrs, ref);
        h2_vcrs = h1_vcrs + (h2s_vcrs - h1_vcrs) / eta_isen;
        
        h3 = py.CoolProp.CoolProp.PropsSI('H', 'P', P_cond, 'Q', 0, ref);
        s3 = py.CoolProp.CoolProp.PropsSI('S', 'P', P_cond, 'Q', 0, ref);
        h4_vcrs = h3;
        
        Q_spec_vcrs = h1_vcrs - h4_vcrs;
        W_spec_vcrs = h2_vcrs - h1_vcrs;
        
        COP_vcrs = Q_spec_vcrs / W_spec_vcrs;
        eta_ex_vcrs(i) = COP_vcrs * (T0 / T_evap - 1);
        EDR_vcrs(i) = (1 / eta_ex_vcrs(i)) - 1;
        
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
        
        COP_ej = Q_target / (m_comp_ej * (h2_ej - h1_ej));
        eta_ex_ej(i) = COP_ej * (T0 / T_evap - 1);
        EDR_ej(i) = (1 / eta_ex_ej(i)) - 1;
    end
    
    results_store.(ref_field).eta_ex_vcrs = eta_ex_vcrs;
    results_store.(ref_field).eta_ex_ej = eta_ex_ej;
    results_store.(ref_field).EDR_vcrs = EDR_vcrs;
    results_store.(ref_field).EDR_ej = EDR_ej;
end


for r = 1:length(refrigerants)
    ref = refrigerants{r};
    ref_field = matlab.lang.makeValidName(ref); 
    
    fprintf('\n%s Results\n', ref);
    fprintf('%-10s | %-12s | %-12s | %-12s | %-12s\n', 'Evap Temp (K)', 'VCRS eta_ex', 'Ej eta_ex', 'VCRS EDR', 'Ej EDR');
    fprintf('%s\n', repmat('-', 1, 65));
    
    for i = 1:n
        fprintf('%-10.2f | %-12.4f | %-12.4f | %-12.4f | %-12.4f\n', ...
            T_evap_array(i), ...
            results_store.(ref_field).eta_ex_vcrs(i), ...
            results_store.(ref_field).eta_ex_ej(i), ...
            results_store.(ref_field).EDR_vcrs(i), ...
            results_store.(ref_field).EDR_ej(i));
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
    results_store.(ref1).eta_ex_vcrs, results_store.(ref1).eta_ex_ej, ...
    results_store.(ref2).eta_ex_vcrs, results_store.(ref2).eta_ex_ej, ...
    'Evaporator Temperature (K)', 'Exergy Efficiency (\eta_{ex})', ...
    'R1234yf Simple VCRS', 'R1234yf Ejector VCRS', ...
    'R1234ze(E) Simple VCRS', 'R1234ze(E) Ejector VCRS', ...
    fixed_param, c_ref1_vcrs, c_ref1_ej, c_ref2_vcrs, c_ref2_ej, 'Combined_Evap_ExergyEff.png');

create_combined_plot(T_evap_array, ...
    results_store.(ref1).EDR_vcrs, results_store.(ref1).EDR_ej, ...
    results_store.(ref2).EDR_vcrs, results_store.(ref2).EDR_ej, ...
    'Evaporator Temperature (K)', 'Exergy Destruction Ratio (EDR)', ...
    'R1234yf Simple VCRS', 'R1234yf Ejector VCRS', ...
    'R1234ze(E) Simple VCRS', 'R1234ze(E) Ejector VCRS', ...
    fixed_param, c_ref1_vcrs, c_ref1_ej, c_ref2_vcrs, c_ref2_ej, 'Combined_Evap_EDR.png');

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