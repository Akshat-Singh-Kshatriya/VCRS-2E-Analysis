clc; clear; close all;

refrigerant = 'R744';
T_evap_array = [275.15, 277.15, 279.15, 281.15, 283.15];
T_cond = 308.15; 
Q_target = 3517;
eta_isen = 0.85;
eta_n = 0.85; 
P_gc = 110e5;
COP_paper = [3.4, 3.55, 3.75, 3.95, 4.15];

n = length(T_evap_array);
COP_code = zeros(1, n);
PLR_calc = zeros(1, n);
omega_calc = zeros(1, n);

h_gc_out = py.CoolProp.CoolProp.PropsSI('H', 'P', P_gc, 'T', T_cond, refrigerant);
s_gc_out = py.CoolProp.CoolProp.PropsSI('S', 'P', P_gc, 'T', T_cond, refrigerant);

for i = 1:n
    T_evap = T_evap_array(i);
    
    P_evap = py.CoolProp.CoolProp.PropsSI('P', 'T', T_evap, 'Q', 1, refrigerant);
    h_evap_out = py.CoolProp.CoolProp.PropsSI('H', 'T', T_evap, 'Q', 1, refrigerant); 
    
    h_n_s = py.CoolProp.CoolProp.PropsSI('H', 'P', P_evap, 'S', s_gc_out, refrigerant);
    h_n_out = h_gc_out - eta_n * (h_gc_out - h_n_s);
    u_n_out = sqrt(2 * (h_gc_out - h_n_out));
    
    mu_low = 0.01; 
    mu_high = 1.0;
    
    for iter = 1:100
        mu = (mu_low + mu_high) / 2;
        % Mixing chamber momentum and energy balance
        u_mix = u_n_out / (1 + mu);
        h_mix = ((h_n_out + 0.5 * u_n_out^2) + mu * h_evap_out) / (1 + mu) - 0.5 * u_mix^2;
        s_mix = py.CoolProp.CoolProp.PropsSI('S', 'P', P_evap, 'H', h_mix, refrigerant);
        
        % Diffuser energy balance
        h_diff_out = (h_gc_out + mu * h_evap_out) / (1 + mu);
        h_diff_s = h_mix + eta_n * (h_diff_out - h_mix);
        
        P_sep = py.CoolProp.CoolProp.PropsSI('P', 'S', s_mix, 'H', h_diff_s, refrigerant);
        
        
        try 
            x_sep = max(0, min(1, py.CoolProp.CoolProp.PropsSI('Q', 'P', P_sep, 'H', h_diff_out, refrigerant)));
        catch 
            if h_diff_out >= py.CoolProp.CoolProp.PropsSI('H', 'P', P_sep, 'Q', 1, refrigerant)
                x_sep = 1.0;
            else
                x_sep = 0.0;
            end
        end
           
        residual = (1 + mu) * x_sep - 1;
        if residual > 0 
            mu_high = mu;
        else 
            mu_low = mu;
        end
            
        if abs(residual) < 1e-4 
            break;
        end
    end
            

    omega = mu;
    omega_calc(i) = omega;
    PLR = P_sep / P_evap;
    PLR_calc(i) = PLR;
    
 
    h6 = py.CoolProp.CoolProp.PropsSI('H', 'P', P_sep, 'Q', 0, refrigerant);
    h7 = h6; 
    
    m_evap = Q_target / (h_evap_out - h7);
    m_comp = m_evap / omega;
    
    h1 = py.CoolProp.CoolProp.PropsSI('H', 'P', P_sep, 'Q', 1, refrigerant);
    s1 = py.CoolProp.CoolProp.PropsSI('S', 'P', P_sep, 'Q', 1, refrigerant);
    
    h2s = py.CoolProp.CoolProp.PropsSI('H', 'P', P_gc, 'S', s1, refrigerant);
    h2 = h1 + (h2s - h1) / eta_isen;
    
    W_total = m_comp * (h2 - h1);
    COP_code(i) = Q_target / W_total;
end


fprintf('%-12s | %-6s | %-6s | %-10s\n', 'Evaporator Temperature', 'COP', 'PLR', 'Omega (mu)');
fprintf('%s\n', repmat('-', 45, 1));
for i = 1:n
    fprintf('%-12.2f | %-6.3f | %-6.3f | %-10.3f\n', T_evap_array(i), COP_code(i), PLR_calc(i), omega_calc(i));
end

%% PLOTTING 
fixed_param = sprintf('Gas Cooler Exit Temp = %.2f K', T_cond);

create_styled_plot(T_evap_array, COP_paper, COP_code, ...
                   'Evaporator Temperature (K)', 'COP', ...
                   'Research Paper Values', 'Code Generated Values', ...
                   fixed_param, '#FF0000', '#0055FF', 'COP_Validation.png');


function create_styled_plot(x, y1, y2, x_lbl, y_lbl, leg1, leg2, fixed_text, c1, c2, filename)
    fig = figure('Position', [100, 100, 800, 600], 'Color', 'w');
    hold on;
    
    
    ax = gca;
    ax.Color = 'w';
    ax.XColor = 'k'; % Force X-axis and tick labels to black
    ax.YColor = 'k'; % Force Y-axis and tick labels to black
    
    % Grid settings
    grid on;
    ax.GridColor = '#D3D3D3';
    ax.GridLineStyle = '-';
    ax.GridAlpha = 1.0; 
    
    % Font and axis styling
    ax.FontName = 'Times New Roman';
    ax.FontSize = 14;
    ax.LineWidth = 1.2;

    % Plot lines
    plot(x, y1, '-o', 'Color', c1, 'DisplayName', leg1, ...
         'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'w', ...
         'MarkerEdgeColor', c1);
     
    plot(x, y2, '-s', 'Color', c2, 'DisplayName', leg2, ...
         'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'w', ...
         'MarkerEdgeColor', c2);

    % Force axis labels to black
    xlabel(x_lbl, 'FontWeight', 'bold', 'FontSize', 16, 'Color', 'k');
    ylabel(y_lbl, 'FontWeight', 'bold', 'FontSize', 16, 'Color', 'k');
    
    % Add padding to x-axis
    x_padding = (max(x) - min(x)) * 0.05;
    xlim([min(x) - x_padding, max(x) + x_padding]);

    % Add dynamic Y-axis padding
    y_min = min(min(y1), min(y2));
    y_max = max(max(y1), max(y2));
    y_range = y_max - y_min;
    if y_range == 0
        y_range = 1; 
    end
    ylim([y_min - y_range * 0.05, y_max + y_range * 0.35]);

    % Legend & Text Box
    leg = legend('Location', 'best', 'EdgeColor', 'k');
    t = title(leg, fixed_text, 'FontWeight', 'normal', 'FontSize', 13);
    t.Color = 'k'; 
    leg.FontSize = 13;
    leg.LineWidth = 1.0;
    leg.Color = 'w';   
    leg.TextColor = 'k';  

    hold off;
   
    exportgraphics(fig, filename, 'Resolution', 600, 'BackgroundColor', 'w');
end