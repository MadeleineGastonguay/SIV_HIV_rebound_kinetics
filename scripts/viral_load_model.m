function dydt = viral_load_model(t, y, lambda, NPs, p, m, ta, EC50)
    % Use effect model to describe ta as function of drug:
    % ta = 0.3 + TLR7/(0.2 + TLR7);
    betas = 0; % represent ART


    % Fixed parameters
    NP = NPs * 10^5;
    beta = betas * 10^(-7);
    dT = 0.05;
    dI = 0.4;
    dP = 0.001;
    dE = 1;
    c = 23;
    k = 50000;
    NE = 10000;
    F = 0.9;
    tw = 0;
    I1 = 4 * 10^(-5);
    rho = 10;
    N = 3;
    
    % Variable values during ART
    T_art = lambda / dT;
    I_art = I1 / (ta * dI);
    V_art = I_art * k / c;
    
    if dP > p * (1 - F) * V_art / (V_art + NP)
        P_art = m / (dP - p * (1 - F) * V_art / (V_art + NP));
    else
        P_art = 1000 * NE;
    end
    
    E_art = p * F * P_art * V_art / (dE * (V_art + NP));
    
    % Composite parameters
    R0 = lambda * beta * k / ((1 + E_art / NE) * dT * dI * c);
    R = (-(c + dI) + sqrt((c - dI)^2 + 4 * c * dI * R0)) / 2;
    
    % Solve for survival probability using approximation to Lambert W
    if R0 > 1
        x_lamW = -R0 * exp(-R0);
        N2 = 3 * sqrt(2) + 6 - sqrt(2 * (1 + exp(1) * x_lamW)) * ((2237 + 1457 * sqrt(2)) * exp(1) - 4108 * sqrt(2) - 5764) / ((215 + 199 * sqrt(2)) * exp(1) - 430 * sqrt(2) - 796);
        N1 = (1 - 1 / sqrt(2)) * (N2 + sqrt(2));
        LamW_x = -1 + sqrt(2 * (1 + exp(1) * x_lamW)) / (1 + (N1 * sqrt(2 * (1 + exp(1) * x_lamW))) / (N2 + sqrt(2 * (1 + exp(1) * x_lamW))));
        p_surv = (LamW_x + R0) / (R0 + rho - 1);
    else
        p_surv = 0;
    end
    
    % Get effective a value
    if R0 > 1
        aeff = I1 * R / (p_surv * (exp(R * ta / p_surv) - 1));
    else
        aeff = I1 / ta;
    end
    
    % Effective weight of non-zero initial condition
    if R0 > 1
        w = 2^(-(ta * dI / p_surv)^N);
    else
        w = 0;
    end
    
    % Time to start numerical integration
    if R0 > 1
        if ta > 200
            t_go = 200 / p_surv - log(w * p_surv * (exp(R * 200 / p_surv) - 1) / (200 * dI) + 1) / R;
        else
            t_go = max(0, ta / p_surv - log(w * p_surv * (exp(R * ta / p_surv) - 1) / (ta * dI) + 1) / R);
        end
    else
        t_go = 0;
    end
    
     
    % Define condition
    t = 0; % Initiation of therapy
    if t < tw + t_go
        cond = 0;
    else
        cond = 1;
    end

    dydt = zeros(5,1);
    
    % Differential equations
    dydt(1) = cond * (lambda - beta * y(1) * y(3) - dT * y(1));
    dydt(2) = cond * ((aeff*y(6)/(EC50 + y(6))) + beta * y(1) * y(3) / (1 + y(5) / NE) - dI * y(2));
    dydt(3) = cond * (k * y(2) - c * y(3));
    dydt(4) = cond * (m + p * (1 - F) * y(3) * y(4) / (y(3) + NP) - dP * y(4));
    dydt(5) = cond * (p * F * y(3) * y(4) / (y(3) + NP) - dE * y(5));
    dydt(6) = -log(2)/10*24*y(6);
    % dydt(6) = 0;
    
    
end
