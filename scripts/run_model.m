lambda = 100;
p = 1;
m = 50;
NPs = 2;
ta = 0.3;
betas = 0; % represent ART
EC50 = 0.2;


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
tw = 3;
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


% Initial conditions
if R0 > 1
    I_0t = max(w * I_art, I1 / p_surv);
else
    I_0t = I_art;
end

T_0 = T_art;
I_0 = I_0t;
V_0 = I_0 * (k / c);
P_0 = P_art;
E_0 = E_art;

y0 = [T_0 I_0 V_0 P_0 E_0, 1];


options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);
[T1,Y1] = ode45(@viral_load_model,[0 3],y0,options,lambda, NPs, p, m, ta, EC50);
