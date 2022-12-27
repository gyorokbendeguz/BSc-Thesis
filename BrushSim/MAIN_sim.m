close all; 
clear;
clc;

%% Felhasználói beállítások / User inputs

calc_fric = 1;      % súrlódási együttható számítás / calculating friction coeff.
anim_plot = 1;      % animáció kirajzolása / showing the animation
anim_SolCalc = 1;   % analitikus közelítés a számított + pontos mu-vel animálva/ 
                    % animation of the analytical approx. w/ calculated + accurate mu

%% Szimuláció bemeneti adatainak megadása / Initializing the simulation

global l k J  b x dx n Fn a mu qcrp eps Mdes qdiff psi_max

l = 0.01;
a = 0.04;
k = 60000;
J = 0.1;
b = 2*sqrt(2*a*k*J*(l^2 + (a^3)/3));
Fn = 180;
mu = 1;
n = 200;
eps = 4e-10;

Mdes = [0 0;
        0.01 0.4*2;
        2 0.8*2;
        4 -0.9*2;
        5 0.2*2;
        5.1 0.2*2];

dx = 2*a/n;
x = linspace(-a,a,n+1);

qcrp = 3*Fn*mu/(4*a^3*k)*(a^2-x.^2);
qdiff = -3*Fn*mu/(2*a^3*k);

% kezdeti ertekek / initial conditions
y0 = [0;0;((1:n+1)==0)'];

% lepeskoz / timestep size
dt = 0.01;%dx/v/2

% idolepesek szama / number of timestep
pont = 2^9;

% kezdeti idopillanat / initial time value
t0 = 0;

% vegso idopillanat / final time value
%%pont = fix((tf-t0)/dt) + 2;
t = t0 + (0:pont-1)*dt;
tf = t0 + pont*dt;

%% Szimuláció / Simulation

%%%%%%%%%%%%%%%%%%%%%%
% SIMULATION 
%%%%%%%%%%%%%%%%%%%%%%
y = zeros(pont, n+3);
y(1, :) = y0;
tic;
ii = 1;
wb = waitbar(0, 'Simulation in process...');
M = zeros(pont, 1);
psi_max = 0;
Psi_max = zeros(pont, 1);
max_index = 1;
q_anal = zeros(pont, n+1);
q_anal(1, :) = ((1:n+1)==0)';

%analitikus megoldáshoz állapotautomata / state machine for analytical sol.
state = "3P"; %kezdeti állapot / initial state
States = string(zeros(1, pont));
States(1, 1) = state;

for i=2:pont
    % numerikus szimuláció / numerical simulation
    szim_index = i;
    [y(i,:), y(i-1,:), M(i, 1)] = solver(t(i),y(i-1,:),dt);
    
    %analitikus / analytical
    States(1, i) = state;
    [q_anal(i,:), state] = anal_sol(y(i, :), y(i-1, :), mu, state);
    Psi_max(i, 1) = psi_max;
    
    if i > (fix(pont/10)*ii + 1)
        ii = ii + 1;
    end
    waitbar(i/pont,wb)
end
simtime = toc
close(wb);

%% Plots / Grafikonok

% elfordulas az idoben
figure(1);
M_le = interp1(Mdes(:,1), Mdes(:,2), t, 'linear', 'extrap');
subplot(3, 1, 1)
hold on;
plot(t, M_le, 'k-', 'LineWidth', 1);
hold off;
xlim([0, t(end)]);
xlabel('Time [s]');
ylabel('Steering torque [Nm]');
box on;
grid on;

subplot(3, 1, 2)
hold on;
plot(t, y(:,1), 'k-', 'LineWidt', 1);
hold off;
xlim([0, t(end)]);
xlabel('Time [s]');
ylabel('\psi [rad]');
grid on;
box on;

subplot(3, 1, 3)
plot(t, y(:,2), 'k-', 'LineWidt', 1);
xlim([0, t(end)]);
xlabel('Time [s]');
ylabel("\psi' [rad/sec]");
grid on;

%% Súrlódási eh. számítása numerikusan / Calculating the friction coeff. numerically

if calc_fric == true
    tic;
    mu_calc = fric_calc_v2(pont, y(:,1), Psi_max, M, States);
    mu_real = ones(1, pont).*mu;
    calctime = toc
end

%% Animáció / Animation

if anim_plot == true
    lim_y = 3*Fn*mu/(4*a*k);
    if calc_fric == false
        fig_sim = figure;
        for i=1:pont
                plot(x(:), qcrp, 'k--');
                hold on;
                plot(x(:), -qcrp, 'k--');
                plot(x(:), y(i,3:n+3), 'b.-');
                plot(x(:), q_anal(i,:), 'r-');
                hold off;
                set(fig_sim,'Color','w')
                grid on;
                ylim([-lim_y, lim_y]);
                legend(["", "", "Numerikus", "Analitikus"]);
                xlabel('x [m]');
                ylabel('q(x,t) [m]');
                drawnow;
        end
    else
        for i=1:pont
            figure(3);
            subplot(2,2,[1,3]);
            plot(x(:), y(i,3:n+3), '.-');
            hold on;
            plot(x(:), q_anal(i,:), 'r-');
            hold off;
            grid on;
            ylim([-lim_y, lim_y]);
            legend(["Numerikus", "Analitikus"]);
                
            subplot(2,2,2);
            plot(t(1:i), mu_real(1:i), 'k');
            hold on;
            plot(t(1:i), mu_calc(1:i), 'r');
            hold off;
            xlabel("Time [s]");
            ylabel("Coeff. of fric. [1]");
            legend("Actual", "Estimated");
            xlim([t(1), t(end)]);
            ylim([0, 2*mu]);
            grid on;
            
            [MuPlot, FmuPlot] = anim_mu(t(i), M(i), y(i,1), Psi_max(i), States(i));
            RightSideOfEq = zeros(size(MuPlot));
            subplot(2,2,4);
            plot(MuPlot, FmuPlot, 'b');
            hold on;
            plot(MuPlot, RightSideOfEq, 'r');
            hold off;
            grid on;
            xlabel("Coeff. of fric [1]");
            ylabel("f(mu) [1]");
            ylim([-4,4]);
            xlim([-2*mu, 2*mu]);
            
            drawnow;
        end
    end
end

if (anim_SolCalc == true) && (calc_fric == true)
    q_calc = zeros(pont, n+1);
    q_calc(1, :) = ((1:n+1)==0)';
    state_calc = "3P";
    for i = 2:pont
       [q_calc(i,:), state_calc] = anal_sol(y(i, :), y(i-1, :), mu_calc(i), state_calc); 
    end
    lim_y = 3*Fn*mu/(4*a*k);
    
    for i=1:pont
        figure(5);
        plot(x(:), q_anal(i,:), '-b');
        hold on;
        plot(x(:), q_calc(i,:), 'r');
        hold off;
        grid on;
        ylim([-lim_y, lim_y]);
        legend(["Accurate", "Calculated"]);
        drawnow;
    end
end

if calc_fric == true
    figure;
    plot(t, mu_real, 'k', 'LineWidth', 1);
    hold on;
    plot(t, 1.1.*mu_real, '--k');
    plot(t, 0.9.*mu_real, '--k');
    plot(t, mu_calc, 'r', 'LineWidth', 1);
    hold off;
    xlabel("Time [s]");
    ylabel("Friction coefficient [1]");
    legend("Actual", "10% error", "", "Estimated");
    xlim([t(1), t(end)]);
end
