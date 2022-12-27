function [q, state_end] = anal_sol(y, y_prev, mu, STATE)%, xLt, xRt)

global x l k Fn a n qcrp psi_max

psi = y(1);
epsilon = y(2);
psi_prev = y_prev(1);
epsilon_prev = y_prev(2);
q = zeros(1, n+1);

switch STATE %---------------------------------------------------------
    case "3P" %3 szakasz, pozitív irányba
        xL = (0.3333*a*((9*Fn^2*mu^2 + sign(psi)*12*l*psi*Fn*a*k*mu + 4*psi^2*a^4*k^2)^(1/2) - sign(psi)*2*a^2*k*psi))/(Fn*mu);
        xR = -(0.3333*a*((9*Fn^2*mu^2 - sign(psi)*12*l*psi*Fn*a*k*mu + 4*psi^2*a^4*k^2)^(1/2) - sign(psi)*2*a^2*k*psi))/(Fn*mu);
        
        for i = 1:n+1
            if x(i) < xR %hátsó csúszó szakasz
                q(i) = sign(psi)*qcrp(i);
            elseif x(i) > xL %első csúszó szakasz
                q(i) = -sign(psi)*qcrp(i);
            else %tapadás
                q(i) = -psi*(x(i) - l);
            end
        end
        
    case "3M" %3 szakasz, negatív irányba
        xL = (0.3333*a*((9*Fn^2*mu^2 + sign(psi)*12*l*psi*Fn*a*k*mu + 4*psi^2*a^4*k^2)^(1/2) - sign(psi)*2*a^2*k*psi))/(Fn*mu);
        xR = -(0.3333*a*((9*Fn^2*mu^2 - sign(psi)*12*l*psi*Fn*a*k*mu + 4*psi^2*a^4*k^2)^(1/2) - sign(psi)*2*a^2*k*psi))/(Fn*mu);
        for i = 1:n+1
            if x(i) < xR %hátsó csúszó szakasz
                q(i) = sign(psi)*qcrp(i);
            elseif x(i) > xL %első csúszó szakasz
                q(i) = -sign(psi)*qcrp(i);
            else %tapadás
                q(i) = -psi*(x(i) - l);
            end
        end
        
    case "5M"
        mL = -(l*psi - (0.3333*a*psi*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 + 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu) + (0.7500*Fn*mu*(a^2 - (0.1111*a^2*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 + 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max))^2)/(Fn^2*mu^2)))/(a^3*k))/(l - (0.3333*a*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 + 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu));
        mR = -(l*psi + (0.3333*a*psi*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 - 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu) - (0.7500*Fn*mu*(a^2 - (0.1111*a^2*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 - 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max))^2)/(Fn^2*mu^2)))/(a^3*k))/(l + (0.3333*a*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 - 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu));
        xL = (0.3333*a^3*k*(((l*psi - (0.3333*a*psi*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 + 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu) + (0.7500*Fn*mu*(a^2 - (0.1111*a^2*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 + 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max))^2)/(Fn^2*mu^2)))/(a^3*k))^2/(l - (0.3333*a*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 + 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu))^2 + (6*Fn*mu*((1.5000*Fn*mu)/(a*k) - (l*(l*psi - (0.3333*a*psi*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 + 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu) + (0.7500*Fn*mu*(a^2 - (0.1111*a^2*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 + 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max))^2)/(Fn^2*mu^2)))/(a^3*k)))/(l - (0.3333*a*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 + 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu))))/(a^3*k))^(1/2) + (l*psi - (0.3333*a*psi*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 + 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu) + (0.7500*Fn*mu*(a^2 - (0.1111*a^2*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 + 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max))^2)/(Fn^2*mu^2)))/(a^3*k))/(l - (0.3333*a*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 + 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu))))/(Fn*mu);
        xR = -(0.3333*a^3*k*(((l*psi + (0.3333*a*psi*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 - 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu) - (0.7500*Fn*mu*(a^2 - (0.1111*a^2*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 - 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max))^2)/(Fn^2*mu^2)))/(a^3*k))^2/(l + (0.3333*a*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 - 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu))^2 + (6*Fn*mu*((1.5000*Fn*mu)/(a*k) + (l*(l*psi + (0.3333*a*psi*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 - 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu) - (0.7500*Fn*mu*(a^2 - (0.1111*a^2*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 - 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max))^2)/(Fn^2*mu^2)))/(a^3*k)))/(l + (0.3333*a*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 - 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu))))/(a^3*k))^(1/2) + (l*psi + (0.3333*a*psi*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 - 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu) - (0.7500*Fn*mu*(a^2 - (0.1111*a^2*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 - 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max))^2)/(Fn^2*mu^2)))/(a^3*k))/(l + (0.3333*a*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 - 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu))))/(Fn*mu);
        xLc = (0.3333*a*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 + 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu);
        xRc = -(0.3333*a*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 - 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu);
        
        for i = 1:n+1
            if x(i) < xR %hátsó csúszó szakasz
                q(i) = -qcrp(i);
            elseif x(i) > xL %első csúszó szakasz
                q(i) = qcrp(i);
            elseif (x(i) > xR) && (x(i) < xRc) %hátsó elfordult parabola
                q(i) = qcrp(i) + mR*x(i) - mR*l;
            elseif (x(i) > xLc) && (x(i) < xL) %első elfordult parabola
                q(i) = -qcrp(i) + mL*x(i) - mL*l;
            else %lineáris szakasz
                q(i) = -psi*(x(i) - l);
            end
        end
        
    case "5P"
        mL = ((0.3333*a*psi*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 + 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu) - l*psi + (0.7500*Fn*mu*(a^2 - (0.1111*a^2*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 + 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max))^2)/(Fn^2*mu^2)))/(a^3*k))/(l - (0.3333*a*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 + 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu));
        mR = -(l*psi + (0.3333*a*psi*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 - 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu) + (0.7500*Fn*mu*(a^2 - (0.1111*a^2*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 - 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max))^2)/(Fn^2*mu^2)))/(a^3*k))/(l + (0.3333*a*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 - 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu));
        xL = (0.3333*a^3*k*((((0.3333*a*psi*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 + 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu) - l*psi + (0.7500*Fn*mu*(a^2 - (0.1111*a^2*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 + 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max))^2)/(Fn^2*mu^2)))/(a^3*k))^2/(l - (0.3333*a*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 + 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu))^2 + (6*Fn*mu*((1.5000*Fn*mu)/(a*k) - (l*((0.3333*a*psi*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 + 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu) - l*psi + (0.7500*Fn*mu*(a^2 - (0.1111*a^2*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 + 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max))^2)/(Fn^2*mu^2)))/(a^3*k)))/(l - (0.3333*a*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 + 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu))))/(a^3*k))^(1/2) + ((0.3333*a*psi*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 + 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu) - l*psi + (0.7500*Fn*mu*(a^2 - (0.1111*a^2*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 + 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max))^2)/(Fn^2*mu^2)))/(a^3*k))/(l - (0.3333*a*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 + 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu))))/(Fn*mu);
        xR = -(0.3333*a^3*k*(((l*psi + (0.3333*a*psi*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 - 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu) + (0.7500*Fn*mu*(a^2 - (0.1111*a^2*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 - 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max))^2)/(Fn^2*mu^2)))/(a^3*k))^2/(l + (0.3333*a*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 - 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu))^2 + (6*Fn*mu*((1.5000*Fn*mu)/(a*k) - (l*(l*psi + (0.3333*a*psi*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 - 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu) + (0.7500*Fn*mu*(a^2 - (0.1111*a^2*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 - 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max))^2)/(Fn^2*mu^2)))/(a^3*k)))/(l + (0.3333*a*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 - 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu))))/(a^3*k))^(1/2) - (l*psi + (0.3333*a*psi*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 - 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu) + (0.7500*Fn*mu*(a^2 - (0.1111*a^2*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 - 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max))^2)/(Fn^2*mu^2)))/(a^3*k))/(l + (0.3333*a*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 - 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu))))/(Fn*mu);
        xLc = (0.3333*a*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 + 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu);
        xRc = -(0.3333*a*((9*Fn^2*mu^2 + 4*a^4*k^2*psi_max^2 - 12*Fn*a*l*k*mu*psi_max*sign(psi_max))^(1/2) - 2*a^2*k*psi_max*sign(psi_max)))/(Fn*mu);
        
        for i = 1:n+1
            if x(i) < xR %hátsó csúszó szakasz
                q(i) = qcrp(i);
            elseif x(i) > xL %első csúszó szakasz
                q(i) = -qcrp(i);
            elseif (x(i) > xR) && (x(i) < xRc) %hátsó elfordult parabola
                q(i) = -qcrp(i) + mR*x(i) - mR*l;
            elseif (x(i) > xLc) && (x(i) < xL) %első elfordult parabola
                q(i) = qcrp(i) + mL*x(i) - mL*l;
            else %lineáris szakasz
                q(i) = -psi*(x(i) - l);
            end
        end
end

%5 szakasz kisimul / 5 sections to three-------------------------------
if (sign(psi)*sign(psi_max) == -1) && (abs(psi) >= abs(psi_max))
        switch STATE
            case "5M"
                STATE = "3M";
                 if sign(epsilon)*sign(epsilon_prev) == -1
                     STATE = "5P";
                 end
            case "5P"
                STATE = "3P";
                if sign(epsilon)*sign(epsilon_prev) == -1
                     STATE = "5M";
                end
        end
end

%fordulás / change of direction----------------------------------------
    if sign(epsilon)*sign(epsilon_prev) == -1 %eps=0
        switch STATE
            case "3P"
                psi_max = (psi + psi_prev)/2;
                STATE = "5M";
            case "3M"
                psi_max = (psi + psi_prev)/2;
                STATE = "5P";
        end
    end
    
state_end = STATE;