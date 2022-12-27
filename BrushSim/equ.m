function [yprime, ymod, M_mod] = equ(t, y, dt)
% Brush tyre with sliding at zero speed
% y = [psi; psid; qi]

global x dx l k J b n qcrp eps Mdes qdiff

psi = y(1);
ome = y(2);
q = y(3:n+3);

qd = zeros(1, n+1);
dq = zeros(1, n+1);
dq(n+1) = sign(q(n+1))*qdiff*x(n+1)*dx;

for i=1:n+1
    if abs(q(i)) >= qcrp(i)
        q(i) = (qcrp(i)-eps)*sign(q(i));
        qd(i) = 0;
    end
    if i ~= n+1
        dq(i) = q(i+1) - q(i);
    end
    if abs(q(i)) < qcrp(i)
        qd(i) = (l-x(i))*ome;
    end
end

defmom = 0;
defmomd=0;
for i=2:n+1
    defmom = defmom + (l - (x(i) + x(i-1))/2)*(q(i) + q(i-1))/2*dx;
    defmomd=defmomd+(l-(x(i)+x(i-1))/2)*(qd(i)+qd(i-1))/2*dx;
end

M = interp1(Mdes(:,1), Mdes(:,2), t, 'linear', 'extrap');

epsilon = -k/J*defmom - b/J*ome+ M/J;
M_mod = M - J*epsilon - b*ome;

ymod = [psi, ome, q];
yprime = [ome, epsilon, qd];


end

