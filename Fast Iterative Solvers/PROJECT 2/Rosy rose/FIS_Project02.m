% Multigrid method for poisson
clear;clc;

% Variables that need to be changed: n, gamma, nu1, nu2.
n = 4;    % Grid size. n = 4 or 7
N = 2^n;  % N = 16 or 128
Nc= N / 2; 
h = 1 / N;
H = 1 / Nc;
l = n -1; 
gamma = 2; % Iteration times for each level
nu1 = 1;   % Pre-smoothing steps
nu2 = 1;   % Post-smoothing steps
u0 = zeros(N + 1); % N-1 interial points and 2 boundary points
ua = zeros(N + 1); % Analytical solution
f = zeros(N + 1);  % Right hand side
xgrid = linspace(0, 1, N + 1); % Fine grid
ygrid = linspace(0, 1, N + 1);
xgridc = linspace(0, 1, Nc + 1); % Coarse grid
ygridc = linspace(0, 1, Nc + 1);
if gamma == 1        % Iteration times
    it = l;
else
    it = gamma^l-1;  % it should be equal to count
end
xgridr = linspace(0,it-1,it); % Grid for residual plot

% Calculate analytical solution ua and right hand side f
for i = 1 : N+1
    for j = 1 : N+1
        f(i, j) = 8 * pi^2 * sin(2*pi*xgrid(i)) * sin(2*pi*ygrid(j));
        ua(i, j) = sin(2*pi*xgrid(i)) * sin(2*pi*ygrid(j));
    end
end

% Calculate the numerical solution
rmaxl = zeros(1,it);
tic
[uout, rmaxl,count] = MG(l, u0, f, gamma, nu1, nu2, 0, rmaxl); 
toc
rmaxlr = rmaxl / rmaxl(1);

% Plots
figure(5)
subplot(1,2,1)
surf(uout)
title('Multigrid Method Solution')
subplot(1,2,2)
surf(ua)
title('Analytical Solution')

figure(2)
surf(abs(ua - uout))
title('Error  |u_{a} - u_{MG}|')

figure(3)
semilogy(xgridr,rmaxlr)

% Calculate root mean square error
RMS = 0; 
for i = 1 : N+1
    for j = 1 : N+1
        err = uout(i,j) - ua(i,j);
        RMS = RMS + err^2;
    end
end
RMS = sqrt(RMS);


%% Test

% Find out nu to make norm(u-up,Inf) < 10^-8 and calculate the converged maximum error
nu = 300;  % N = 10: 192； N = 100: 12360； N=16: 458. N=128:19021
up = GS(u0, f, nu - 1); % u for the previous nu
uc = GS(up, f, 1);      % u for the current nu
infnorm = norm(uc-up, Inf);
CME = 0; % Converged maximum error
for i = 1 : N+1
    for j = 1 : N+1
        err = abs(uc(i,j) - ua(i,j));
        if ( CME < err )
            CME = err;
        end
    end
end
CME;

% Restriction operator
uh1 = zeros(N + 1);
for i = 1 : N + 1
    for j = 1 : N + 1
        uh1(i, j) = sin(2*pi*xgrid(i)) * sin(2*pi*ygrid(j));
    end
end
uH1 = RESTR(uh1, Nc);
e2h = 0;
for i = 2 : Nc
    ii = 2 * i - 1;
    for j = 2 : Nc
        jj = 2 * j - 1;
        e = abs(uH1(i,j) - uh1(ii,jj));
        if e2h < e
            e2h = e;
        end
    end
end
e2h;

% Prolongation operator
uH2 = zeros(Nc + 1);
for i = 1 : Nc + 1
    for j = 1 : Nc + 1
        uH2(i, j) = sin(2*pi*xgridc(i)) * sin(2*pi*ygridc(j));
    end
end
uh2 = PROLONG(uH2, Nc);
eh = 0;
for i = 2 : Nc
    ii = 2 * i - 1;
    for j = 2 : Nc
        jj = 2 * j - 1;
        e = abs(uh2(ii,jj) - uH2(i,j));
        if eh < e
            eh = e;
        end
    end
end
eh;


