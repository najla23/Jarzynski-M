% coordination of paclitaxel
z = 0;

% Virtual point
zv = 0;

% COM of membrane
zcom = 0;

% Constants
forceconstant = 1000; % kJ/mol/nm^2
kT = 2.48;            % kJ/mol
bn = 100;             % number of bins

% Initialize arrays
c = zeros(bn, 1);
B = zeros(bn, 1);
n = zeros(bn, 1);

% Distance from membrane center
x = z - zcom;
hh = size(x);
fz = forceconstant * (zv - z);

% bin population
for i = 1:hh(1, 1)
    p = floor((x(i, 1) * 10) + bn * 0.5);
    if p >= 1 && p <= bn
        c(p, 1) = 1 + c(p, 1);
        B(p, 1) = fz(i, 1) + B(p, 1);
        n(p, 1) = x(i, 1) + n(p, 1);
    end
end

% define q
q = ((1:99)' - 50) / 10;  % unit: nm

% Symmetrize and calculate BB and NN
BB = zeros(99, 1);
NN = zeros(99, 1);

for i = 1:99
    if (c(100 - i, 1) ~= 0 && c(i, 1) ~= 0)
        BB(i, 1) = (B(i, 1) - B(100 - i, 1)) / (c(100 - i, 1) + c(i, 1));
        NN(i, 1) = (n(i, 1) - n(100 - i, 1)) / (c(100 - i, 1) + c(i, 1));
    else
        BB(i, 1) = 0;
        NN(i, 1) = 0;
    end
end

% Integrate to get free energy
w = cumtrapz(q, BB) / kT;
w1 = cumtrapz(z, fz) / kT;

% Boltzmann weight
ew = exp(-w);
ew1 = exp(-w1);

