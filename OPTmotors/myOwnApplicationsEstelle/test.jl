using Neurthino
using CairoMakie
osc = OscillationParameters(3);

setθ!(osc, 1=>2, 0.59);
setθ!(osc, 1=>3, 0.15);
setθ!(osc, 2=>3, 0.84);
setδ!(osc, 1=>3, 3.86);

setΔm²!(osc, 2=>3, -2.523e-3);
setΔm²!(osc, 1=>2, -7.39e-5);

U = PMNSMatrix(osc)
H = Hamiltonian(osc)

zenith = acos.(range(-1,stop=0,length=200));
paths = Neurthino.prempath(zenith, 2.5, samples=100, discrete_densities=0:0.1:14);
energies = 10 .^ range(0, stop=2, length=200);
prob = Pνν(U, H, energies, paths);
probs = prob[:,:,1,1]
matprobs=parent(probs)
