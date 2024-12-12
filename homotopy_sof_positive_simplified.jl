using HomotopyContinuation, DelimitedFiles, Symbolics, LinearAlgebra, Arblib, BenchmarkTools

# select radius
r = 40;
# initialize outcome
outcome = 0;

# import matrices from csv
Anum = readdlm("./A.csv", ',', Float64);
Bnum = readdlm("./B.csv", ',', Float64);
Cnum = readdlm("./C.csv", ',', Float64);

# matrices dimensions
n_x = size(Anum,1);
n_u = size(Bnum,2);
n_y = size(Cnum,1);

# rationalize the input matrices
A = Anum;
B = Bnum;
C = Cnum;

# define the symbolic variables
@var K[1:n_u, 1:n_y] s
KK = Symbolics.scalarize(K);

# compute closed-loop matrix
Acl = A+B*KK*C;

# compute the characteristic polynomial of the closed loop
p = HomotopyContinuation.expand(det(s*Matrix{Int128}(I, n_x, n_x)-Acl));
listcoeffp = coefficients(p,s)

# construct the Hurwitz matrix
Hp = 0*Acl;
for i = 1:n_x
    for j = 1:n_x
        if (2*j-i >= 0) & (2*j-i <= n_x)
            Hp[i,j] = listcoeffp[2*j-i+1]
        end
    end
end

# define the border polynomial
hp = HomotopyContinuation.expand(det(Hp))

# construct the perimetrical polynomial
var = vcat(K...);
bpol = -r^2;
for i in 1:n_u*n_y
    global bpol
    bpol = bpol + var[i]^2;
end
bpr = hp*bpol;

# to enforce the positivity constraint multiply bpr times the off diagonal terms of Acl
for i = 1:nx
    for j = 1:nx
        if i != j
            global bpr = bpr*Acl[i,j]
        end
    end
end

# construct the map Fepr
Fepr = 0*var;
for i = 1:n_u*n_y
    Fepr[i] = HomotopyContinuation.expand(1.0*HomotopyContinuation.MultivariatePolynomials.differentiate(bpr,var[i]));
end

# solve the system
F = System(Fepr);
solutions = solve(F; compile = true, threading = true);
timing = @benchmark solve(F; compile = true, threading = true);
res = results(solutions)

# keep only real solutions
real_sols = real_solutions(res)

# Use Hurwitz determinants to check if a solution has been found
nsol = size(real_sols,1);
solutionFound = false;
for i = 1:nsol
    global Hps = Hp;
    for j = 1:n_u*n_y
        Hps = HomotopyContinuation.MultivariatePolynomials.subs(Hps, var[j] => real_sols[i][j])
    end
    positiveHps = true
    for k = 1:n_x
        detHn = det(Hps[1:k, 1:k])
        if Float64(detHn) < 0
            positiveHps = false
            break;
        end
    end
    if positiveHps == true
        global solutionFound = true;
        println(real_sols[i][1:n_u*n_y])
        writedlm("solution.csv", real_sols[i][1:n_u*n_y])
    end
end

# save the outcome of the procedure
if outcome == 0
    writedlm("Juliatiming.csv",mean(timing).time*1e-9)
else
    writedlm("Juliatiming.csv",-1)
end