using HomotopyContinuation, DelimitedFiles, Symbolics, LinearAlgebra, Arblib, BenchmarkTools

# select radius
r = 100;
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
A = Array{Rational}(undef, n_x, n_x)
B = Array{Rational}(undef, n_x, n_u)
C = Array{Rational}(undef, n_y, n_x)
for i = 1:n_x
    for j = 1:n_x
        A[i,j] = rationalize(Anum[i,j],tol=1e-6);
    end
    for j = 1:n_u
        B[i,j] = rationalize(Bnum[i,j],tol=1e-6);
    end
    for j = 1:n_y
        C[j,i] = rationalize(Cnum[j,i],tol=1e-6);
    end
end

# define the symbolic variables
@var K[1:n_u, 1:n_y] s l
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

# construct the map Fepr
Fepr = 0*var;
for i = 1:n_u*n_y
    Fepr[i] = HomotopyContinuation.expand(1.0*HomotopyContinuation.MultivariatePolynomials.differentiate(bpr,var[i]));
end
push!(Fepr,HomotopyContinuation.expand(1+1.0*l*bpr));

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
        Hps = HomotopyContinuation.MultivariatePolynomials.subs(Hps, var[j] => rationalize(real_sols[i][j]))
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

# if a soluion has not been found, then check if a solution exists
isSolutionSetCertifiedEmpty = true;
numberComplexSolutions = 0;
if solutionFound == false
    numberRsol = 0;
    # certification procedure
    certs = certificates(certify(F, solutions))
    hInt = Matrix(undef, length(certs), 1);
    candSol = [];
    # check if all the conditions hold
    for i in eachindex(certs)
        if is_certified(certs[i])
            global numberComplexSolutions = numberComplexSolutions + 1;
            if is_real(certs[i])
                interSol = real(certs[i].I);
                global numberRsol = numberRsol+1;
                if Arblib.is_finite(interSol[end]) == false
                    global isSolutionSetCertifiedEmpty = false;
                    println("The algorithm failed since the interval is too large")
                    global outcome = 1;
                    break
                end
            end
        else
            global isSolutionSetCertifiedEmpty = false
            println("The algorithm failed due to certification issues")
            global outcome = 2;
            break
        end
    end

    if isSolutionSetCertifiedEmpty == true
        open("bpr.txt", "w") do file
            write(file, string(bpr))
        end
        print("If the cardinality is ")
        print(numberRsol)
        println(",then there does not exist a solution")
        global outcome = 3;
    end
end

# save the outcome of the procedure
if outcome == 0
    writedlm("Juliatiming.csv",mean(timing).time*1e-9)
elseif outcome == 1
    writedlm("Juliatiming.csv",-1)
elseif outcome == 2
    writedlm("Juliatiming.csv",-2)
elseif outcome == 3
    writedlm("Juliatiming.csv",mean(timing).time*1e-9)
    writedlm("NoSol.csv",numberRsol)
    writedlm("numberComplexSolutions.csv",numberComplexSolutions)
end