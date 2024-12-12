using HomotopyContinuation, DelimitedFiles, Symbolics, LinearAlgebra, Arblib, BenchmarkTools

# select radius
r = 10;
# initialize outcome
outcome = 0;


# define the symbolic variables
@var b1 s l

# define the two polynomials
p1 = s^3+s^2+3*s+b1-5;
p2 = 2*s^3-s^2+11*s+2*b1-12;
# compute their product
p = p1*p2
# determine its coefficients
listcoeffp = coefficients(p,s)
n_x = degree(p);

# construct the Hurwitz matrix
Hp = 0*s*Matrix{Int128}(I, n_x, n_x);
for i = 1:n_x
    for j = 1:n_x
        if (2*j-i >= 0) & (2*j-i <= n_x)
            Hp[i,j] = listcoeffp[2*j-i+1]
        else
            Hp[i,j] = 0*listcoeffp[1]
        end
    end
end

# define the border polynomial
hp = HomotopyContinuation.expand(det(Hp))

# construct the perimetrical polynomial
var = [b1];
bpol = -r^2;
for i in 1:1
    global bpol
    bpol = bpol + var[i]^2;
end
bpr = hp*bpol;

# construct the map Fepr
Fepr = 0*var;
for i = 1:1
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
    for j = 1:1
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