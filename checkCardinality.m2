-- define the ambient ring
R = QQ[K11,K12,K21,K22,l]

-- load the polynomial from the workspace
bpr = value get("bpr.txt");

-- define the ideals
I = ideal diff(matrix{{K11,K12,K21,K22}}, bpr)
J = I + ideal(1+l*bpr)

-- load the numeric solutions package
loadPackage("NumericSolutions")

-- compute the trace form
T = traceForm(J)

-- compute the eigenvalues of the trace form
eigT = eigenvalues(substitute(T,QQ))

-- return the eigenvalues
g = openOut("eigenvalues_trace.txt")
g << toString(eigT) << close