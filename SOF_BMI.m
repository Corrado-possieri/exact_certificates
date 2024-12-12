clear all
close all
clc

% name of the example to be tested
exName =  'AC4';
% – Load COMPleib test example 
[A,B1,B,C1,C,D11,D12,D21,nx,nw,nu,nz,ny] = COMPleib(exName);
% Solve a BMI for each COMPleib test example
P = sdpvar(nx,nx);
K = sdpvar(nu,ny);
BMI = [P>=0,(A+B*K*C)'*P+P*(A+B*K*C)>=0];
options = sdpsettings('solver','bmibnb','verbose',1);
options.bmibnb.lpsolver = 'mosek';
% solve the BMIs
sol = optimize(BMI,[],options);
% alternatively, one may use 
% KK = sof(exName);
% Check if the solution is stabilizing
KK = value(K);
if max(real(eig(A+B*KK*C))) < 0
    BMItiming(i) = sol.solvertime;
else
    BMItiming(i) = -2;
end
