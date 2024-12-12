clear all
close all
clc

% name of the example to be tested
exName =  'AC4';
epsilon = 1e-3;
% – Load COMPleib test example 
[A,B1,B,C1,C,D11,D12,D21,nx,nw,nu,nz,ny] = COMPleib(exName);
% define Cz and Dz
Cz = C;
Dz = eye(ny,nu);
% check if the methid can be applied
W = sdpvar(nx,nx);
M = sdpvar(ny,ny);
N = sdpvar(nu,ny);
LMI = [A*W + W*A' - B*N*C - C'*N'*B' <= -epsilon*eye(nx),...
    W >= epsilon*eye(nx), M*C==C*W];
options = sdpsettings('solver','mosek','verbose',0);
% solve the LMIs
sol = optimize(LMI,[],options);
% check if the determined solution is stabilizing
NN = value(N);
MM = value(M);
KK = - NN/MM;
if max(real(eig(A+B*KK*C))) < 0
    LMItiming = sol.solvertime;
else
    partialTime = sol.solvertime;
    P = sdpvar(nx,nx);
    M = sdpvar(nu,nu);
    N = sdpvar(nu,ny);
    LMI = [P*A + A'*P - B*N*C - C'*N'*B' <= -epsilon*eye(nx),...
        P >= epsilon*eye(nx), B*M==W*B];
    sol = optimize(LMI,[],options);
    NN = value(N);
    MM = value(M);
    KK = - MM\NN;
    if max(real(eig(A+B*KK*C))) < 0
        LMItiming(i) = partialTime + sol.solvertime;
    else
        LMItiming(i) = -2;
    end
end
