clear all
close all
clc

% name of the example to be tested
exName =  'AC4';
% – Load COMPleib test example 
[A,B1,B,C1,C,D11,D12,D21,nx,nw,nu,nz,ny] = COMPleib(exName);
% define Cz and Dz
Cz = C;
Dz = eye(ny,nu);
% check if the methid can be applied
if norm(Dz'*Dz - eye(nu)) < 1e-3
    disp(exName)
    % define the LMIs
    P = sdpvar(nx,nx);
    K = sdpvar(nu,ny);
    M = [P*A+A'*P+Cz'*Dz*K*C+C'*K'*Dz'*Cz, P*B+C'*K';
        B'*P+K*C, -eye(nu)];
    LMI = [P >= 0; M <= 0];
    options = sdpsettings('solver','mosek','verbose',0);
    % solve the LMIs
    sol = optimize(LMI,[],options);
    % check if the determined solution is stabilizing
    KK = value(K);
    if max(real(eig(A+B*KK*C))) < 0
        LMItiming = sol.solvertime;
    else
        LMItiming = -2;
    end
else
    LMItiming = -1;
end