function [R_est,t_est,E_est,f_sdp,f_round,Z_sdp,solverOutput] = estimateEssentialSDP(bearing1,bearing2,varargin)
% implement the essential matrix estimation algorithm as proposed by Ji
% Zhao in the paper:
% An Efficient Solution to Non-minimal Case Essential Matrix Estimation
% Use SDP relaxation to solve the QCQP
% Author: Heng Yang
% Last update: 06/02/2019
% Massachusetts Institute of Technology


%% Build cost matrix C from bearing1 and bearing2
dim_e = 9;
dim_t = 3;
dim_x = dim_e + dim_t;
C = zeros(dim_e,dim_e);
% first calculate C based on Eq. (12)
% TODO

% Build Q from C by adding zeros
Q = zeros(dim_x,dim_x);
Q = [C, zeros(dim_e,dim_t); zeros(dim_t,dim_e), zeros(dim_t,dim_t)];

%% Generate constraint matrices A (should be 7 constant matrices)
% x = [e\tran, t\tran]\tran, dimension 12, each Q should also be dimension
% 12 by 12
numCons = 7;
A = zeros(dim_x,dim_x,numCons);
% Manually build each of the 7 sparse constraint matrices
% build A1 (15a)
A1 = zeros(dim_x,dim_x);
A1(1,1)=1; A1(4,4)=1; A1(7,7)=1; A1(11,11)=-1; A1(12,12)=-1;
A(:,:,1) = A1;
% build A2 (15b)




%% Solve using CVX interface
cvx_solver mosek % sedumi, sdpt3
cvx_begin sdp
cvx_precision best
variable Z(dim_x,dim_x) symmetric
% the cost function
minimize( trace(Q * Z) )
% the constraints
subject to
Z == semidefinite(dim_x);
for i=1:numCons
    trace(A(:,:,i) * Z) == 0;
end
cvx_end

%% Obtain solution from SDP result, calcualte rank 
Z_sdp = Z;
rank_Z_sdp = rank(Z_sdp,1e-3);
if rank_Z_sdp == 1
    fprintf('Rank = 1, SDP solves the original QCQP.\n');
else
    warning('Rank = %g, possibly wrong solution.\n',rank_Z_sdp);
end
% if rank=1, do eigen decomposition and recover x, e, and t


% recover R from e and t (should have 4 pairs of solutions)
