function [RHat,tHat,RStar,tStar,ZStar,recoverStatus,cvxStatus,pStarRelax,tightness,rankZ] = ...
    qcqprelaxation(bearing1,bearing2,options)
%QCQPRELAXATION Solve Two view geometry using SDP relaxation
% This script solves two view geometry using QCQP.
% From paper:
% "A Certifiably Globally Optimal Solution to                                
% the Non-Minimal Relative Pose Problem"
% Author: Heng Yang           
% Last update: 1/31/2019

if ~isfield(options,'CVXDisp')
    options.CVXDisp=false;
end

if ~isfield(options,'QCQPDisp')
    options.QCQPDisp=false;
end

%% useful constants
DIM_r = 9;
DIM_t = 3;
DIM_r_HOMO = DIM_r + 1;
DIM_t_HOMO = DIM_t + 1;
DIM_z_HOMO = DIM_r_HOMO * DIM_t_HOMO;
NUM_POINTS = size(bearing1, 2);

%% construct lifted sphere constraints \hat{C}_t
PHomo = [-1, zeros(1,3); zeros(3,1), eye(3)];
Q_sphere = {};
for i=1:DIM_r_HOMO
    for j=i:DIM_r_HOMO
        e_ij = zeros(DIM_r_HOMO, DIM_r_HOMO);
        e_ij(i,j) = 1;
        Q_sphere{end+1} = kron(PHomo, (e_ij + e_ij')/2);
    end
end
if options.QCQPDisp
    fprintf('Constructed %g lifted sphere constraints Q_sphere.\n', size(Q_sphere,2));
    fprintf('Dimension of each matrix inside Q_sphere is %g by %g.\n',...
        size(Q_sphere{1,1},1), size(Q_sphere{1,1},2));
    fprintf('Dimension of matrix homogeneous P is %g by %g.\n',size(PHomo,1),size(PHomo,2));
end

%% construct lifted rotation constraints \hat{C}_R
% construct matrices in C_R
PRotHomo = {};
% Orthonormality of rotation columns
for i=1:DIM_t
    for j=i:DIM_t
        e_ij = zeros(DIM_t, DIM_t);
        e_ij(i,j) = 1;
        PRotHomo{end+1} = [-kronDelta(i,j), zeros(1, DIM_t^2);
            zeros(DIM_t^2, 1), kron((e_ij + e_ij')/2, eye(DIM_t))];
    end
end
% Orthonormality of rotation rows
for i=1:DIM_t
    for j=i:DIM_t
        e_ij = zeros(DIM_t, DIM_t);
        e_ij(i,j) = 1;
        PRotHomo{end+1} = [-kronDelta(i,j), zeros(1, DIM_t^2);
            zeros(DIM_t^2, 1), kron(eye(DIM_t), (e_ij + e_ij')/2)];
    end
end
% right-hand rule on rotation cols (same as rotation rows)
ijk_cycle = {[1,2,3], [2,3,1], [3,1,2]};
for beta=1:size(ijk_cycle, 2)
    for alpha=1:3
        ijk = ijk_cycle{beta};
        i = ijk(1); j=ijk(2); k= ijk(3);
        e_k = zeros(3,1); e_k(k) = 1;
        e_alpha = zeros(3,1); e_alpha(alpha) = 1;
        e_ij = zeros(3,3); e_ij(i,j) = 1;
        e_ij_skew = 1/2 * (e_ij - e_ij');
        PRotHomo{end+1} = [0, 1/2 * (kron(e_k, e_alpha))';
            1/2 * kron(e_k, e_alpha), kron(e_ij_skew, hatmap(e_alpha))];
    end
end
if options.QCQPDisp
    fprintf('Constructed %g rotation constraint matrices PRotHomo, each with dimension %g by %g.\n',...
        size(PRotHomo,2),size(PRotHomo{1,1},1),size(PRotHomo{1,1},2));
end
% construct lifted rotation constraints
Q_rot = {};
numPRot = size(PRotHomo, 2);
for i=1:DIM_t_HOMO
    for j=i:DIM_t_HOMO
        for k = 1:numPRot
            e_ij = zeros(DIM_t_HOMO, DIM_t_HOMO);
            e_ij(i,j) = 1;
            P_k = PRotHomo{k};
            Q_rot{end+1} = kron((e_ij+e_ij')/2, P_k);
        end
    end
end
if options.QCQPDisp
    fprintf('Constructed %g lifted rotation constraints Q_rot, each with size %g by %g.\n',...
        size(Q_rot, 2),size(Q_rot{1,1},1),size(Q_rot{1,1},2));
end   
        
%% construct lifted auxiliary X constraints \hat{C}_X
Q_X = {};
for i=1:DIM_r_HOMO
    for i_=(i+1):DIM_r_HOMO
        for j=1:DIM_t_HOMO
            for j_=(j+1):DIM_t_HOMO
                e_jj_ = zeros(DIM_t_HOMO, DIM_t_HOMO);
                e_jj_(j, j_) = 1;
                e_ii_ = zeros(DIM_r_HOMO, DIM_r_HOMO);
                e_ii_(i, i_) = 1;
                e_jj_skew = 1/2 * (e_jj_ - e_jj_');
                e_ii_skew = 1/2 * (e_ii_ - e_ii_');
                Q_X{end+1} = kron(e_jj_skew, e_ii_skew);
            end
        end
    end
end
if options.QCQPDisp
    fprintf('Constructed %g lifted auxiliary constraints Q_X, each with size %g by %g.\n',...
        size(Q_X,2),size(Q_X{1,1},1),size(Q_X{1,1},2));
end

%% construct cost function matrix Q_0
% Construct C_extend corresponding to vec(X)
C = {3,3};
for i=1:3
    for j=1:3
        e_i = zeros(3,1); e_i(i) = 1;
        e_j = zeros(3,1); e_j(j) = 1;
        C_ij = zeros(9,9);
        for k=1:NUM_POINTS
            f1 = bearing1(:, k);
            f2 = bearing2(:, k);
            C_ij = C_ij + ...
                kron(f2, cross(f1, e_j))*(kron(f2, cross(f1, e_i)))';
        end
        C{i,j} = C_ij;
    end
end
C_extend = zeros(DIM_r*DIM_t, DIM_r*DIM_t);
for i=1:3
    for j=1:3
        e_ij = zeros(3,3); e_ij(i,j) = 1;
        C_extend = C_extend + kron(e_ij, C{i,j});
    end
end
if options.QCQPDisp
    fprintf('Constructed cost matrix C_extend with size %g by %g.\n',...
        size(C_extend,1),size(C_extend,2));
end

% permutation C_extend in the right place of Q_0, using
% z_homo = vec([1, t'; r, X])
xIndex = [(DIM_r_HOMO+2):(2*DIM_r_HOMO), ...
    (2*DIM_r_HOMO+2):(3*DIM_r_HOMO), (3*DIM_r_HOMO+2):(4*DIM_r_HOMO)];
Q_0 = zeros(DIM_z_HOMO, DIM_z_HOMO);
for i=1:size(C_extend, 1)
    for j=1:size(C_extend, 2)
        Q_0(xIndex(i), xIndex(j)) = C_extend(i,j);
    end
end
if options.QCQPDisp
    fprintf('Constructed cost matrix Q_0 with size %g by %g.\n',...
        size(Q_0,1),size(Q_0,2));
end

%% Use cvx to squietolve the QCQP
% cvx_solver sedumi
cvx_solver mosek
if ~options.CVXDisp
    if options.QCQPDisp
        disp('CVX set quiet.')
    end
    cvx_begin quiet sdp
else
    cvx_begin sdp
end
cvx_precision best
variable ZHomo(DIM_z_HOMO, DIM_z_HOMO) symmetric
minimize( trace(Q_0 * ZHomo) )
subject to
ZHomo >= 0
ZHomo(1,1) == 1
for i=1:size(Q_sphere, 2)
    trace(Q_sphere{i} * ZHomo) == 0
end
for i=1:size(Q_rot, 2)
    trace(Q_rot{i} * ZHomo) == 0
end
for i=1:size(Q_X, 2)
    trace(Q_X{i} * ZHomo) == 0
end
cvx_end

%% recover solution from cvx solution
cvxStatus = cvx_status;
pStarRelax = cvx_optval;
ZStar = ZHomo;

[V_unsorted, D_unsorted, W_unsorted] = eig(ZStar);
[~,ind] = sort(diag(D_unsorted));
D = D_unsorted(ind,ind);
V = V_unsorted(:,ind);
W = W_unsorted(:,ind);

eigenVals = diag(D);

% disp(eigenVals)
lambdas = eigenVals((DIM_z_HOMO-4):end);
lambda_first = eigenVals(1);
epsilon = 1e-3;

rankZ=sum(abs(eigenVals)>epsilon);

if options.QCQPDisp
    fprintf('Rank of Z=%g.\n',rankZ);
    fprintf('The largest four eigen values of Z are: ');
    fprintf('%g ',lambdas(2:end)); fprintf('\n');
    fprintf('The fifth largest eigen values is %g, and the smallest eigen value is %g.\n',...
        lambdas(1),lambda_first);
end

RStar = {};
RHat = {};
tStar = {};
tHat = {};
for i=1:4
    % recover R_hat and t_hat
    v = V(:, DIM_z_HOMO-i+1);
    v_1 = v(1);
    v_r = v(2:10) * sign(v_1);
    v_t = v([11, 21, 31]) * sign(v_1);
    R_hat = reshape(sqrt(3)*v_r/norm(v_r), [3,3]);
    RHat{end+1} = R_hat;
    % project R_hat to SO(3)
    RStar{end+1} = projToSO3(R_hat);
    tStar{end+1} = v_t / norm(v_t);
    tHat{end+1} = v_t / norm(v_t);
end     

if abs(lambdas(1)) < epsilon && abs(lambda_first) < epsilon ...
        && abs(lambdas(2)-lambdas(3)) < epsilon ...
        && abs(lambdas(4)-lambdas(5)) < epsilon
    if options.QCQPDisp    
        fprintf('Recover status true: there are only four nonzeros eigen values with two distinctive values.\n');
    end
    recoverStatus = true;
else
    warning('Rank of solution matrix Z is not 4. Obtain rounded solutions.\n');
    recoverStatus = false;
end

%% check tightness
gap = zeros(size(RStar,2),1);
tightness = true;
tolerance = 1e-6;
for i=1:size(RStar,2)
    rHomo = [1; reshape(RHat{i}, [9,1])];
    tHomo = [1; tStar{i}];
    XHomo = rHomo * tHomo';
    zHomo = reshape(XHomo, [DIM_z_HOMO, 1]);
    fStar = zHomo' * Q_0 * zHomo;
    if abs(pStarRelax - fStar) > tolerance
        tightness = false;
    end
    gap(i) = pStarRelax - fStar;
end
if options.QCQPDisp
    fprintf('CVX return status: %s, CVX optimal cost: %g.', cvxStatus,pStarRelax); 
    fprintf('Tightness gaps are: '); fprintf('%g ',gap); fprintf('\n');
    fprintf('Relaxation tightness: %s.\n',string(tightness));
end
end

