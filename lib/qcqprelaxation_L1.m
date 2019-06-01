function [RStar,tStar,ZStar,recoverStatus,cvxStatus,pStarRelax,tightness] = ...
    qcqprelaxation_L1(bearing1,bearing2,quiet)
%QCQPRELAXATION Summary of this function goes here
% This script solves two view geometry using QCQP.
% It first samples groud truth rotation and translations, as well as pairs
% of bearing vectors, then relaxes the problem into an QCQP. From paper:
% "A Certifiably Globally Optimal Solution to                                
% the Non-Minimal Relative Pose Problem"
% Author: Heng Yang           
% Last update: 11/19/2018

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
disp(['Constructed ', num2str(size(Q_sphere, 2)),...
    ' lifted sphere constriants Q_sphere.'])
if ~quiet
    disp(size(Q_sphere{1,1}))
    disp(size(PHomo))
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
disp(['Constructed ', num2str(size(PRotHomo, 2)),...
    ' rotation constriant matrices PRotHomo.'])
if ~quiet
    disp(size(PRotHomo{1,1}))
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
disp(['Constructed ', num2str(size(Q_rot, 2)),...
    ' lifted rotation constriants Q_rot.'])
if ~quiet
    disp(size(Q_rot{1,1}))
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
disp(['Constructed ', num2str(size(Q_X, 2)),...
    ' lifted X constriants Q_X.'])
if ~quiet
    disp(size(Q_X{1,1}))
end

% r = sym('r', [9,1], 'real');
% r_homo = [1; r];
% t = sym('t', [3,1], 'real');
% t_homo = [1; t];
% X_homo = r_homo * t_homo';
% x_homo = reshape(X_homo, [numel(X_homo), 1]);
% test = sym('test', [40,1], 'real');
% for i=1:length(Q_X)
%     disp(x_homo' * Q_X{i} * x_homo)
% end

%% construct cost function matrix Q_0
tIndex = [11, 21, 31];
rIndex = 2:10;
Q_0 = {};
for i=1:NUM_POINTS
    f1 = bearing1(:, k);
    f2 = bearing2(:, k);
    A = kron(f2', hatmap(f1));
    Q_i = zeros(DIM_z_HOMO, DIM_z_HOMO);
    Q_i(tIndex, rIndex) = A/2;
    Q_i(rIndex, tIndex) = A'/2;
    Q_0{end+1} = Q_i;
end

%% Use cvx to solve the QCQP
if quiet
    disp('set quiet')
    cvx_begin quiet sdp
else
    cvx_begin sdp
end
cvx_precision best
variable ZHomo(DIM_z_HOMO, DIM_z_HOMO) symmetric
variable y(NUM_POINTS)
% minimize( trace(Q_0 * ZHomo) )
minimize( sum(y) + norm_nuc(ZHomo(2:4,2:4)) + norm_nuc(ZHomo(5:13,5:13)))
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
for i=1:NUM_POINTS
    y(i) >= abs(trace(Q_0{i} * ZHomo))
end
cvx_end

%% recover solution from cvx solution
cvxStatus = cvx_status;
pStarRelax = cvx_optval;
ZStar = ZHomo;

% [V_unsorted, D_unsorted, W_unsorted] = eig(ZStar);
% [~,ind] = sort(diag(D_unsorted));
% D = D_unsorted(ind,ind);
% V = V_unsorted(:,ind);
% W = W_unsorted(:,ind);

figure
spy(ZHomo)
figure
bar(eig(full(ZHomo)))

RStar=0,
tStar=0,
recoverStatus=false;
tightness=0;

%{
eigenVals = diag(D);
% disp(eigenVals)
lambdas = eigenVals((DIM_z_HOMO-4):end);
lambda_first = eigenVals(1);
epsilon = 1e-3;
disp('--------Eigen Values-----------------')
disp(lambdas)
disp(lambda_first)
disp('-------------------------------------')
if abs(lambdas(1)) < epsilon && abs(lambda_first) < epsilon ...
        && abs(lambdas(2)-lambdas(3)) < epsilon ...
        && abs(lambdas(4)-lambdas(5)) < epsilon
    disp('low-rank eigen-decomposition is correct.')
    recoverStatus = true;
    RStar = {};
    RHat = {};
    tStar = {};
    for i=1:4
        % recover R_hat and t_hat
        v = V(:, DIM_z_HOMO-i+1);
        v_1 = v(1);
        v_r = v(2:10) * sign(v_1);
        v_t = v([11, 21, 31]) * sign(v_1);
%         disp('v_t')
%         disp(v_t)
        R_hat = reshape(sqrt(3)*v_r/norm(v_r), [3,3]);
        RHat{end+1} = R_hat;
        % project R_hat to SO(3)
        RStar{end+1} = projToSO3(R_hat);
        tStar{end+1} = v_t / norm(v_t);
    end     
else
    disp('Fail to recover rotation and translation.')
    recoverStatus = false;
    RStar = {};
    tStar = {};
end

%% check tightness
if recoverStatus
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
    disp('Tightness gap:')
    disp(gap)
else
    tightness = false;
end

%}
end

