function [RStar,tStar,recoverStatus,cvxStatus,pStarRelax,tightness,ZStar,yStar] = ...
    qcqprelaxation_l1norm_no_binary(bearing1,bearing2,quiet)
%QCQPRELAXATION Summary of this function goes here
% This script solves two view geometry using QCQP.
% It first samples groud truth rotation and translations, as well as pairs
% of bearing vectors, then relaxes the problem into an QCQP. From paper:
% "A Certifiably Globally Optimal Solution to                                
% the Non-Minimal Relative Pose Problem"
% Author: Heng Yang           
% Last update: 11/26/2018
% corresponds to the math in discussion.pdf, 
% "convert L1 norm to reduent y_i's"

%% useful constants
DIM_r = 9;
DIM_t = 3;
DIM_z_HOMO = DIM_r + DIM_t + 1; % 13
NUM_POINTS = size(bearing1, 2);

%% construct sphere constraints \tilde{Q}^{S^2}
Q_sphere = {};
Q_t = zeros(DIM_z_HOMO,DIM_z_HOMO);
Q_t(1,1) = -1;
Q_t(2:4, 2:4) = eye(3);
Q_sphere{1} = Q_t;
disp(['Constructed ', num2str(size(Q_sphere, 2)),...
    ' sphere constriants Q_sphere.'])
if ~quiet
    disp(size(Q_sphere{1,1}))
end


%% construct rotation constraints \tilde{Q}_rot
Q_rot = {};
% Orthonormality of rotation rows
for i=1:3
    for j=i:3
        e_ij = zeros(3, 3);
        e_ij(i,j) = 1;
        Q_rot{end+1} = [-kronDelta(i,j), zeros(1, DIM_t), zeros(1, DIM_r);
            zeros(DIM_t, 1), zeros(DIM_t,DIM_t), zeros(DIM_t,DIM_r);
            zeros(DIM_r, 1), zeros(DIM_r, DIM_t), kron(eye(3), (e_ij + e_ij')/2)];
    end
end
% Orthonormality of rotation cols
for i=1:3
    for j=i:3
        e_ij = zeros(3, 3);
        e_ij(i,j) = 1;
        Q_rot{end+1} = [-kronDelta(i,j), zeros(1, DIM_t), zeros(1, DIM_r);
            zeros(DIM_t, 1), zeros(DIM_t,DIM_t), zeros(DIM_t,DIM_r);
            zeros(DIM_r, 1), zeros(DIM_r, DIM_t), kron((e_ij + e_ij')/2, eye(3))];
    end
end
% % right-hand rule on rotation cols (same as rotation rows)
ijk_cycle = {[1,2,3], [2,3,1], [3,1,2]};
% for beta=1:size(ijk_cycle, 2)
%     for alpha=1:3
%         ijk = ijk_cycle{beta};
%         i = ijk(1); j=ijk(2); k= ijk(3);
%         e_k = zeros(3,1); e_k(k) = 1;
%         e_alpha = zeros(3,1); e_alpha(alpha) = 1;
%         e_ij = zeros(3,3); e_ij(i,j) = 1;
%         e_ij_skew = 1/2 * (e_ij - e_ij');
%         Q_rot{end+1} = [0, zeros(1, DIM_t), 1/2 * (kron(e_k, e_alpha))';
%             zeros(DIM_t, 1), zeros(DIM_t, DIM_t), zeros(DIM_t, DIM_r);
%             1/2 * kron(e_k, e_alpha), zeros(DIM_r, DIM_t), kron(e_ij_skew, hatmap(e_alpha))];
%     end
% end
% right-hand rule on rotation rows (same as rotation rows)
for beta=1:size(ijk_cycle, 2)
    for alpha=1:3
        ijk = ijk_cycle{beta};
        i = ijk(1); j=ijk(2); k= ijk(3);
        e_k = zeros(3,1); e_k(k) = 1;
        e_alpha = zeros(3,1); e_alpha(alpha) = 1;
        e_ij = zeros(3,3); e_ij(i,j) = 1;
        e_ij_skew = 1/2 * (e_ij - e_ij');
        Q_rot{end+1} = [0, zeros(1, DIM_t), 1/2 * (kron(e_alpha, e_k))';
            zeros(DIM_t, 1), zeros(DIM_t, DIM_t), zeros(DIM_t, DIM_r);
            1/2 * kron(e_alpha, e_k), zeros(DIM_r, DIM_t), kron(hatmap(e_alpha), e_ij_skew)];
    end
end
disp(['Constructed ', num2str(size(Q_rot, 2)),...
    ' rotation constriant matrices Q_rot.'])
if ~quiet
    disp(size(Q_rot{1,1}))
end

%% form gi(t,r)
Q_A = {};
for i=1:NUM_POINTS
    f2 = bearing2(:, i);
    f1 = bearing1(:, i);
    F1 = hatmap(f1);
    A = kron(f2', F1);
    Q_A{i} = [0, zeros(1, DIM_t), zeros(1, DIM_r);
        zeros(DIM_t, 1), zeros(DIM_t,DIM_t), 1/2 * A;
        zeros(DIM_r, 1), 1/2 * A', zeros(DIM_r, DIM_r)];
end
disp(['Constructed ', num2str(size(Q_A, 2)),...
    ' data matrices Q_A.'])
if ~quiet
    disp(size(Q_A{1,1}))
end
%% Use cvx to solve the QCQP relaxation
if quiet
    disp('set quiet')
    cvx_begin quiet sdp
else
    cvx_begin sdp
end
cvx_precision best
variable ZHomo(DIM_z_HOMO, DIM_z_HOMO) symmetric
variable y(NUM_POINTS)
minimize( sum(y) )
subject to
ZHomo >= 0
ZHomo(1,1) == 1
for i=1:size(Q_sphere, 2)
    trace(Q_sphere{i} * ZHomo) == 0
end
for i=1:size(Q_rot, 2)
    trace(Q_rot{i} * ZHomo) == 0
end
for i=1:NUM_POINTS
    y(i) - trace(Q_A{i} * ZHomo) >= 0
    y(i) + trace(Q_A{i} * ZHomo) >= 0
end
cvx_end

%% recover solution from cvx solution
cvxStatus = cvx_status;
pStarRelax = cvx_optval;
ZStar = ZHomo;
yStar = y;

%{
[V_unsorted, D_unsorted, W_unsorted] = eig(ZStar);
[~,ind] = sort(diag(D_unsorted));
D = D_unsorted(ind,ind);
V = V_unsorted(:,ind);
W = W_unsorted(:,ind);

eigenVals = diag(D);
disp(eigenVals(abs(eigenVals)>1e-6))

%}

RStar=0;
tStar=0;
recoverStatus=0;
tightness=0;



%{

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
%}

%% check tightness
% if recoverStatus
%     gap = zeros(size(RStar,2),1);
%     tightness = true;
%     tolerance = 1e-6;
%     for i=1:size(RStar,2)
%         rHomo = [1; reshape(RHat{i}, [9,1])];
%         tHomo = [1; tStar{i}];
%         XHomo = rHomo * tHomo';
%         zHomo = reshape(XHomo, [DIM_z_HOMO, 1]);
%         fStar = zHomo' * Q_0 * zHomo;
%         if abs(pStarRelax - fStar) > tolerance
%             tightness = false;
%         end
%         gap(i) = pStarRelax - fStar;
%     end
%     disp('Tightness gap:')
%     disp(gap)
% else
%     tightness = false;
% end

end

