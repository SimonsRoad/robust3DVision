function [RStar,tStar,recoverStatus,cvxStatus,pStarRelax,tightness,WStar] = ...
    qcqprelaxation_l1norm(bearing1,bearing2,quiet)
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
DIM_r_HOMO = DIM_r + 1;  % 10
DIM_t_HOMO = DIM_t + 1;  % 4
DIM_z_HOMO = DIM_r_HOMO * DIM_t_HOMO; % 40
NUM_POINTS = size(bearing1, 2);
DIM_w_HOMO = DIM_z_HOMO + NUM_POINTS * 2;

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

Q_sphere_L1 = {};
for i=1:length(Q_sphere)
    temp = zeros(DIM_w_HOMO, DIM_w_HOMO);
    temp(1:DIM_z_HOMO, 1:DIM_z_HOMO) = Q_sphere{i};
    Q_sphere_L1{i} = temp;
end

if ~quiet
    disp(size(Q_sphere_L1{1,1}))
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

Q_rot_L1 = {};
for i=1:length(Q_rot)
    temp = zeros(DIM_w_HOMO, DIM_w_HOMO);
    temp(1:DIM_z_HOMO, 1:DIM_z_HOMO) = Q_rot{i};
    Q_rot_L1{i} = temp;
end

if ~quiet
    disp(size(Q_rot_L1{1,1}))
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

Q_X_L1 = {};
for i=1:length(Q_X)
    temp = zeros(DIM_w_HOMO, DIM_w_HOMO);
    temp(1:DIM_z_HOMO, 1:DIM_z_HOMO) = Q_X{i};
    Q_X_L1{i} = temp;
end

if ~quiet
    disp(size(Q_X_L1{1,1}))
end


%% form constriants on binary variable b1,b2,...,bN, bi^2 - 1 = 0
Q_b = {};
for i=1:NUM_POINTS
    b_c = zeros(DIM_w_HOMO, DIM_w_HOMO);
    b_c(1,1) = -1;
    b_c(DIM_z_HOMO+NUM_POINTS+i, DIM_z_HOMO+NUM_POINTS+i) = 1;
    Q_b{i} = b_c;
end

disp(['Constructed ', num2str(size(Q_b, 2)),...
    ' lifted binary variable constriants Q_b.'])


%% form constaints on yi = bi*gi(x), -yi + bi*gi(x) = 0
% w_homo = [1; r; t1; t1*r; t2; t2*r; t3; t3*r; y; b]
xIndex = [(DIM_r_HOMO+2):(2*DIM_r_HOMO), ...
    (2*DIM_r_HOMO+2):(3*DIM_r_HOMO), (3*DIM_r_HOMO+2):(4*DIM_r_HOMO)];
Q_abs = {};
for i=1:NUM_POINTS
    f2 = bearing2(:, i);
    f1 = bearing1(:, i);
    F1 = hatmap(f1);
    C = kron(f2', F1);
    A = reshape(C', [numel(C), 1]);
    A = A';
    abs_c = zeros(DIM_w_HOMO, DIM_w_HOMO);
    abs_c(1, DIM_z_HOMO+i) = -1/2;
    abs_c(DIM_z_HOMO+i, 1) = -1/2;
    abs_c(DIM_z_HOMO+NUM_POINTS+i, xIndex) = A/2;
    abs_c(xIndex, DIM_z_HOMO+NUM_POINTS+i) = A'/2;
    Q_abs{i} = abs_c;
end

disp(['Constructed ', num2str(size(Q_abs, 2)),...
    ' absolute value constriants Q_abs.'])

%% form positivibity constraint yi>=0
Q_positive = {};
for i=1:NUM_POINTS
    pos_c = zeros(DIM_w_HOMO, DIM_w_HOMO);
    pos_c(DIM_z_HOMO+i, 1) = 1/2;
    pos_c(1, DIM_z_HOMO+i) = 1/2;
    Q_positive{i} = pos_c;
end
disp(['Constructed ', num2str(size(Q_positive, 2)),...
    ' positive constriants Q_positive.'])

%% form cost function cost = y1+y2+...+yN
Q_0 = zeros(DIM_w_HOMO, DIM_w_HOMO);
for i=1:length(Q_positive)
    Q_0 = Q_0 + Q_positive{i};
end

disp('Constructed cost function Q_0.')

%% Use cvx to solve the QCQP relaxation
if quiet
    disp('set quiet')
    cvx_begin quiet sdp
else
    cvx_begin sdp
end
cvx_precision best
variable WLiftedHomo(DIM_w_HOMO, DIM_w_HOMO) symmetric
minimize( trace(Q_0 * WLiftedHomo) )
subject to
WLiftedHomo >= 0
WLiftedHomo(1,1) == 1
for i=1:size(Q_sphere_L1, 2)
    trace(Q_sphere_L1{i} * WLiftedHomo) == 0
end
for i=1:size(Q_rot_L1, 2)
    trace(Q_rot_L1{i} * WLiftedHomo) == 0
end
for i=1:size(Q_X_L1, 2)
    trace(Q_X_L1{i} * WLiftedHomo) == 0
end
for i=1:size(Q_b, 2)
    trace(Q_b{i} * WLiftedHomo) == 0
end
for i=1:size(Q_abs, 2)
    trace(Q_abs{i} * WLiftedHomo) == 0
end
for i=1:size(Q_positive, 2)
    trace(Q_positive{i} * WLiftedHomo) >= 0
end
cvx_end

%% recover solution from cvx solution
cvxStatus = cvx_status;
pStarRelax = cvx_optval;
WStar = WLiftedHomo;

%{
[V_unsorted, D_unsorted, W_unsorted] = eig(WStar);
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

