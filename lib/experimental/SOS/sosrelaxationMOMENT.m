function [sol,xoptimal,momentdata,sos] = sosrelaxationMOMENT(bearing1,bearing2,R_gt,t_gt)
%UNTITLED Summary of this function goes here
%   solve two-view geometry using SOS
nSamples = size(bearing1, 2);
r = sdpvar(9,1);
t = sdpvar(3,1);
R = reshape(r, [3,3]);
X = r * t';
x = reshape(X, [numel(X), 1]);
C = {3,3};

test = sdpvar(40,1);

for i=1:3
    for j=1:3
        e_i = zeros(3,1); e_i(i) = 1;
        e_j = zeros(3,1); e_j(j) = 1;
        C_ij = zeros(9,9);
        for k=1:nSamples
            f1 = bearing1(:, k);
            f2 = bearing2(:, k);
            C_ij = C_ij + ...
                kron(f2, cross(f1, e_j))*(kron(f2, cross(f1, e_i)))';
        end
        C{i,j} = C_ij;
    end
end
C_extend = zeros(numel(X), numel(X));
for i=1:3
    for j=1:3
        e_ij = zeros(3,3); e_ij(i,j) = 1;
        C_extend = C_extend + kron(e_ij, C{i,j});
    end
end

cost = x' * C_extend * x;
disp('---------cost----------')
disp(cost)

%% copy from QCQP relax
DIM_r = 9;
DIM_t = 3;
DIM_r_HOMO = DIM_r + 1;
DIM_t_HOMO = DIM_t + 1;
DIM_z_HOMO = DIM_r_HOMO * DIM_t_HOMO;
NUM_POINTS = size(bearing1, 2);

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
disp(size(Q_X{1}))


r_homo = [1; r];
t_homo = [1; t];
X_homo = r_homo * t_homo';
x_homo = reshape(X_homo, [numel(X_homo), 1]);

equalities = [];
for i=1:length(Q_sphere)
    equalities = [equalities, test' * Q_sphere{i} * test == 0];
end
for i=1:length(Q_rot)
    equalities = [equalities, test' * Q_rot{i} * test == 0];
end
% for i=1:length(Q_X)
%     equalities = [equalities, test' * Q_X{i} * test == 0];
% end

xIndex = [(DIM_r_HOMO+2):(2*DIM_r_HOMO), ...
    (2*DIM_r_HOMO+2):(3*DIM_r_HOMO), (3*DIM_r_HOMO+2):(4*DIM_r_HOMO)];
Q_0 = zeros(DIM_z_HOMO, DIM_z_HOMO);
for i=1:size(C_extend, 1)
    for j=1:size(C_extend, 2)
        Q_0(xIndex(i), xIndex(j)) = C_extend(i,j);
    end
end
disp(['Constructed cost matrix Q_0 with size: '...
    num2str(size(Q_0, 1)), ', ', num2str(size(Q_0, 2)), '.'])

new_cost = test' * Q_0 * test;

%{
% define SOSP
equalities = [ t'*t - 1 == 0' ];
for i=1:3
    for j=i:3
        if j==i  % Unit norm
            equalities = [equalities; R(:,i)'*R(:,i) - 1 == 0];
        else
            equalities = [equalities; R(:,i)'*R(:,j) == 0];
        end
    end
end
disp(equalities)
%}
% sdisplay(x_homo)
disp(equalities)
degree = 2;

[sol,xoptimal,momentdata,sos] = solvemoment(equalities, new_cost, [], degree);
disp('---------Lower bound------------')
disp(sol)
value(new_cost)

r_gt = reshape(R_gt, [9,1]);
X_gt = r_gt * t_gt';
x_gt = reshape(X_gt, [numel(X_gt), 1]);
cost_gt = x_gt' * C_extend * x_gt;
disp('---------ground truth cost----------')
disp(cost_gt)
end

