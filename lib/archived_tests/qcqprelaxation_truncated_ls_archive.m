function [cvxStatus,pStarRelax,ZStar] = ...
    qcqprelaxation_truncated_ls(bearing1,bearing2,R_gt,t_gt,quiet)
% See truncated_LS.pdf for the corresponding math derivation
%% useful constants
NUM_POINTS = size(bearing1, 2);
DIM_r = 9;
DIM_t = 3;
DIM_tr = DIM_t * DIM_r;
DIM_theta = NUM_POINTS;
DIM_z_HOMO = 28*DIM_theta + 40;

%% construct lifted rotation constraints Q_rot
% first construct unlifted rotational constraints P_rot
P_rot = {};
% orthonormal rows
for i=1:3
    for j=i:3
        e_ij = zeros(3, 3); e_ij(i,j) = 1;
        P_rot{end+1} = [-kronDelta(i,j), zeros(1, 3*3);
            zeros(3*3, 1), kron( eye(3), (e_ij + e_ij')/2 )];
    end
end
% orthonormal cols
for i=1:3
    for j=i:3
        e_ij = zeros(3, 3); e_ij(i,j) = 1;
        P_rot{end+1} = [-kronDelta(i,j), zeros(1, 3*3);
            zeros(3*3, 1), kron((e_ij + e_ij')/2, eye(3))];
    end
end
% right-hand
ijk_cycle = {[1,2,3], [2,3,1], [3,1,2]};
for beta=1:size(ijk_cycle, 2)
    for alpha=1:3
        ijk = ijk_cycle{beta};
        i = ijk(1); j=ijk(2); k= ijk(3);
        e_k = zeros(3,1); e_k(k) = 1;
        e_alpha = zeros(3,1); e_alpha(alpha) = 1;
        e_ij = zeros(3,3); e_ij(i,j) = 1;
        e_ij_skew = 1/2 * (e_ij - e_ij');
        P_rot{end+1} = [0, 1/2 * (kron(e_k, e_alpha))';
            1/2 * kron(e_k, e_alpha), kron(e_ij_skew, hatmap(e_alpha))];
    end
end
% then create lifted constraints
Q_rot = {};
% first, create lifted version by multiplying t_i and t_j
% create appropriate indexing from r*t_i to z
tr2zIndex = {};
for i=1:4
    if i==1 % 1*\tilde{r} = 1, r
        tr2zIndex{i} = [1, (1+DIM_t+1):(1+DIM_t+DIM_r)];
    else % t_i * \tilde{r} = t_i, t_i*r
        tr2zIndex{i} = [i,...
            (1+DIM_t+DIM_r+DIM_theta+(i-2)*DIM_r+1):...
            (1+DIM_t+DIM_r+DIM_theta+(i-2)*DIM_r+DIM_r)];
    end
end        
for i=1:4
    for j=i:4
        for k=1:length(P_rot)
            rowIndex = tr2zIndex{i};
            colIndex = tr2zIndex{j};
            temp = zeros(DIM_z_HOMO, DIM_z_HOMO);
            temp(rowIndex, colIndex) = P_rot{k};
            Q_rot{end+1} = temp;
        end
    end
end
% second, create lifted version by multiplying t_i and theta_j
% create appropriate indexing from t_i * theta_j * r to z
tthetar2zIndex = {};
for i=1:DIM_t
    for j=1:DIM_theta
        temp = 1+DIM_t+DIM_r+DIM_theta+DIM_tr+(j-1)*DIM_tr+(i-1)*DIM_r+1:...
            1+DIM_t+DIM_r+DIM_theta+DIM_tr+(j-1)*DIM_tr+(i-1)*DIM_r+DIM_r;
        tthetar2zIndex{i,j} = temp;
    end
end           
for i=1:DIM_t
    for j=1:DIM_theta
        for k=1:length(P_rot)
            P_temp = P_rot{k};
            P_topleft = P_temp(1,1);
            P_topright = P_temp(1,2:end);
            P_bottomleft = P_temp(2:end,1);
            P_bottomright = P_temp(2:end, 2:end);
            Q_temp = zeros(DIM_z_HOMO, DIM_z_HOMO);
            Q_temp(1+DIM_t+1:1+DIM_t+DIM_r, tthetar2zIndex{i,j}) = P_bottomright;
            Q_temp(1, tthetar2zIndex{i,j}) = P_topright;
            Q_temp(tthetar2zIndex{i,j}, 1) = P_bottomleft;
            Q_temp(1+i, 1+DIM_t+DIM_r+j) = P_topleft;
            Q_rot{end+1} = Q_temp;
        end
    end
end
% third, create lifted version by multiplying t_i t_j and theta_u
for i=1:DIM_t
    for j=i+1:DIM_t
        for u=1:DIM_theta
            for k=1:length(P_rot)
                P_temp = P_rot{k};
                P_topleft = P_temp(1,1);
                if P_topleft == 0
                    P_topright = P_temp(1,2:end);
                    P_bottomleft = P_temp(2:end,1);
                    P_bottomright = P_temp(2:end,2:end);
                    Q_temp = zeros(DIM_z_HOMO, DIM_z_HOMO);
                    Q_temp(tr2zIndex{i+1}(2:end), tthetar2zIndex{j,u}) = P_bottomright;
                    Q_temp(1+i, tthetar2zIndex{j,u}) = P_topright;
                    Q_temp(tthetar2zIndex{j,u}, 1+i) = P_bottomleft;
                    Q_rot{end+1} = Q_temp;
                end
            end
        end
    end
end
% fourth, create lifted version by multiplying t_i t_j and theta_u theta_v
for i=1:DIM_t
    for j=i+1:DIM_t
        for u=1:DIM_theta
            for v=u+1:DIM_theta
                for k=1:length(P_rot)
                    P_temp = P_rot{k};
                    P_topleft = P_temp(1,1);
                    P_topright = P_temp(1,2:end);
                    P_bottomleft = P_temp(2:end,1);
                    if P_topleft==0 && norm(P_topright)==0 && norm(P_bottomleft)==0
                        P_bottomright = P_temp(2:end,2:end);
                        Q_temp = zeros(DIM_z_HOMO, DIM_z_HOMO);
                        Q_temp(tthetar2zIndex{i,u}, tthetar2zIndex{j,v})=P_bottomright;
                        Q_rot{end+1} = Q_temp;
                    end
                end
            end
        end
    end
end
                
%% construct lifted sphere constraints Q_sphere
% first construct unlifted sphere constraints
P_sphere = [-1, zeros(1,3);zeros(3,1), eye(3)];
% then build lifted constraints
Q_sphere = {};
% first, lift by multiplying elements in r_homo
% create appropriate indexing from t_homo * r_i to z
rt2zIndex = {};
for i=1:1+DIM_r
    if i==1 % 1
        rt2zIndex{i} = 1:1+DIM_t;
    else % r_i
        rt2zIndex{i} = [DIM_t+i,...
            (1+DIM_t+DIM_r+DIM_theta+i-1):DIM_r:...
            (1+DIM_t+DIM_r+DIM_theta+i-1+2*DIM_r)];
    end
end
for i=1:1+DIM_r
    for j=i:1+DIM_r
        temp = zeros(DIM_z_HOMO, DIM_z_HOMO);
        temp(rt2zIndex{i}, rt2zIndex{j}) = P_sphere;
        Q_sphere{end+1} = temp;
    end
end
% second, lift by multiplying elements from r_homo and theta_homo
rthetat2zIndex = {};
for i=1:DIM_r
    for j=1:DIM_theta
        temp = (1+DIM_t+DIM_r+DIM_theta+DIM_tr+(j-1)*DIM_tr+i):DIM_r:...
            (1+DIM_t+DIM_r+DIM_theta+DIM_tr+(j-1)*DIM_tr+i+2*DIM_r);
        rthetat2zIndex{i,j} = temp;
    end
end
for i=1:DIM_r
    for j=1:DIM_theta
        P_temp = P_sphere;
        P_topleft = P_temp(1,1);
        P_topright = P_temp(1,2:end);
        P_bottomleft = P_temp(2:end,1);
        P_bottomright = P_temp(2:end,2:end);
        Q_temp = zeros(DIM_z_HOMO, DIM_z_HOMO);
        Q_temp(2:1+DIM_t, rthetat2zIndex{i,j}) = P_bottomright;
        Q_temp(1,rthetat2zIndex{i,j}) = P_topright;
        Q_temp(rthetat2zIndex{i,j},1) = P_bottomleft;
        Q_temp(1+DIM_t+i, 1+DIM_t+DIM_r+j) = P_topleft;
        Q_sphere{end+1} = Q_temp;
    end
end

%% construct binary constraints on theta
% first create unlifted binary constraints
P_binary = [-1, 0; 0, 1];
Q_binary = {};
for i=1:DIM_theta
    temp = zeros(DIM_z_HOMO, DIM_z_HOMO);
    temp([1,1+DIM_t+DIM_r+i], [1,1+DIM_t+DIM_r+i]) = P_binary;
    Q_binary{end+1} = temp;
end
% then create lifted binary constraints by multiplying t_i and r_j
% create appropriate indexing
trtheta2zIndex = {};
for i=1:DIM_t
    for j=1:DIM_r
        for k=1:DIM_theta
            temp = [1+DIM_t+DIM_r+DIM_theta+(i-1)*DIM_r+j,...
                1+DIM_t+DIM_r+DIM_theta+DIM_tr+(k-1)*DIM_tr+(i-1)*DIM_r+j];
            trtheta2zIndex{i,j,k} = temp;
        end
    end
end
for i=1:DIM_t
    for j=1:DIM_r
        for k=1:DIM_theta
            temp = zeros(DIM_z_HOMO, DIM_z_HOMO);
            temp([1,1+DIM_t+DIM_r+k], trtheta2zIndex{i,j,k}) = P_binary;
            Q_binary{end+1} = temp;
        end
    end
end

%% construct auxiliary constraints on t \otimes r, X^1
Q_aux_1 = {};
% create appropriate indexing from vec(X^1) to z
X12zIndex = [1, (1+DIM_t+1):(1+DIM_t+DIM_r)];
for i=1:DIM_t
    temp = [1+i,...
        1+DIM_t+DIM_r+DIM_theta+(i-1)*DIM_r+1:...
        1+DIM_t+DIM_r+DIM_theta+(i-1)*DIM_r+DIM_r];
    X12zIndex = horzcat(X12zIndex, temp);
end
for i=1:DIM_r+1
    for i_=i+1:DIM_r+1
        for j=1:DIM_t+1
            for j_=j+1:DIM_t+1
                e_ii_ = zeros(DIM_r+1, DIM_r+1);
                e_ii_(i,i_) = 1;
                e_jj_ = zeros(DIM_t+1, DIM_t+1);
                e_jj_(j,j_) = 1;
                P = kron((e_jj_ - e_jj_'), e_ii_);
                temp = zeros(DIM_z_HOMO, DIM_z_HOMO);
                temp(X12zIndex, X12zIndex) = P;
                Q_aux_1{end+1} = temp;
            end
        end
    end
end

%% construct auxiliary constraints on \theta \otimes t \otimes r, X^2
Q_aux_2 = {};
% create appropriate indexing from vec(X^2) to z
X22zIndex = [1, 1+DIM_t+DIM_r+1:1+DIM_t+DIM_r+DIM_theta];
for i=1:DIM_t
    for j=1:DIM_r
        temp = [1+DIM_t+DIM_r+DIM_theta+(i-1)*DIM_r+j,...
            1+DIM_t+DIM_r+DIM_theta+DIM_tr+(i-1)*DIM_r+j:DIM_tr:...
            1+DIM_t+DIM_r+DIM_theta+DIM_tr+(i-1)*DIM_r+j+(DIM_theta-1)*DIM_tr];
        X22zIndex = horzcat(X22zIndex, temp);
    end
end
for i=1:1
    for i_=i+1:DIM_theta+1
        for j=1
            for j_=j+1:DIM_tr+1
                e_ii_ = zeros(DIM_theta+1, DIM_theta+1);
                e_ii_(i,i_) = 1;
                e_jj_ = zeros(DIM_tr+1, DIM_tr+1);
                e_jj_(j,j_) = 1;
                P = kron((e_jj_ - e_jj_'), e_ii_);
                temp = zeros(DIM_z_HOMO, DIM_z_HOMO);
                temp(X22zIndex, X22zIndex) = P;
                Q_aux_2{end+1} = temp;
            end
        end
    end
end
%% construct combined tr constraints
Q_comb_tr = {};
for j=1:DIM_theta
    for k=1:length(P_rot)
        P_temp = P_rot{k};
        if P_temp(1,1) == -1
            P_bottomright = P_temp(2:end,2:end);
            Q_temp = zeros(DIM_z_HOMO, DIM_z_HOMO);
            for i=1:3
                temp = zeros(DIM_z_HOMO, DIM_z_HOMO);
                temp(tr2zIndex{i+1}(2:end), tthetar2zIndex{i,j}) = P_bottomright;
                Q_temp = Q_temp + temp;
            end
            Q_temp(1, 1+DIM_t+DIM_r+j) = -1;
            Q_comb_tr{end+1} = Q_temp;
        end
    end
end
for u=1:DIM_theta
    for v=u+1:DIM_theta
        for k=1:length(P_rot)
            P_temp = P_rot{k};
            if P_temp(1,1) == -1
                P_bottomright = P_temp(2:end,2:end);
                Q_temp = zeros(DIM_z_HOMO, DIM_z_HOMO);
                for i=1:3
                    temp = zeros(DIM_z_HOMO, DIM_z_HOMO);
                    temp(tthetar2zIndex{i,u}, tthetar2zIndex{i,v}) = P_bottomright;
                    Q_temp = Q_temp + temp;
                end
                Q_temp(1+DIM_t+DIM_r+u, 1+DIM_t+DIM_r+v) = -1;
                Q_comb_tr{end+1} = Q_temp;
            end
        end
    end
end        
            
%% construct cost matrix Q_cost
Q_cost_1 = {};
Q_cost_2 = {};
Q_cost_3 = {};
Q_cost_4 = {};
Q_cost = {};
c = 1;
tr_cost_zIndex = (1+DIM_t+DIM_r+DIM_theta+1):DIM_r:...
    (1+DIM_t+DIM_r+DIM_theta+1+2*DIM_r);
for i=2:DIM_r
    temp = (1+DIM_t+DIM_r+DIM_theta+i):DIM_r:...
        (1+DIM_t+DIM_r+DIM_theta+i+2*DIM_r);
    tr_cost_zIndex = horzcat(tr_cost_zIndex, temp);
end
rt_cost_zIndex = (1+DIM_t+DIM_r+DIM_theta+1):...
    (1+DIM_t+DIM_r+DIM_theta+DIM_tr);
for i=1:NUM_POINTS
    f2 = bearing2(:,i);
    f1 = bearing1(:,i);
    G = kron(f2', hatmap(f1));
    A = kron(G', G);
    temp_1 = zeros(DIM_z_HOMO, DIM_z_HOMO);
    temp_2 = zeros(DIM_z_HOMO, DIM_z_HOMO);
    temp_3 = zeros(DIM_z_HOMO, DIM_z_HOMO);
    temp_4 = zeros(DIM_z_HOMO, DIM_z_HOMO);
    temp_1(tr_cost_zIndex, rt_cost_zIndex) = A/2;
    theta_i_rt_cost_zIndex = (1+DIM_t+DIM_r+DIM_theta+DIM_tr+(i-1)*DIM_tr+1):...
    (1+DIM_t+DIM_r+DIM_theta+DIM_tr+(i-1)*DIM_tr+DIM_tr);
    temp_2(tr_cost_zIndex, theta_i_rt_cost_zIndex) = A/2;
    temp_3(1, 1+DIM_t+DIM_r+i) = -c^2/2;
    temp_4(1,1) = c^2/2;
    Q_cost_1{end+1} = temp_1;
    Q_cost_2{end+1} = temp_2;
    Q_cost_3{end+1} = temp_3;
    Q_cost_4{end+1} = temp_4;
    temp = temp_1 + temp_2 + temp_3 + temp_4;
    Q_cost{end+1} = temp;
end
Q_cost_total = zeros(DIM_z_HOMO, DIM_z_HOMO);
for i=1:length(Q_cost)
    Q_cost_total = Q_cost_total + Q_cost{i};
end

%% print a few information
disp('Size of decision matrix:')
disp(DIM_z_HOMO)
disp('Number of lifted constraints:')
disp(length(Q_rot))
disp(length(Q_sphere))
disp(length(Q_binary))
disp(length(Q_aux_1))
disp(length(Q_aux_2))
disp(length(Q_comb_tr))
disp(length(Q_cost))

%% DO sanity check using ground truth
errorMatrix = {};
tolerance = 1e-12;
r_gt = reshape(R_gt, [9,1]);
t_gt = t_gt/norm(t_gt);
theta_gt = ones(NUM_POINTS,1); % accept all the correspondences
z_gt = [1;t_gt;r_gt;theta_gt;kron(t_gt,r_gt);kron(theta_gt, kron(t_gt,r_gt))];
ZHomo_gt = z_gt * z_gt';
cost_gt = trace(Q_cost_total * ZHomo_gt);
disp('-----------------Sanity Check--------------')
if abs(cost_gt) > tolerance
    disp(cost_gt)
    error('GT cost is nonzero.')
end
for i=1:length(Q_rot)
    if abs(trace(Q_rot{i} * ZHomo_gt)) > tolerance
        disp(i)
        disp(trace(Q_rot{i} * ZHomo_gt))
        assignin('base','Q_rot',Q_rot{i});
        error('Rotation costraint not met.')
    end
end

for i=1:length(Q_sphere)
    if abs(trace(Q_sphere{i} * ZHomo_gt)) > tolerance
        disp(i)
        disp(trace(Q_sphere{i} * ZHomo_gt))
        assignin('base','Q_sphere', Q_sphere{i});
        error('Translation costraint not met.')
    end
end

for i=1:length(Q_binary)
    if abs(trace(Q_binary{i} * ZHomo_gt)) > tolerance
        disp(trace(Q_binary{i} * ZHomo_gt))
        error('Binary costraint not met.')
    end
end

for i=1:length(Q_aux_1)
    if abs(trace(Q_aux_1{i} * ZHomo_gt)) > tolerance
        disp(trace(Q_aux_1{i} * ZHomo_gt))
        error('AUX 1 costraint not met.')
    end
end

for i=1:length(Q_aux_2)
    if abs(trace(Q_aux_2{i} * ZHomo_gt)) > tolerance
        disp(trace(Q_aux_2{i} * ZHomo_gt))
        error('AUX 2 costraint not met.')
    end
end

for i=1:length(Q_comb_tr)
    if abs(trace(Q_comb_tr{i} * ZHomo_gt)) > tolerance
        disp(trace(Q_comb_tr{i} * ZHomo_gt))
        error('Combined tr costraint not met.')
    end
end
disp('-----------------END Sanity Check--------------')

%% begin CVX solve
if quiet
    disp('set quiet')
    cvx_begin quiet sdp
else
    cvx_begin sdp
end
cvx_precision best
variable ZHomo(DIM_z_HOMO, DIM_z_HOMO) symmetric
minimize( trace(Q_cost_total * ZHomo) )
subject to
ZHomo >= 0
ZHomo(1,1) == 1
disp('---ADD ROT CONSTRAINTS---')
for i=1:length(Q_rot)
    trace(Q_rot{i} * ZHomo) == 0;
end
disp('---ADD TRANS CONSTRAINTS---')
for i=1:length(Q_sphere)
    trace(Q_sphere{i} * ZHomo) == 0;
end
disp('---ADD BINARY CONSTRAINTS---')
for i=1:length(Q_binary)
    trace(Q_binary{i} * ZHomo) == 0;
end
disp('---ADD AUX 1 CONSTRAINTS---')
for i=1:length(Q_aux_1)
    trace(Q_aux_1{i} * ZHomo) == 0;
end
disp('---ADD AUX 2 CONSTRAINTS---')
for i=1:length(Q_aux_2)
    trace(Q_aux_2{i} * ZHomo) == 0;
end
disp('---ADD COMB TR CONSTRAINTS---')
for i=1:length(Q_comb_tr)
    trace(Q_comb_tr{i} * ZHomo) == 0;
end
cvx_end

%% get results
cvxStatus = cvx_status;
pStarRelax = cvx_optval;
ZStar = ZHomo;

figure
spy(ZStar)
figure
bar(eig(full(ZStar)))

end
    
    
            
        
        

            