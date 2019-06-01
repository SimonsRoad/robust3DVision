function [cvxStatus,pStarRelax,WStar] = ...
    qcqprelaxation_l1norm_binary_variable(bearing1,bearing2,R_gt,t_gt,quiet)
% This function writes the L1 cost function using binary variables.
% See discussion.pdf for corresponding math derivations.
% To date, this relaxation is not tight.
% Last update: 12/03/2018.
%% useful constants
NUM_POINTS = size(bearing1, 2);
DIM_r = 9;
DIM_t = 3;
DIM_b = NUM_POINTS;
DIM_tr = DIM_r*DIM_t;
DIM_bt = DIM_b*DIM_t;
DIM_br = DIM_b*DIM_r;
DIM_w_HOMO = 1+DIM_t+DIM_r+DIM_b+DIM_tr+DIM_bt+DIM_br;

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
% w=[1,t,r,b,txr,bxt,bxr]  r_c = [1,t,b]
dim_r_c = 1+DIM_t+DIM_b;
rc2wIndex = {};
for i=1:dim_r_c
    if i==1 % 1
        rc2wIndex{i} = [1,(1+DIM_t+1):(1+DIM_t+DIM_r)];
    elseif i<=1+DIM_t % t
        rc2wIndex{i} = [i,...
            (1+DIM_t+DIM_r+DIM_b+(i-2)*DIM_r+1):...
            (1+DIM_t+DIM_r+DIM_b+(i-2)*DIM_r+DIM_r)];
    else % b
        rc2wIndex{i} = [i+DIM_r,...
            (1+DIM_t+DIM_r+DIM_b+DIM_tr+DIM_bt+(i-5)*DIM_r+1):...
            (1+DIM_t+DIM_r+DIM_b+DIM_tr+DIM_bt+(i-5)*DIM_r+DIM_r)];
    end
end
Q_rot = {};
for i=1:dim_r_c
    for j=i:dim_r_c
        for k=1:length(P_rot)
            temp = zeros(DIM_w_HOMO, DIM_w_HOMO);
            temp(rc2wIndex{i}, rc2wIndex{j}) = P_rot{k};
            Q_rot{end+1} = temp;
        end
    end
end
%% construct lifted sphere constraints Q_sphere
% first construct unlifted sphere constraints
P_sphere = [-1, zeros(1,3);zeros(3,1), eye(3)];
% then build lifted constraints
% w=[1,t,r,b,txr,bxt,bxr]  t_c = [1,r,b]
dim_t_c = 1+DIM_r+DIM_b;
tc2wIndex = {};
for i=1:dim_t_c
    if i==1 % 1
        tc2wIndex{i} = 1:(1+DIM_t);
    elseif i <= 1+DIM_r % r
        tc2wIndex{i} = [i+DIM_t,...
            (1+DIM_t+DIM_r+DIM_b+i-1):DIM_r:...
            (1+DIM_t+DIM_r+DIM_b+i-1+2*DIM_r)];
    else % b
        tc2wIndex{i} = [i+DIM_t,...
            (1+DIM_t+DIM_r+DIM_b+DIM_tr+(i-11)*DIM_t+1):...
            (1+DIM_t+DIM_r+DIM_b+DIM_tr+(i-11)*DIM_t+DIM_t)];
    end
end
Q_sphere = {};
for i=1:dim_t_c
    for j=i:dim_t_c
        temp = zeros(DIM_w_HOMO, DIM_w_HOMO);
        temp(tc2wIndex{i}, tc2wIndex{j}) = P_sphere;
        Q_sphere{end+1} = temp;
    end
end
%% construct lifted binary constraints Q_binary
% first construct unlifted binary constraints
P_binary = {};
for i=1:NUM_POINTS
    temp = zeros(NUM_POINTS, NUM_POINTS);
    temp(i,i) = 1;
    P_binary{end+1} = [-1, zeros(1,NUM_POINTS); zeros(NUM_POINTS,1), temp];
end
% then construct lifted binary constraints
% w=[1,t,r,b,txr,bxt,bxr]  b_c = [1,t,r]
dim_b_c = 1+DIM_t+DIM_r;
bc2wIndex = {};
for i=1:dim_b_c
    if i==1 % 1
        bc2wIndex{i} = [1, 1+DIM_t+DIM_r+1:...
            1+DIM_t+DIM_r+DIM_b];
    elseif i<=1+DIM_t % t
        bc2wIndex{i} = [i,...
            (1+DIM_t+DIM_r+DIM_b+DIM_tr+(i-1)):DIM_t:...
            (1+DIM_t+DIM_r+DIM_b+DIM_tr+(i-1)+(DIM_b-1)*DIM_t)];
    else % r
        bc2wIndex{i} = [i,...
            (1+DIM_t+DIM_r+DIM_b+DIM_tr+DIM_bt+(i-4)):DIM_r:...
            (1+DIM_t+DIM_r+DIM_b+DIM_tr+DIM_bt+(i-4)+(DIM_b-1)*DIM_r)];
    end
end
Q_binary = {};
for i=1:dim_b_c
    for j=i:dim_b_c
        for u=1:length(P_binary)
            temp = zeros(DIM_w_HOMO, DIM_w_HOMO);
            temp(bc2wIndex{i}, bc2wIndex{j}) = P_binary{u};
            Q_binary{end+1} = temp;
        end
    end
end
%% constuct auxiliary constraints on X_k
Q_aux = {};
% first do X_rt X_rt = [1,t^T; r, rt^T]
% w=[1,t,r,b,txr,bxt,bxr]
% X_rt size: DIM_r+1 by DIM_t+1
XrtIndex = [1, (1+DIM_t+1):(1+DIM_t+DIM_r)];
for i=1:DIM_t
    temp = [1+i, (1+DIM_t+DIM_r+DIM_b+(i-1)*DIM_r+1):...
        (1+DIM_t+DIM_r+DIM_b+(i-1)*DIM_r+DIM_r)];
    XrtIndex = horzcat(XrtIndex, temp);
end
for i=1:(DIM_r+1)
    for i_=i:(DIM_r+1)
        for j=1:(DIM_t+1)
            for j_=j:(DIM_t+1)
                e_ii_ = zeros(DIM_r+1, DIM_r+1);
                e_ii_(i,i_) = 1;
                e_jj_ = zeros(DIM_t+1, DIM_t+1);
                e_jj_(j,j_) = 1;
                P = kron((e_jj_ - e_jj_'), e_ii_);
                temp = zeros(DIM_w_HOMO, DIM_w_HOMO);
                temp(XrtIndex, XrtIndex) = P;
                Q_aux{end+1} = temp;
            end
        end
    end
end
% then do X_rb X_rb = [1, b^T; r, rb^T]
% w=[1,t,r,b,txr,bxt,bxr]
% X_rb size: DIM_r+1 by DIM_b+1
XrbIndex = [1, (1+DIM_t+1):(1+DIM_t+DIM_r)];
for i=1:DIM_b
    temp = [1+DIM_t+DIM_r+i,...
        1+DIM_t+DIM_r+DIM_b+DIM_tr+DIM_bt+(i-1)*DIM_r+1:...
        1+DIM_t+DIM_r+DIM_b+DIM_tr+DIM_bt+(i-1)*DIM_r+DIM_r];
    XrbIndex = horzcat(XrbIndex, temp);
end
for i=1:(DIM_r+1)
    for i_=i:(DIM_r+1)
        for j=1:(DIM_b+1)
            for j_=j:(DIM_b+1)
                e_ii_ = zeros(DIM_r+1, DIM_r+1);
                e_ii_(i,i_) = 1;
                e_jj_ = zeros(DIM_b+1, DIM_b+1);
                e_jj_(j,j_) = 1;
                P = kron((e_jj_ - e_jj_'), e_ii_);
                temp = zeros(DIM_w_HOMO, DIM_w_HOMO);
                temp(XrbIndex, XrbIndex) = P;
                Q_aux{end+1} = temp;
            end
        end
    end
end
% then do X_tb X_tb = [1, b^T; t, rb^T]
% w=[1,t,r,b,txr,bxt,bxr]
% X_tb size: DIM_t+1 by DIM_b+1
XtbIndex = [1, (1+1):(1+DIM_t)];
for i=1:DIM_b
    temp = [1+DIM_t+DIM_r+i,...
        (1+DIM_t+DIM_r+DIM_b+DIM_tr+(i-1)*DIM_t+1):...
        (1+DIM_t+DIM_r+DIM_b+DIM_tr+(i-1)*DIM_t+DIM_t)];
    XtbIndex = horzcat(XtbIndex, temp);
end
for i=1:(DIM_t+1)
    for i_=i:(DIM_t+1)
        for j=1:(DIM_b+1)
            for j_=j:(DIM_b+1)
                e_ii_ = zeros(DIM_t+1, DIM_t+1);
                e_ii_(i,i_) = 1;
                e_jj_ = zeros(DIM_b+1, DIM_b+1);
                e_jj_(j,j_) = 1;
                P = kron((e_jj_ - e_jj_'), e_ii_);
                temp = zeros(DIM_w_HOMO, DIM_w_HOMO);
                temp(XtbIndex, XtbIndex) = P;
                Q_aux{end+1} = temp;
            end
        end
    end
end
%% construct cost matrix Q_A
Q_A = {};
for i=1:NUM_POINTS
    f2 = bearing2(:,i);
    f1 = bearing1(:,i);
    A = kron(f2', hatmap(f1));
    temp = zeros(DIM_w_HOMO, DIM_w_HOMO);
    tIndex = 2:(1+DIM_t);
    birIndex = (1+DIM_t+DIM_r+DIM_b+DIM_tr+DIM_bt+(i-1)*DIM_r+1):...
        (1+DIM_t+DIM_r+DIM_b+DIM_tr+DIM_bt+(i-1)*DIM_r+DIM_r);
    temp(tIndex, birIndex) = A;
    Q_A{end+1} = temp;
end
%% construct total cost Q_cost
Q_cost = zeros(DIM_w_HOMO, DIM_w_HOMO);
for i=1:length(Q_A)
    Q_cost = Q_cost + Q_A{i};
end


%% print a few information
disp('Size of decision matrix:')
disp(DIM_w_HOMO)
disp('Number of lifted constraints:')
disp('-----Q_rot------')
disp(length(Q_rot))
disp('-----Q_sphere------')
disp(length(Q_sphere))
disp('-----Q_binary------')
disp(length(Q_binary))
disp('-----Q_aux------')
disp(length(Q_aux))
disp('-----Q_A------')
disp(length(Q_A))

%% sanity check
errorMatrix = {};
tolerance = 1e-12;
r_gt = reshape(R_gt, [9,1]);
t_gt = t_gt/norm(t_gt);
b_gt = ones(NUM_POINTS,1); % since noise level is 0, b does not matter
w_gt = [1;t_gt;r_gt;b_gt;...
    kron(t_gt,r_gt);kron(b_gt,t_gt);kron(b_gt,r_gt)];
WHomo_gt = w_gt * w_gt';
cost_gt = trace(Q_cost * WHomo_gt);
disp('-----------------Sanity Check--------------')
if abs(cost_gt) > tolerance
    disp(cost_gt)
    error('GT cost is nonzero.')
end
for i=1:length(Q_rot)
    if abs(trace(Q_rot{i} * WHomo_gt)) > tolerance
        disp(i)
        disp(trace(Q_rot{i} * WHomo_gt))
        assignin('base','Q_rot',Q_rot{i});
        error('Rotation costraint not met.')
    end
end

for i=1:length(Q_sphere)
    if abs(trace(Q_sphere{i} * WHomo_gt)) > tolerance
        disp(i)
        disp(trace(Q_sphere{i} * WHomo_gt))
        assignin('base','Q_sphere', Q_sphere{i});
        error('Translation costraint not met.')
    end
end

for i=1:length(Q_binary)
    if abs(trace(Q_binary{i} * WHomo_gt)) > tolerance
        disp(trace(Q_binary{i} * WHomo_gt))
        error('Binary costraint not met.')
    end
end

for i=1:length(Q_aux)
    if abs(trace(Q_aux{i} * WHomo_gt)) > tolerance
        disp(trace(Q_aux{i} * WHomo_gt))
        error('AUX costraint not met.')
    end
end

for i=1:length(Q_A)
    if  trace(Q_A{i} * WHomo_gt) < -tolerance
        disp(trace(Q_A{i} * WHomo_gt))
        error('Cost positive costraint not met.')
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
variable WHomo(DIM_w_HOMO, DIM_w_HOMO) symmetric
minimize( trace(Q_cost*WHomo) + norm_nuc(WHomo) )
subject to
WHomo >= 0
WHomo(1,1) == 1
disp('-------ADD ROT CONSTRAINT-----')
for i=1:length(Q_rot)
    if rem(i,100)==1
        disp(i)
    end
    trace(Q_rot{i} * WHomo) == 0;
end

disp('-------ADD TRANS CONSTRAINT-----')
for i=1:length(Q_sphere)
    if rem(i,100)==1
        disp(i)
    end
    trace(Q_sphere{i} * WHomo) == 0;
end

disp('-------ADD BINARY CONSTRAINT-----')
for i=1:length(Q_binary)
    if rem(i,100)==1
        disp(i)
    end
    trace(Q_binary{i} * WHomo) == 0;
end

disp('-------ADD AUX CONSTRAINT-----')
for i=1:length(Q_aux)
    if rem(i,100)==1
        disp(i)
    end
    trace(Q_aux{i} * WHomo) == 0;
end

disp('-------ADD COST CONSTRAINT-----')
for i=1:length(Q_A)
    if rem(i,100)==1
        disp(i)
    end
    trace(Q_A{i} * WHomo) >= 0;
end
cvx_end

%% get results
cvxStatus = cvx_status;
pStarRelax = cvx_optval;
WStar = WHomo;

figure
spy(WStar)
figure
bar(eig(full(WStar)))

end
    
    
            
        
        

            
