function [cvxStatus,pStarRelax,ZStar,Q_aux,Q_sphere,Q_rot,Q_A] = ...
    qcqprelaxation_l1norm_full_monomial(bearing1,bearing2,quiet)

%% useful constants
NUM_POINTS = size(bearing1, 2);
DIM_r = 9;
DIM_t = 3;
DIM_z_HOMO = 1+DIM_t+DIM_r+DIM_t*DIM_r+DIM_t^2+DIM_r^2;

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
% z_homo=[1,t,r,txr,txt,rxr]  v = [1,t,r]
dim_v = 1+DIM_t+DIM_r;
vr2zIndex = {};
for i=1:dim_v
    if i==1 % 1
        vr2zIndex{i} = [1,(1+DIM_t+1):(1+DIM_t+DIM_r)];
    elseif i<=1+DIM_t % t
        vr2zIndex{i} = [i,...
            (1+DIM_t+DIM_r+(i-2)*DIM_r+1):...
            (1+DIM_t+DIM_r+(i-2)*DIM_r+DIM_r)];
    else % r
        vr2zIndex{i} = [i,...
            (1+DIM_t+DIM_r+DIM_t*DIM_r+DIM_t^2+(i-5)*DIM_r+1):...
            (1+DIM_t+DIM_r+DIM_t*DIM_r+DIM_t^2+(i-5)*DIM_r+DIM_r)];
    end
end
Q_rot = {};
for i=1:dim_v
    for j=i:dim_v
        for k=1:length(P_rot)
            temp = zeros(DIM_z_HOMO, DIM_z_HOMO);
            temp(vr2zIndex{i}, vr2zIndex{j}) = P_rot{k};
            Q_rot{end+1} = temp;
        end
    end
end
%% construct lifted sphere constraints Q_sphere
% first construct unlifted sphere constraints
P_sphere = [-1, zeros(1,3);zeros(3,1), eye(3)];
% then build lifted constraints
% z_homo=[1,t,r,txr,txt,rxr]  v = [1,t,r]
vt2zIndex = {};
for i=1:dim_v
    if i==1 % 1
        vt2zIndex{i} = 1:(1+DIM_t);
    elseif i <= 1+DIM_t % t
        vt2zIndex{i} = [i,...
            (1+DIM_t+DIM_r+DIM_t*DIM_r+(i-2)*DIM_t+1):...
            (1+DIM_t+DIM_r+DIM_t*DIM_r+(i-2)*DIM_t+DIM_t)];
    else % r
        vt2zIndex{i} = [i,...
            (1+DIM_t+DIM_r+(i-4)):DIM_r:...
            (1+DIM_t+DIM_r+(i-4)+2*DIM_r)];
    end
end
Q_sphere = {};
for i=1:dim_v
    for j=i:dim_v
        temp = zeros(DIM_z_HOMO, DIM_z_HOMO);
        temp(vt2zIndex{i}, vt2zIndex{j}) = P_sphere;
        Q_sphere{end+1} = temp;
    end
end

%% constuct auxiliary constraints on X_u
Q_aux = {};

% first do X^tt X^tt = [1,t^T; t, tt^T]
% X^tt size: DIM_t+1 by DIM_t+1
XttIndex = 1:(1+DIM_t);
for i=1:DIM_t
    temp = (1+DIM_t+DIM_r+DIM_t*DIM_r+(i-1)*DIM_t+1):...
           (1+DIM_t+DIM_r+DIM_t*DIM_r+(i-1)*DIM_t+DIM_t);
    XttIndex = horzcat(XttIndex, temp);
end
real_t_index = [1:4,2,5:7,3,8:10,4,11:13];
for i=1:(DIM_t+1)
    for i_=(i+1):(DIM_t+1)
        for j=1:(DIM_t+1)
            for j_=(j+1):(DIM_t+1)
                e_ii_ = zeros(DIM_t+1, DIM_t+1);
                e_ii_(i,i_) = 1;
                e_jj_ = zeros(DIM_t+1, DIM_t+1);
                e_jj_(j,j_) = 1;
                P = kron((e_jj_ - e_jj_'), e_ii_);
                % extra care needed here
                [rows, cols] = find(P~=0);
                if (length(rows) ~= 2)
                    disp(length(rows))
                end
                P_real = zeros(length(XttIndex), length(XttIndex));
                for idx=1:length(rows)
                    row = rows(idx);
                    col = cols(idx);
                    real_row = real_t_index(row);
                    real_col = real_t_index(col);
                    P_real(real_row, real_col) = P(row, col);
                end 
                temp = zeros(DIM_z_HOMO, DIM_z_HOMO);
                temp(XttIndex, XttIndex) = P_real;
                Q_aux{end+1} = temp;
            end
        end
    end
end
% then do X^rr X^rr = [1, r^T; r, rr^T]
% X_rr size: DIM_r+1 by DIM_r+1
XrrIndex = [1, (1+DIM_t+1):(1+DIM_t+DIM_r)];
for i=1:DIM_r
    temp = (1+DIM_t+DIM_r+DIM_t*DIM_r+DIM_t^2+(i-1)*DIM_r+1):...
           (1+DIM_t+DIM_r+DIM_t*DIM_r+DIM_t^2+(i-1)*DIM_r+DIM_r);
    XrrIndex = horzcat(XrrIndex, temp);
end
real_r_index = [1:10,2,11:19,3,20:28,4,29:37,5,38:46,6,47:55,7,56:64,8,65:73,9,74:82,10,83:91];
for i=1:(DIM_r+1)
    for i_=(i+1):(DIM_r+1)
        for j=1:(DIM_r+1)
            for j_=(j+1):(DIM_r+1)
                e_ii_ = zeros(DIM_r+1, DIM_r+1);
                e_ii_(i,i_) = 1;
                e_jj_ = zeros(DIM_r+1, DIM_r+1);
                e_jj_(j,j_) = 1;
                P = kron((e_jj_ - e_jj_'), e_ii_);
                % extra care needed here
                [rows, cols] = find(P~=0);
                if (length(rows) ~= 2)
                    disp(length(rows))
                end
                P_real = zeros(length(XrrIndex), length(XrrIndex));
                for idx=1:length(rows)
                    row = rows(idx);
                    col = cols(idx);
                    real_row = real_r_index(row);
                    real_col = real_r_index(col);
                    P_real(real_row, real_col) = P(row, col);
                end 
                temp = zeros(DIM_z_HOMO, DIM_z_HOMO);
                temp(XrrIndex, XrrIndex) = P_real;
                Q_aux{end+1} = temp;
            end
        end
    end
end
% then do X^rt X^rt = [1, t^T; r, rt^T]
% X^rt size: DIM_r+1 by DIM_t+1
XrtIndex = [1,(1+DIM_t+1):(1+DIM_t+DIM_r)];
for i=1:DIM_t
    temp = [1+i,...
        (1+DIM_t+DIM_r+(i-1)*DIM_r+1):...
        (1+DIM_t+DIM_r+(i-1)*DIM_r+DIM_r)];
    XrtIndex = horzcat(XrtIndex, temp);
end
for i=1:(DIM_r+1)
    for i_=(1+i):(DIM_r+1)
        for j=1:(DIM_t+1)
            for j_=(1+j):(DIM_t+1)
                e_ii_ = zeros(DIM_r+1, DIM_r+1);
                e_ii_(i,i_) = 1;
                e_jj_ = zeros(DIM_t+1, DIM_t+1);
                e_jj_(j,j_) = 1;
                P = kron((e_jj_ - e_jj_'), e_ii_);
                temp = zeros(DIM_z_HOMO, DIM_z_HOMO);
                temp(XrtIndex, XrtIndex) = P;
                Q_aux{end+1} = temp;
            end
        end
    end
end
for i=1:DIM_t
    temp = zeros(DIM_z_HOMO, DIM_z_HOMO);
    temp(i+1, i+1) = 1;
    temp(1, 1+DIM_t+DIM_r+DIM_t*DIM_r+1+(i-1)*4) = -1;
    Q_aux{end+1} = temp;
end
for i=1:DIM_r
    temp = zeros(DIM_z_HOMO, DIM_z_HOMO);
    temp(1+DIM_t+i, 1+DIM_t+i) = 1;
    temp(1, 1+DIM_t+DIM_r+DIM_t*DIM_r+DIM_t^2+1+(i-1)*10) = -1;
    Q_aux{end+1} = temp;
end
%% construct cost matrix Q_A
Q_A = {};
tIndex = 2:(1+DIM_t);
rIndex = (1+DIM_t+1):(1+DIM_t+DIM_r);
Q_A = {};
for i=1:NUM_POINTS
    f1 = bearing1(:, i);
    f2 = bearing2(:, i);
    A = kron(f2', hatmap(f1));
    Q_i = zeros(DIM_z_HOMO, DIM_z_HOMO);
    Q_i(tIndex, rIndex) = A/2;
    Q_i(rIndex, tIndex) = A'/2;
    Q_A{end+1} = Q_i;
end

%% print a few information
disp('Size of decision matrix:')
disp(DIM_z_HOMO)
disp('Number of lifted constraints:')
disp(length(Q_sphere))
disp(length(Q_rot))
disp(length(Q_aux))
disp(length(Q_A))

%% begin CVX solve
if quiet
    disp('set quiet')
    cvx_begin quiet sdp
else
    cvx_begin sdp
end
cvx_precision best
variable ZHomo(DIM_z_HOMO, DIM_z_HOMO) symmetric
variable y(NUM_POINTS)
minimize( sum(y) + norm_nuc(ZHomo(2:4,2:4)) + norm_nuc(ZHomo(5:13,5:13)) ) % use nuclear norm to penalize rank
subject to
ZHomo >= 0
ZHomo(1,1) == 1
ZHomo(41:49,1) == reshape(ZHomo(tIndex,tIndex), [9,1])
ZHomo(50:130,1) == reshape(ZHomo(rIndex,rIndex), [81,1])
for i=1:length(Q_rot)
    trace(Q_rot{i} * ZHomo) == 0;
end
for i=1:length(Q_sphere)
    trace(Q_sphere{i} * ZHomo) == 0;
end
for i=1:length(Q_aux)
    trace(Q_aux{i} * ZHomo) == 0;
end
for i=1:length(Q_A)
    y(i) >= abs(trace(Q_A{i} * ZHomo));
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
    
    
            
        
        

            