function [cvxStatus,pStarRelax,ZStar] = ...
    qcqprelaxation_truncated_ls_outlier(bearing1,bearing2,R_gt,t_gt,sanityCheck,quiet)
% See truncated_LS_update.pdf for the corresponding math derivation
% This is the default truncated LS formuation.
% Last update 12/03/2018
%% useful constants
NUM_POINTS = size(bearing1, 2);
if NUM_POINTS < 7
    error('QCQP needs at least 7 points to work with outliers.')
end
DIM_theta = NUM_POINTS;
DIM_r = 9;
DIM_t = 3;
DIM_tr = DIM_t * DIM_r;
DIM_thetar = DIM_theta * DIM_r;
DIM_thetat = DIM_theta * DIM_t;
DIM_thetatr = DIM_theta * DIM_t * DIM_r;
DIM_z_HOMO = 1+DIM_t+DIM_r+DIM_theta+DIM_tr+DIM_thetat+DIM_thetar+DIM_thetatr;

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
% create appropriate indexing from t_i t_j theta_u theta_v r to z
t_theta_r_zIndex = {};
for i=1:1+DIM_t
    for u=1:1+DIM_theta
        if i==1 && u==1
            temp = [1, (1+DIM_t+1):(1+DIM_t+DIM_r)];    
        elseif i==1 && u>1
            temp = [DIM_t+DIM_r+u,...
                1+DIM_t+DIM_r+DIM_theta+DIM_tr+DIM_thetat+(u-2)*DIM_r+1:...
                1+DIM_t+DIM_r+DIM_theta+DIM_tr+DIM_thetat+(u-2)*DIM_r+DIM_r];    
        elseif i>1 && u==1
            temp = [i,...
                1+DIM_t+DIM_r+DIM_theta+(i-2)*DIM_r+1:...
                1+DIM_t+DIM_r+DIM_theta+(i-2)*DIM_r+DIM_r];
        else
            temp = [1+DIM_t+DIM_r+DIM_theta+DIM_tr+(u-2)*DIM_t+i-1,...
                1+DIM_t+DIM_r+DIM_theta+DIM_tr+DIM_thetat+DIM_thetar+(u-2)*DIM_tr+(i-2)*DIM_r+1:...
                1+DIM_t+DIM_r+DIM_theta+DIM_tr+DIM_thetat+DIM_thetar+(u-2)*DIM_tr+(i-2)*DIM_r+DIM_r];    
        end
        t_theta_r_zIndex{i,u}=temp;
    end
end
for i=1:1+DIM_t
    for j=i:1+DIM_t
        for u=1:1+DIM_theta
            for  v=u:1+DIM_theta
                for k=1:length(P_rot)
                    rowIndex=t_theta_r_zIndex{i,u};
                    colIndex=t_theta_r_zIndex{j,v};
                    temp=zeros(DIM_z_HOMO,DIM_z_HOMO);
                    temp(rowIndex,colIndex)=P_rot{k};
                    Q_rot{end+1}=sparse(temp);
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
% create appropriate indexing from r_i r_j theta_u theta_v t to z
r_theta_t_zIndex = {};
for i=1:1+DIM_r
    for u=1:1+DIM_theta
        if i==1 && u==1
            temp = 1:1+DIM_t;
        elseif i==1 && u>1
            temp=[DIM_t+DIM_r+u,...
                1+DIM_t+DIM_r+DIM_theta+DIM_tr+(u-2)*DIM_t+1:...
                1+DIM_t+DIM_r+DIM_theta+DIM_tr+(u-2)*DIM_t+DIM_t];
        elseif i>1 && u==1
            temp=[DIM_t+i,...
                1+DIM_t+DIM_r+DIM_theta+(i-1):DIM_r:...
                1+DIM_t+DIM_r+DIM_theta+(i-1)+2*DIM_r];    
        else
            temp=[1+DIM_t+DIM_r+DIM_theta+DIM_tr+DIM_thetat+(u-2)*DIM_r+(i-1),...
                1+DIM_t+DIM_r+DIM_theta+DIM_tr+DIM_thetat+DIM_thetar+(u-2)*DIM_tr+(i-1):DIM_r:...
                1+DIM_t+DIM_r+DIM_theta+DIM_tr+DIM_thetat+DIM_thetar+(u-2)*DIM_tr+(i-1)+2*DIM_r];
        end
        r_theta_t_zIndex{i,u} = temp;
    end
end
for i=1:1+DIM_r
    for j=i:1+DIM_r
        for u=1:1+DIM_theta
            for v=u:1+DIM_theta
                rowIndex=r_theta_t_zIndex{i,u};
                colIndex=r_theta_t_zIndex{j,v};
                temp=zeros(DIM_z_HOMO, DIM_z_HOMO);
                temp(rowIndex,colIndex)=P_sphere;
                Q_sphere{end+1}=sparse(temp);
            end
        end
    end
end

%% construct binary constraints on theta
% first create unlifted binary constraints
P_binary = [-1, 0; 0, 1];
% create lifted constraints
Q_binary = {};
% create appropriate indexing from t_i t_j r_u r_v theta_k to z
t_r_theta_zIndex={};
for i=1:1+DIM_t
    for u=1:1+DIM_r
        for k=1:DIM_theta
            if i==1 && u==1
                temp=[1,1+DIM_t+DIM_r+k];
            elseif i==1 && u>1
                temp=[DIM_t+u,...
                    1+DIM_t+DIM_r+DIM_theta+DIM_tr+DIM_thetat+(k-1)*DIM_r+u-1];
            elseif i>1 && u==1
                temp=[i,...
                    1+DIM_t+DIM_r+DIM_theta+DIM_tr+(k-1)*DIM_t+i-1];
            else
                temp=[1+DIM_t+DIM_r+DIM_theta+(i-2)*DIM_r+u-1,...
                    1+DIM_t+DIM_r+DIM_theta+DIM_tr+DIM_thetat+DIM_thetar+(k-1)*DIM_tr+(i-2)*DIM_r+u-1];
            end
            t_r_theta_zIndex{i,u,k}=temp;
        end
    end
end
for i=1:1+DIM_t
    for j=i:1+DIM_t
        for u=1:1+DIM_r
            for v=u:1+DIM_r
                for k=1:DIM_theta
                    rowIndex=t_r_theta_zIndex{i,u,k};
                    colIndex=t_r_theta_zIndex{j,v,k};
                    temp=zeros(DIM_z_HOMO,DIM_z_HOMO);
                    temp(rowIndex,colIndex)=P_binary;
                    Q_binary{end+1}=sparse(temp);
                end
            end
        end
    end
end

%% construct auxiliary constraints degree 1
Q_aux_1 = {};
% first, do X^1_rt
% create appropriate indexing from X^1_rt to z
X_1_rt_zIndex = [1, 1+DIM_t+1:1+DIM_t+DIM_r];
for i=1:DIM_t
    temp=[1+i,...
        1+DIM_t+DIM_r+DIM_theta+(i-1)*DIM_r+1:...
        1+DIM_t+DIM_r+DIM_theta+(i-1)*DIM_r+DIM_r];
    X_1_rt_zIndex=horzcat(X_1_rt_zIndex,temp);
end
for i=1:1+DIM_r
    for i_=i+1:1+DIM_r
        for j=1:1+DIM_t
            for j_=j+1:1+DIM_t
                e_ii_ = zeros(DIM_r+1, DIM_r+1);
                e_ii_(i,i_) = 1;
                e_jj_ = zeros(DIM_t+1, DIM_t+1);
                e_jj_(j,j_) = 1;
                P = kron((e_jj_ - e_jj_'), e_ii_);
                temp = zeros(DIM_z_HOMO, DIM_z_HOMO);
                temp(X_1_rt_zIndex, X_1_rt_zIndex) = P;
                Q_aux_1{end+1} = sparse(temp);
            end
        end
    end
end
% second, do X^1_ttheta
% create appropriate indexing from X^1_ttheta to z
X_1_ttheta_zIndex = [1,1+1:1+DIM_t];
for i=1:DIM_theta
    temp=[1+DIM_t+DIM_r+i,...
        1+DIM_t+DIM_r+DIM_theta+DIM_tr+(i-1)*DIM_t+1:...
        1+DIM_t+DIM_r+DIM_theta+DIM_tr+(i-1)*DIM_t+DIM_t];
    X_1_ttheta_zIndex=horzcat(X_1_ttheta_zIndex,temp);
end
for i=1:1+DIM_t
    for i_=i+1:1+DIM_t
        for j=1:1+DIM_theta
            for j_=j+1:1+DIM_theta
                e_ii_ = zeros(DIM_t+1, DIM_t+1);
                e_ii_(i,i_) = 1;
                e_jj_ = zeros(DIM_theta+1, DIM_theta+1);
                e_jj_(j,j_) = 1;
                P = kron((e_jj_ - e_jj_'), e_ii_);
                temp = zeros(DIM_z_HOMO, DIM_z_HOMO);
                temp(X_1_ttheta_zIndex,X_1_ttheta_zIndex)=P;
                Q_aux_1{end+1}=sparse(temp);
            end
        end
    end
end
% third, do X^1_rtheta
% create appropriate indexing from X^1_rtheta to z
X_1_rtheta_zIndex = [1,1+DIM_t+1:1+DIM_t+DIM_r];
for i=1:DIM_theta
    temp=[1+DIM_t+DIM_r+i,...
        1+DIM_t+DIM_r+DIM_theta+DIM_tr+DIM_thetat+(i-1)*DIM_r+1:...
        1+DIM_t+DIM_r+DIM_theta+DIM_tr+DIM_thetat+(i-1)*DIM_r+DIM_r];
    X_1_rtheta_zIndex=horzcat(X_1_rtheta_zIndex,temp);
end
for i=1:1+DIM_r
    for i_=i+1:1+DIM_r
        for j=1:1+DIM_theta
            for j_=j+1:1+DIM_theta
                e_ii_ = zeros(DIM_r+1, DIM_r+1);
                e_ii_(i,i_) = 1;
                e_jj_ = zeros(DIM_theta+1, DIM_theta+1);
                e_jj_(j,j_) = 1;
                P = kron((e_jj_ - e_jj_'), e_ii_);
                temp = zeros(DIM_z_HOMO, DIM_z_HOMO);
                temp(X_1_rtheta_zIndex,X_1_rtheta_zIndex)=P;
                Q_aux_1{end+1}=sparse(temp);
            end
        end
    end
end

%% construct auxiliary constraints degree 2
Q_aux_2={};
% first, do X^2_theta_t_r
% create appropriate indexing from X^2_theta_t_r to z
X_2_thetatr_zIndex=[1,1+DIM_t+DIM_r+1:1+DIM_t+DIM_r+DIM_theta];
for i=1:DIM_t
    for j=1:DIM_r
        temp=[1+DIM_t+DIM_r+DIM_theta+(i-1)*DIM_r+j,...
            1+DIM_t+DIM_r+DIM_theta+DIM_tr+DIM_thetat+DIM_thetar+(i-1)*DIM_r+j:DIM_tr:...
            1+DIM_t+DIM_r+DIM_theta+DIM_tr+DIM_thetat+DIM_thetar+(i-1)*DIM_r+j+(DIM_theta-1)*DIM_tr];
        X_2_thetatr_zIndex=horzcat(X_2_thetatr_zIndex,temp);
    end
end
% fix first row and first col
for i=1:1
    for i_=i+1:DIM_theta+1
        for j=1:1
            for j_=j+1:DIM_tr+1
                e_ii_ = zeros(DIM_theta+1, DIM_theta+1);
                e_ii_(i,i_) = 1;
                e_jj_ = zeros(DIM_tr+1, DIM_tr+1);
                e_jj_(j,j_) = 1;
                P = kron((e_jj_ - e_jj_'), e_ii_);
                temp = zeros(DIM_z_HOMO, DIM_z_HOMO);
                temp(X_2_thetatr_zIndex,X_2_thetatr_zIndex)=P;
                Q_aux_2{end+1}=sparse(temp);
            end
        end
    end
end
% Only fix first col
for i=2:1+DIM_theta
    for i_=i+1:1+DIM_theta
        for j=1:1
            for j_=j+1:DIM_tr
                e_ii_ = zeros(DIM_theta+1, DIM_theta+1);
                e_ii_(i,i_) = 1;
                e_jj_ = zeros(DIM_tr+1, DIM_tr+1);
                e_jj_(j,j_) = 1;
                P = kron((e_jj_ - e_jj_'), e_ii_);
                temp = zeros(DIM_z_HOMO, DIM_z_HOMO);
                temp(X_2_thetatr_zIndex,X_2_thetatr_zIndex)=P;
                Q_aux_2{end+1}=sparse(temp);
            end
        end
    end
end
% only fix first row
for i=1:1
    for i_=i+1:1+DIM_theta
        for j=2:1+DIM_tr
            for j_=j+1:1+DIM_tr
                e_ii_ = zeros(DIM_theta+1, DIM_theta+1);
                e_ii_(i,i_) = 1;
                e_jj_ = zeros(DIM_tr+1, DIM_tr+1);
                e_jj_(j,j_) = 1;
                P = kron((e_jj_ - e_jj_'), e_ii_);
                temp = zeros(DIM_z_HOMO, DIM_z_HOMO);
                temp(X_2_thetatr_zIndex,X_2_thetatr_zIndex)=P;
                Q_aux_2{end+1}=sparse(temp);
            end
        end
    end
end

% second, do X^2_t_theta_r
% create appropriate indexing from X^2_t_theta_r to z
X_2_tthetar_zIndex=1:1+DIM_t;
for i=1:DIM_theta
    for j=1:DIM_r
        temp=[1+DIM_t+DIM_r+DIM_theta+DIM_tr+DIM_thetat+(i-1)*DIM_r+j,...
            1+DIM_t+DIM_r+DIM_theta+DIM_tr+DIM_thetat+DIM_thetar+(i-1)*DIM_tr+j:DIM_r:...
            1+DIM_t+DIM_r+DIM_theta+DIM_tr+DIM_thetat+DIM_thetar+(i-1)*DIM_tr+j+2*DIM_r];
        X_2_tthetar_zIndex=horzcat(X_2_tthetar_zIndex,temp);
    end
end
% Fix first row and first col
for i=1:1
    for i_=i+1:1+DIM_t
        for j=1:1
            for j_=j+1:1+DIM_thetar
                e_ii_ = zeros(DIM_t+1, DIM_t+1);
                e_ii_(i,i_) = 1;
                e_jj_ = zeros(DIM_thetar+1, DIM_thetar+1);
                e_jj_(j,j_) = 1;
                P = kron((e_jj_ - e_jj_'), e_ii_);
                temp = zeros(DIM_z_HOMO, DIM_z_HOMO);
                temp(X_2_tthetar_zIndex,X_2_tthetar_zIndex)=P;
                Q_aux_2{end+1}=sparse(temp);
            end
        end
    end
end
% Only fix first col
for i=2:1+DIM_t
    for i_=i+1:1+DIM_t
        for j=1:1
            for j_=j+1:DIM_thetar
                e_ii_ = zeros(DIM_t+1, DIM_t+1);
                e_ii_(i,i_) = 1;
                e_jj_ = zeros(DIM_thetar+1, DIM_thetar+1);
                e_jj_(j,j_) = 1;
                P = kron((e_jj_ - e_jj_'), e_ii_);
                temp = zeros(DIM_z_HOMO, DIM_z_HOMO);
                temp(X_2_tthetar_zIndex,X_2_tthetar_zIndex)=P;
                Q_aux_2{end+1}=sparse(temp);
            end
        end
    end
end
% only fix first row
for i=1:1
    for i_=i+1:1+DIM_t
        for j=2:1+DIM_thetar
            for j_=j+1:1+DIM_thetar
                e_ii_ = zeros(DIM_t+1, DIM_t+1);
                e_ii_(i,i_) = 1;
                e_jj_ = zeros(DIM_thetar+1, DIM_thetar+1);
                e_jj_(j,j_) = 1;
                P = kron((e_jj_ - e_jj_'), e_ii_);
                temp = zeros(DIM_z_HOMO, DIM_z_HOMO);
                temp(X_2_tthetar_zIndex,X_2_tthetar_zIndex)=P;
                Q_aux_2{end+1}=sparse(temp);
            end
        end
    end
end
% third, do X^2_r_theta_t
% create appropriate indexing from X^2_r_theta_t to z
X_2_rthetat_zIndex=[1,1+DIM_t+1:1+DIM_t+DIM_r];
for i=1:DIM_theta
    for j=1:DIM_t
        temp=[1+DIM_t+DIM_r+DIM_theta+DIM_tr+(i-1)*DIM_t+j,...
            1+DIM_t+DIM_r+DIM_theta+DIM_tr+DIM_thetat+DIM_thetar+(i-1)*DIM_tr+(j-1)*DIM_r+1:...
            1+DIM_t+DIM_r+DIM_theta+DIM_tr+DIM_thetat+DIM_thetar+(i-1)*DIM_tr+(j-1)*DIM_r+DIM_r];
        X_2_rthetat_zIndex=horzcat(X_2_rthetat_zIndex,temp);
    end
end
% fix first row and first col
for i=1:1
    for i_=i+1:1+DIM_r
        for j=1:1
            for j_=j+1:1+DIM_thetat
                e_ii_ = zeros(DIM_r+1, DIM_r+1);
                e_ii_(i,i_) = 1;
                e_jj_ = zeros(DIM_thetat+1, DIM_thetat+1);
                e_jj_(j,j_) = 1;
                P = kron((e_jj_ - e_jj_'), e_ii_);
                temp = zeros(DIM_z_HOMO, DIM_z_HOMO);
                temp(X_2_rthetat_zIndex,X_2_rthetat_zIndex)=P;
                Q_aux_2{end+1}=sparse(temp);
            end
        end
    end
end
% Only fix first col
for i=2:1+DIM_r
    for i_=i+1:1+DIM_r
        for j=1:1
            for j_=j+1:DIM_thetat
                e_ii_ = zeros(DIM_r+1, DIM_r+1);
                e_ii_(i,i_) = 1;
                e_jj_ = zeros(DIM_thetat+1, DIM_thetat+1);
                e_jj_(j,j_) = 1;
                P = kron((e_jj_ - e_jj_'), e_ii_);
                temp = zeros(DIM_z_HOMO, DIM_z_HOMO);
                temp(X_2_rthetat_zIndex,X_2_rthetat_zIndex)=P;
                Q_aux_2{end+1}=sparse(temp);
            end
        end
    end
end
% only fix first row
for i=1:1
    for i_=i+1:1+DIM_r
        for j=2:1+DIM_thetat
            for j_=j+1:1+DIM_thetat
                e_ii_ = zeros(DIM_r+1, DIM_r+1);
                e_ii_(i,i_) = 1;
                e_jj_ = zeros(DIM_thetat+1, DIM_thetat+1);
                e_jj_(j,j_) = 1;
                P = kron((e_jj_ - e_jj_'), e_ii_);
                temp = zeros(DIM_z_HOMO, DIM_z_HOMO);
                temp(X_2_rthetat_zIndex,X_2_rthetat_zIndex)=P;
                Q_aux_2{end+1}=sparse(temp);
            end
        end
    end
end

            
%% construct cost matrix Q_cost
Q_cost_1 = {};
Q_cost_2 = {};
Q_cost_3 = {};
Q_cost_4 = {};
Q_cost_5 = {};
Q_cost = {};
% try some numerical tricks
c_real = 1e-3;
numeri_scale=1;
c=sqrt(numeri_scale*c_real^2);

% convert tr^T to z index
tr_cost_zIndex = (1+DIM_t+DIM_r+DIM_theta+1):DIM_r:...
    (1+DIM_t+DIM_r+DIM_theta+1+2*DIM_r);
for i=2:DIM_r
    temp = (1+DIM_t+DIM_r+DIM_theta+i):DIM_r:...
        (1+DIM_t+DIM_r+DIM_theta+i+2*DIM_r);
    tr_cost_zIndex = horzcat(tr_cost_zIndex, temp);
end
% convert theta_i * tr^T to z index for every theta_i
theta_tr_cost_zIndex={};
for i=1:DIM_theta
    temp=(1+DIM_t+DIM_r+DIM_theta+DIM_tr+DIM_thetat+DIM_thetar+(i-1)*DIM_tr+1):DIM_r:...
        (1+DIM_t+DIM_r+DIM_theta+DIM_tr+DIM_thetat+DIM_thetar+(i-1)*DIM_tr+1+2*DIM_r);
    for j=2:DIM_r
        temp_temp=(1+DIM_t+DIM_r+DIM_theta+DIM_tr+DIM_thetat+DIM_thetar+(i-1)*DIM_tr+j):DIM_r:...
        (1+DIM_t+DIM_r+DIM_theta+DIM_tr+DIM_thetat+DIM_thetar+(i-1)*DIM_tr+j+2*DIM_r);
        temp=horzcat(temp,temp_temp);
    end
    theta_tr_cost_zIndex{i}=temp;
end
% convert rt^T to z index for every theta_i
rt_cost_zIndex = (1+DIM_t+DIM_r+DIM_theta+1):...
    (1+DIM_t+DIM_r+DIM_theta+DIM_tr);
for i=1:NUM_POINTS
    f2 = bearing2(:,i);
    f1 = bearing1(:,i);
    G = kron(f2', hatmap(f1));
    A = numeri_scale * kron(G', G);
    temp_1 = zeros(DIM_z_HOMO, DIM_z_HOMO);
    temp_2 = zeros(DIM_z_HOMO, DIM_z_HOMO);
    temp_3 = zeros(DIM_z_HOMO, DIM_z_HOMO);
    temp_4 = zeros(DIM_z_HOMO, DIM_z_HOMO);
    temp_5 = zeros(DIM_z_HOMO, DIM_z_HOMO);
    temp_1(tr_cost_zIndex, rt_cost_zIndex) = A/2;
    theta_i_rt_cost_zIndex = (1+DIM_t+DIM_r+DIM_theta+DIM_tr+DIM_thetat+DIM_thetar+(i-1)*DIM_tr+1):...
    (1+DIM_t+DIM_r+DIM_theta+DIM_tr+DIM_thetat+DIM_thetar+(i-1)*DIM_tr+DIM_tr);

    theta_i_tr_cost_zIndex=theta_tr_cost_zIndex{i};

    temp_2(tr_cost_zIndex, theta_i_rt_cost_zIndex) = A/2;
    temp_3(1, 1+DIM_t+DIM_r+i) = -c^2/2;
    temp_4(1,1) = c^2/2;
    temp_5(theta_i_tr_cost_zIndex,theta_i_rt_cost_zIndex)=A/2;
    Q_cost_1{end+1} = sparse(temp_1);
    Q_cost_2{end+1} = sparse(temp_2);
    Q_cost_3{end+1} = sparse(temp_3);
    Q_cost_4{end+1} = sparse(temp_4);
    Q_cost_5{end+1} = sparse(temp_5);
    temp = temp_1 + temp_2 + temp_3 + temp_4;
    Q_cost{end+1} = sparse(temp);
end
Q_cost_total = zeros(DIM_z_HOMO, DIM_z_HOMO);
for i=1:length(Q_cost)
    Q_cost_total = Q_cost_total + Q_cost{i};
end
Q_cost_total=sparse(Q_cost_total);

%% print a few information
disp('Size of decision matrix:')
disp(DIM_z_HOMO)
disp('Number of lifted constraints:')
disp('-----Q_rot-------')
disp(length(Q_rot))
disp('-----Q_sphere-------')
disp(length(Q_sphere))
disp('-----Q_binary-------')
disp(length(Q_binary))
disp('-----Q_aux_1-------')
disp(length(Q_aux_1))
disp('-----Q_aux_2-------')
disp(length(Q_aux_2))

%% DO sanity check using ground truth
errorMatrix = {};
tolerance = 1e-12;
r_gt = reshape(R_gt, [9,1]);
t_gt = t_gt/norm(t_gt);
theta_gt = ones(NUM_POINTS,1); % accept all the correspondences
theta_gt(1)=-1;
z_gt = [1;t_gt;r_gt;theta_gt;...
    kron(t_gt,r_gt);kron(theta_gt,t_gt);kron(theta_gt,r_gt);...
    kron(theta_gt, kron(t_gt,r_gt))];
ZHomo_gt = z_gt * z_gt';
cost_gt = trace(Q_cost_total * ZHomo_gt);
disp('------Ground truth Cost-----')
disp(cost_gt)
if sanityCheck
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
    disp('-----------------END Sanity Check--------------')
end

%% begin CVX solve
thetaIndex=1+DIM_t+DIM_r+1:1+DIM_t+DIM_r+NUM_POINTS;
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
sum(ZHomo(1,thetaIndex)) >= (12-NUM_POINTS)

disp('---ADD Cost CONSTRAINTS---')
for i=1:length(Q_cost)
    disp(i)
    
    temp=Q_cost{i};
    temp1=Q_cost_1{i};
    temp2=Q_cost_2{i};
    temp3=Q_cost_3{i};
    temp4=Q_cost_4{i};
    temp5=Q_cost_5{i};
    
    [row,col,val]=find(temp);
    expression x
    x=0;
    for nonzero=1:length(row)
        x = x + val(nonzero) * ZHomo(col(nonzero), row(nonzero));
    end
    x>=0;
    
    [row,col,val]=find(temp1);
    expression x
    x=0;
    for nonzero=1:length(row)
        x = x + val(nonzero) * ZHomo(col(nonzero), row(nonzero));
    end
    x>=0;
    
    [row,col,val]=find(temp1+temp2);
    expression x
    x=0;
    for nonzero=1:length(row)
        x = x + val(nonzero) * ZHomo(col(nonzero), row(nonzero));
    end
    x>=0;
    
    [row,col,val]=find(temp3+temp4);
    expression x
    x=0;
    for nonzero=1:length(row)
        x = x + val(nonzero) * ZHomo(col(nonzero), row(nonzero));
    end
    x>=0;
    
    [row,col,val]=find(temp5-temp1);
    expression x
    x=0;
    for nonzero=1:length(row)
        x=x+val(nonzero)*ZHomo(col(nonzero), row(nonzero));
    end
    x==0;
    
%     trace(Q_cost{i} * ZHomo) >= 0
%     trace(Q_cost_1{i} * ZHomo) >= 0
%     trace((Q_cost_1{i}+Q_cost_2{i})*ZHomo) >= 0
%     trace((Q_cost_3{i}+Q_cost_4{i})*ZHomo) >= 0
end

disp('---ADD ROT CONSTRAINTS---')
for i=1:length(Q_rot)
    if rem(i,1000)==1
        disp(i)
    end
    temp=Q_rot{i};
    [row,col,val]=find(temp);
    expression x
    x=0;
    for nonzero=1:length(row)
        x = x + val(nonzero) * ZHomo(col(nonzero), row(nonzero));
    end
    x==0;
    
%     trace(Q_rot{i} * sparse(ZHomo)) == 0;
end
disp('---ADD TRANS CONSTRAINTS---')
for i=1:length(Q_sphere)
    if rem(i,1000)==1
        disp(i)
    end
    temp=Q_sphere{i};
    [row,col,val]=find(temp);
    expression x
    x=0;
    for nonzero=1:length(row)
        x=x+val(nonzero) * ZHomo(col(nonzero), row(nonzero));
    end
    x==0;
    
%     trace(Q_sphere{i} * ZHomo) == 0;
end
disp('---ADD BINARY CONSTRAINTS---')
for i=1:length(Q_binary)
    if rem(i,1000)==1
        disp(i)
    end
    temp=Q_binary{i};
    [row,col,val]=find(temp);
    expression x
    x=0;
    for nonzero=1:length(row)
        x=x+val(nonzero) * ZHomo(col(nonzero), row(nonzero));
    end
    x==0;
    
%     trace(Q_binary{i} * ZHomo) == 0;
end
disp('---ADD AUX 1 CONSTRAINTS---')
for i=1:length(Q_aux_1)
    if rem(i,1000)==1
        disp(i)
    end
    temp=Q_aux_1{i};
    [row,col,val]=find(temp);
    expression x
    x=0;
    for nonzero=1:length(row)
        x=x+val(nonzero) * ZHomo(col(nonzero), row(nonzero));
    end
    x==0;
    
%     trace(Q_aux_1{i} * ZHomo) == 0;
end
disp('---ADD AUX 2 CONSTRAINTS---')
for i=1:length(Q_aux_2)
    if rem(i,1000)==1
        disp(i)
    end
    temp=Q_aux_2{i};
    [row,col,val]=find(temp);
    expression x
    x=0;
    for nonzero=1:length(row)
        x=x+val(nonzero) * ZHomo(col(nonzero), row(nonzero));
    end
    x==0;
%     trace(Q_aux_2{i} * ZHomo) == 0;
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
    
    
            
        
        

            
