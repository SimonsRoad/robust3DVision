% This script benchmarks the following two view geometry methods
% 5-point Nister + nonlinear optimization refine, from OPENGV
% 8-point + nonlinear optimization refine, from OPENGV
% N-point QCQP relaxation, from Briales
% Test condition: 30 runs, noise level from 0 to high, no outlier
% benchmark comparison: pose error, rotation error and translation error

%% start clean
clc
clear
close all
cvx_solver mosek
cvx_save_prefs

addpath(genpath('../lib'))

%% benchmark conditions
NUM_ITERATIONS = 30; % run different methods for multiple iterations
NOISE_LEVEL = 10; % in pixel
FIELD_OF_VIEW = 100; % in degress
NUM_POINTS = 100;
MAX_PARALLAX = 2; % in meters
OUTLIER_RATIO = 0.0;
DEBUG_SYNTHETIC = false; % true to turn plot on
QCQP_QUIET = true;

five_pt_idx = 1:5; % use only the first 5 observations
eight_pt_idx = 1:8; % use only the first 8 observations
nonlin_opt_idx = 1:NUM_POINTS; % use all points to refine 5/8-pt
qcqp_idx = 1:NUM_POINTS; % use all points for QCQP

for i=1:length(NOISE_LEVEL)
    noise_level = NOISE_LEVEL(i);
    for itr = 1:NUM_ITERATIONS
        %% first generate bearing vectors and ground truth poses
        [bearing1, bearing2, t_gt, R_gt] = create2ViewCentralExperiment(...
                noise_level, NUM_POINTS, FIELD_OF_VIEW, MAX_PARALLAX,...
                OUTLIER_RATIO, DEBUG_SYNTHETIC);
        t_gt_norm = norm(t_gt);
        t_gt = t_gt / t_gt_norm;
        disp('--------Ground Truth---------')
        disp([R_gt,t_gt])
        
        %{
        %% use 5-point + nonlinear refinement
        % first use 5-point to get initial estimation
        X_five_pt=opengv('fivept_nister',five_pt_idx,bearing1,bearing2);
        % transform the essential matrices to rotation matrices
        R_five_pt=transformEssentials(X_five_pt);
        % take the R that is closest to R_gt
        temp=zeros(1,size(R_five_pt,3));
        for k=1:size(R_five_pt,3)
            temp(k)=evaluateRotationError(R_gt,R_five_pt(:,:,k));
        end
        [~,minIdx]=min(temp);
        R_five_pt_best=R_five_pt(:,:,minIdx);
        E_five_pt_best=X_five_pt(:,:,ceil(minIdx/2));
        % use best R to get t
        t_five_pt_best=vmap(E_five_pt_best*R_five_pt_best');
        t_five_pt_best=t_five_pt_best/norm(t_five_pt_best);
        if norm(-t_five_pt_best-t_gt) < norm(t_five_pt_best-t_gt)
            t_five_pt_best=-t_five_pt_best;
        end
        % feed result from 5-point to nonlin opt
        X_5pt_nonlin_opt=opengv('rel_nonlin_central',nonlin_opt_idx,...
            bearing1,bearing2,[R_five_pt_best,t_five_pt_best]);
        X_5pt_nonlin_opt(:,4)=...
            X_5pt_nonlin_opt(:,4)/norm(X_5pt_nonlin_opt(:,4));
        disp('-------five pt--------')
        disp(X_5pt_nonlin_opt)
        % calculate rotation and translation error, and store
        rot_error(i,itr,1)=evaluateRotationError(R_gt,...
            X_5pt_nonlin_opt(:,1:3));
        trans_error(i,itr,1)=norm(t_gt-X_5pt_nonlin_opt(:,4));     
        
        %% use 8-point + nonlinear refinement
        % first use 8-point to get initial estimation
        X_eight_pt=opengv('eightpt',eight_pt_idx,bearing1,bearing2);
        % transform essential to rotation matrices
        R_eight_pt=transformEssentials(X_eight_pt);
        % take the best R compared to gt
        temp=zeros(1,size(R_eight_pt,3));
        for k=1:size(R_eight_pt,3)
            temp(k)=evaluateRotationError(R_gt,R_eight_pt(:,:,k));
        end
        [~,minIdx]=min(temp);
        R_eight_pt_best=R_eight_pt(:,:,minIdx);
        E_eight_pt_best=X_eight_pt(:,:,ceil(minIdx/2));
        % use best R to get t
        t_eight_pt_best=vmap(E_eight_pt_best*R_eight_pt_best');
        t_eight_pt_best=t_eight_pt_best/norm(t_eight_pt_best);
        if norm(-t_eight_pt_best-t_gt) < norm(t_eight_pt_best-t_gt)
            t_eight_pt_best=-t_eight_pt_best;
        end
        % feed result from 5-point to nonlin opt
        X_8pt_nonlin_opt=opengv('rel_nonlin_central',nonlin_opt_idx,...
            bearing1,bearing2,[R_eight_pt_best,t_eight_pt_best]);
        X_8pt_nonlin_opt(:,4)=...
            X_8pt_nonlin_opt(:,4)/norm(X_8pt_nonlin_opt(:,4));
        disp('-------eight pt--------')
        disp(X_8pt_nonlin_opt)
        % calculate rotation and translation error, and store
        rot_error(i,itr,2)=evaluateRotationError(R_gt,...
            X_8pt_nonlin_opt(:,1:3));
        trans_error(i,itr,2)=norm(t_gt-X_8pt_nonlin_opt(:,4));
        %}
        %% use Briales' QCQP relaxation solver
        [RStarQCQP,tStarQCQP,ZStarQCQP,recoverStatusQCQP,...
            cvxStatusQCQP,pStarRelaxQCQP,tightnessQCQP] = ...
            qcqprelaxation(bearing1, bearing2, QCQP_QUIET);
        disp(cvxStatusQCQP)
        disp('QCQP relax tightness:')
        disp(tightnessQCQP)
        if recoverStatusQCQP
            rotErrorQCQP = zeros(1,size(RStarQCQP, 2));
            transErrorQCQP = zeros(1,size(RStarQCQP, 2));
            combErrorQCQP = zeros(1,size(RStarQCQP, 2));
            for k=1:size(RStarQCQP, 2)
                rotErrorQCQP(k) = norm(RStarQCQP{k} - R_gt, 'fro');
                transErrorQCQP(k) = norm(tStarQCQP{k} - t_gt);
                combErrorQCQP(k) = sqrt(rotErrorQCQP(k)^2+transErrorQCQP(k)^2);
            end
            [~, bestIdx] = min(combErrorQCQP);
            bestRQCQP = RStarQCQP{bestIdx};
            besttQCQP = tStarQCQP{bestIdx};
            disp('----------QCQP Relaxation----------')
            disp([bestRQCQP,besttQCQP]);
            rot_error(i,itr,3)=evaluateRotationError(R_gt,...
                bestRQCQP);
            trans_error(i,itr,3)=transErrorQCQP(bestIdx);
        end     
    end  
end

%{
%% plot translational and rotational error box plot
rot_error_plot=zeros(NUM_ITERATIONS,3*length(NOISE_LEVEL));
trans_error_plot=zeros(NUM_ITERATIONS,3*length(NOISE_LEVEL));
for i=1:length(NOISE_LEVEL)
    rot_error_plot(:,(i-1)*3+1:(i-1)*3+3)=rot_error(i,:,:);
    trans_error_plot(:,(i-1)*3+1:(i-1)*3+3)=trans_error(i,:,:);
end
g1=kron(ones(1,length(NOISE_LEVEL)), 1:3);
g2=kron(1:length(NOISE_LEVEL), ones(1,3));
% rotational error
rot_plot=figure(1);
boxplot(rot_error_plot,{g2,g1},...
    'Colors','rgb','FactorGap',[2,0],'Widths',[0.8,0.8,0.8]);
set(findobj(gca(rot_plot),'type','line'),'linew',1.5)
grid(gca(rot_plot),'on')
set(gca(rot_plot),'xTickLabel','')
xlabel(gca(rot_plot), 'Noise [pix]')
ylabel(gca(rot_plot), 'Rotational Error')
h1=gca(rot_plot);
ylim(gca(rot_plot), [0, h1.YLim(2)])
legend(findobj(gca(rot_plot),'Tag','Box'),...
    'QCQP','8-point+Nonlinear Opt.','5-point+Nonlinear Opt.')

% translational error plot
trans_plot=figure(2);
boxplot(trans_error_plot,{g2,g1},...
    'Colors','rgb','FactorGap',[2,0],'Widths',[0.8,0.8,0.8]);
set(findobj(gca(trans_plot),'type','line'),'linew',1.5)
grid(gca(trans_plot), 'on')
set(gca(trans_plot),'xTickLabel','')
xlabel(gca(trans_plot), 'Noise [pix]')
ylabel(gca(trans_plot), 'Translational Error')
h2=gca(trans_plot);
ylim(gca(trans_plot), [0, h2.YLim(2)])
legend(findobj(gca(trans_plot),'Tag','Box'),...
    'QCQP','8-point+Nonlinear Opt.','5-point+Nonlinear Opt.')
%}        
        
        
        
