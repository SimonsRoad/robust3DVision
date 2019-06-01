% This script solves two view geometry using convex relaxation.
% Redundant constraints test, corresponding to Q_aux_2 degree 2
% having redundant (maybe) constraints
% Author: Heng Yang           
% Last update: 12/3/2018
%% start clean
clc
clear all
close all
cvx_solver mosek
cvx_save_prefs

NUM_ITERATIONS = 1; % run SDP multiple times
NOISE_LEVEL = 0; % in pixel
FIELD_OF_VIEW = 100; % in degress
NUM_POINTS = 8;
MAX_PARALLAX = 2; % in meters
OUTLIER_RATIO = 0.2;
DEBUG_SYNTHETIC = false; % true to turn plot on
SDP_QUIET = true; % false to turn cvx output on
QCQP_QUIET = true;
TLS_QUIET = false;
if NOISE_LEVEL==0 && OUTLIER_RATIO==0
    SANITY_CHECK=true;
else
    SANITY_CHECK=false;
end
if NUM_POINTS < 7
    error('QCQP needs at least 7 points to work with outliers.')
end

finalErrorQCQP = zeros(NUM_ITERATIONS, 3);
for itr=1:NUM_ITERATIONS
    %% generates relative rotation, tranlsation and bearing vectors
    [bearing1, bearing2, t_gt, R_gt] = create2ViewCentralExperiment(...
        NOISE_LEVEL, NUM_POINTS, FIELD_OF_VIEW, MAX_PARALLAX,...
        OUTLIER_RATIO, DEBUG_SYNTHETIC);   
    disp('true R')
    disp(R_gt)
    disp('true t')
    disp(t_gt/norm(t_gt))
    t_gt=t_gt/norm(t_gt);
    r_gt=reshape(R_gt,[9,1]);
    f1=bearing1(:,1);
    f2=bearing2(:,1);
    A=kron(f2',hatmap(f1));
    residual=(t_gt'*A*r_gt)^2;
    disp('--------OUTLIER RESIDUAL-------------')
    disp(residual);

    %% solver using original QCQP
    disp('--------FIRST SOLVE USING QCQP---------')
    [RStar_QCQP,tStar_QCQP,ZStar_QCQP,recoverStatus_QCQP,cvxStatus_QCQP,pStarRelax_QCQP,tightness_QCQP]=...
        qcqprelaxation(bearing1,bearing2,QCQP_QUIET);
    if isnan(pStarRelax_QCQP)
        error('Cost returned by QCQP is NAN.')
    end
    disp('------Cost obtained by QCQP-----')
    disp(pStarRelax_QCQP)
    if recoverStatus_QCQP
        rotErrorQCQP = zeros(1,size(RStar_QCQP, 2));
        transErrorQCQP = zeros(1,size(RStar_QCQP, 2));
        combErrorQCQP = zeros(1,size(RStar_QCQP, 2));
%         disp('------------recovered R and t----------------')
%         for i=1:size(RStar_QCQP, 2)
%             disp(RStar_QCQP{i})
%             disp(tStar_QCQP{i})
%         end
        for i=1:size(RStar_QCQP, 2)
            rotErrorQCQP(i) = norm(RStar_QCQP{i} - R_gt, 'fro');
            transErrorQCQP(i) = norm(tStar_QCQP{i} - t_gt/norm(t_gt));
            combErrorQCQP(i) = sqrt(rotErrorQCQP(i)^2+transErrorQCQP(i)^2);
        end
        [~, bestIdx] = min(combErrorQCQP);
        bestRQCQP = RStar_QCQP{bestIdx};
        besttQCQP = tStar_QCQP{bestIdx};
        bestErrorQCQP = [rotErrorQCQP(bestIdx), ...
            transErrorQCQP(bestIdx), combErrorQCQP(bestIdx)];
        disp('------------Best R estimation:------------')
        disp(bestRQCQP)
        disp('------------Best t estimation:------------')
        disp(besttQCQP)
        disp('------------Best error:-------------------')
        disp(bestErrorQCQP)
    end
    %% solve without outliers using QCQP
    disp('-----------Check sample quality using QCQP-------------')
    [RStar_QCQP_noOtl,tStar_QCQP_noOtl,ZStar_QCQP_noOtl,...
        recoverStatus_QCQP_noOtl,cvxStatus_QCQP_noOtl,...
        pStarRelax_QCQP_noOtl,tightness_QCQP_noOtl]=...
        qcqprelaxation(bearing1(:,2:end),bearing2(:,2:end),QCQP_QUIET);
    if isnan(pStarRelax_QCQP_noOtl)
        error('This case is unsolvable with QCQP, resample please.')
    end
    
    %% solve using truncated LS relax
    disp('-------SOLVE USING TRUNCATED LS--------')
    [cvxStatus_TLS,pStarRelax_TLS,ZStar_TLS] = ...
    qcqprelaxation_truncated_ls_outlier(bearing1,bearing2,R_gt,t_gt,SANITY_CHECK,TLS_QUIET);
end



