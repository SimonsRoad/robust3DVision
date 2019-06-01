% Default truncated LS cost function with binary variables implementation.
% Author: Heng Yang           
% Last update: 12/03/2018
%% start clean
clc
clear all
close all
cvx_solver mosek
cvx_save_prefs

NUM_ITERATIONS = 1; % run SDP multiple times
NOISE_LEVEL = 0; % in pixel
FIELD_OF_VIEW = 100; % in degress
NUM_POINTS = 6;
MAX_PARALLAX = 2; % in meters
OUTLIER_RATIO = 0.0;
DEBUG_SYNTHETIC = false; % true to turn plot on
SDP_QUIET = true; % false to turn cvx output on
QCQP_QUIET = true;
TLS_QUIET = false;
if NOISE_LEVEL==0 && OUTLIER_RATIO==0
    SANITY_CHECK=true;
else
    SANITY_CHECK=false;
end
if NUM_POINTS < 6
    error('QCQP needs at least 6 points to work without outliers.')
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
    
     %% solver using original QCQP
    disp('--------FIRST SOLVE USING QCQP---------')
    [RStar_QCQP,tStar_QCQP,ZStar_QCQP,recoverStatus_QCQP,cvxStatus_QCQP,pStarRelax_QCQP,tightness_QCQP]=...
        qcqprelaxation(bearing1,bearing2,QCQP_QUIET);
    if isnan(pStarRelax_QCQP)
        error('This case is unsolvable, please resample!')
    end
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
    
    %% solve using QCQP TLS relax
    [cvxStatus_TLS,pStarRelax_TLS,ZStar_TLS] = ...
    qcqprelaxation_truncated_ls(bearing1,bearing2,R_gt,t_gt,SANITY_CHECK,TLS_QUIET);
 
end




