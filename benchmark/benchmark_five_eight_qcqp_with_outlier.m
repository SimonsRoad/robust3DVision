% This script benchmarks the following two view geometry methods
% 5-point Nister + RANSAC, from OPENGV
% 8-point + RANSAC, from OPENGV
% N-point QCQP relaxation + RANSAC
% Test condition: 30 runs, noise level fixed (low or high), outlier level
% increasing
% benchmark comparison: pose error, rotation error and translation error

%% start clean
clc
clear
close all
cvx_solver mosek
cvx_save_prefs

%% benchmark conditions
NUM_ITERATIONS = 20; % run different methods for multiple iterations
NOISE_LEVEL = 20; % in pixel
FIELD_OF_VIEW = 100; % in degress
NUM_POINTS = 100;
MAX_PARALLAX = 2; % in meters
OUTLIER_RATIO = 0:0.1:0.3;
DEBUG_SYNTHETIC = false; % true to turn plot on
QCQP_QUIET = true;
NUM_QCQP=10;
% ransac threshold
INLIER_THRESHOLD=1e-4;
MAX_RANSAC=200;

%% run benchmark
rot_error=zeros(length(OUTLIER_RATIO),NUM_ITERATIONS,3);
trans_error=zeros(length(OUTLIER_RATIO),NUM_ITERATIONS,3);
ransac_steps=zeros(length(OUTLIER_RATIO),NUM_ITERATIONS,3);
for i=1:length(OUTLIER_RATIO)
    outlier_ratio = OUTLIER_RATIO(i);
    for itr = 1:NUM_ITERATIONS
        %% first generate bearing vectors and ground truth poses
        [bearing1, bearing2, t_gt, R_gt] = create2ViewCentralExperiment(...
                NOISE_LEVEL, NUM_POINTS, FIELD_OF_VIEW, MAX_PARALLAX,...
                outlier_ratio, DEBUG_SYNTHETIC);
        t_gt_norm = norm(t_gt);
        t_gt = t_gt / t_gt_norm;
        disp('--------Ground Truth---------')
        disp([R_gt,t_gt])
        
         %% use 5-point + RANSAC
        [R_5pt_ransac,t_5pt_ransac,nItrs_5pt]=fivept_ransac(...
            bearing1,bearing2,R_gt,t_gt,INLIER_THRESHOLD,MAX_RANSAC);
        disp('------FIVE point RANSAC------')
        disp([R_5pt_ransac,t_5pt_ransac])
        % calculate rotation and translation error, and store
        rot_error(i,itr,1)=evaluateRotationError(R_gt,R_5pt_ransac);
        trans_error(i,itr,1)=norm(t_gt-t_5pt_ransac);
        ransac_steps(i,itr,1)=nItrs_5pt;
        %% use 8-point + RANSAC
        [R_8pt_ransac,t_8pt_ransac,nItrs_8pt]=eightpt_ransac(...
            bearing1,bearing2,R_gt,t_gt,INLIER_THRESHOLD,MAX_RANSAC);
        disp('------Eight point RANSAC------')
        disp([R_8pt_ransac,t_8pt_ransac])
        % calculate rotation and translation error, and store
        rot_error(i,itr,2)=evaluateRotationError(R_gt,R_8pt_ransac);
        trans_error(i,itr,2)=norm(t_gt-t_8pt_ransac);
        ransac_steps(i,itr,2)=nItrs_8pt;
        %% use Briales' QCQP relaxation + RANSAC
        [R_qcqp_ransac,t_qcqp_ransac,nItrs_qcqp]=qcqprelaxation_ransac(...
            bearing1,bearing2,NUM_QCQP,...
            R_gt,t_gt,INLIER_THRESHOLD,MAX_RANSAC);
        disp('--------QCQP RANSAC-------')
        disp([R_qcqp_ransac,t_qcqp_ransac])
        % calculate rotation and translation error, and store
        rot_error(i,itr,3)=evaluateRotationError(R_gt,R_qcqp_ransac);
        trans_error(i,itr,3)=norm(t_gt-t_qcqp_ransac);
        ransac_steps(i,itr,3)=nItrs_qcqp;  
    end  
end

%% plot translational and rotational error box plot
rot_error_plot=zeros(NUM_ITERATIONS,3*length(OUTLIER_RATIO));
trans_error_plot=zeros(NUM_ITERATIONS,3*length(OUTLIER_RATIO));
ransac_steps_plot=zeros(NUM_ITERATIONS,3*length(OUTLIER_RATIO));
for i=1:length(OUTLIER_RATIO)
    rot_error_plot(:,(i-1)*3+1:(i-1)*3+3)=rot_error(i,:,:);
    trans_error_plot(:,(i-1)*3+1:(i-1)*3+3)=trans_error(i,:,:);
    ransac_steps_plot(:,(i-1)*3+1:(i-1)*3+3)=ransac_steps(i,:,:);
end
g1=kron(ones(1,length(OUTLIER_RATIO)), 1:3);
g2=kron(1:length(OUTLIER_RATIO), ones(1,3));
% rotational error
rot_plot=figure(1);
boxplot(rot_error_plot,{g2,g1},...
    'Colors','rgb','FactorGap',[2,0],'Widths',[0.8,0.8,0.8]);
set(findobj(gca(rot_plot),'type','line'),'linew',1.5)
grid(gca(rot_plot),'on')
set(gca(rot_plot),'xTickLabel','')
xlabel(gca(rot_plot), 'Outlier Ratio')
ylabel(gca(rot_plot), 'Rotational Error')
h1=gca(rot_plot);
ylim(gca(rot_plot), [0, h1.YLim(2)])
legend(findobj(gca(rot_plot),'Tag','Box'),...
    'QCQP+RANSAC','8-point+RANSAC','5-point+RANSAC')

% translational error
trans_plot=figure(2);
boxplot(trans_error_plot,{g2,g1},...
    'Colors','rgb','FactorGap',[2,0],'Widths',[0.8,0.8,0.8]);
set(findobj(gca(trans_plot),'type','line'),'linew',1.5)
grid(gca(trans_plot),'on')
set(gca(trans_plot),'xTickLabel','')
xlabel(gca(trans_plot), 'Outlier Ratio')
ylabel(gca(trans_plot), 'Translational Error')
h1=gca(trans_plot);
ylim(gca(trans_plot), [0, h1.YLim(2)])
legend(findobj(gca(trans_plot),'Tag','Box'),...
    'QCQP+RANSAC','8-point+RANSAC','5-point+RANSAC')

% number of steps taken for RANSAC
steps_plot=figure(3);
boxplot(ransac_steps_plot,{g2,g1},...
    'Colors','rgb','FactorGap',[2,0],'Widths',[0.8,0.8,0.8]);
set(findobj(gca(steps_plot),'type','line'),'linew',1.5)
grid(gca(steps_plot),'on')
set(gca(steps_plot),'xTickLabel','')
xlabel(gca(steps_plot), 'Outlier Ratio')
ylabel(gca(steps_plot), 'Number of RANSAC Steps')
h1=gca(steps_plot);
ylim(gca(steps_plot), [0, h1.YLim(2)])
legend(findobj(gca(steps_plot),'Tag','Box'),...
    'QCQP+RANSAC','8-point+RANSAC','5-point+RANSAC')

        
        
        
        
