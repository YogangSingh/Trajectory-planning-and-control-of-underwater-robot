close all;
clear;
addpath('utils');
Path = @Trajectory;
waypoints = [0    0   0;
             1    2   3;
             3    3   1;
             4    4   5]'; 
Path([],[],waypoints);
Stabilizer = @PIDController;
[t, state, QP] = Simulator(Path, Stabilizer);


%% Functions
%Trajectory function
function [ desired_state ] = Trajectory(t, state, waypoints)
persistent waypoints0 traj_time allCoeffs multFact time d0
if nargin > 2
    waypoints0 = waypoints';
    n = size(waypoints0, 1) - 1;
    multFact = 8;
    avgSpeed = 2; %m/sec
    d = waypoints(:,2:end) - waypoints(:,1:end-1);
    d0 = 1/(avgSpeed) * sqrt(d(1,:).^2 + d(2,:).^2 + d(3,:).^2);
    traj_time = [0, cumsum(d0)];
    A = zeros(multFact * n); 
    b = zeros(multFact * n, 3);
    numPts = size(waypoints0,1);
    b(1,:) = waypoints0(1,:);
    rowP = 1;
    for iPosPt = 2:numPts-1
        row  = rowP + 1;
        rowP = row + 1;
        b(row,:)  = waypoints0(iPosPt,:);
        b(rowP,:) = waypoints0(iPosPt,:);
    end
    b(rowP+1,:) = waypoints0(end,:);    
    % Constraint values for A matrix
    posCon = [1 0 0 0 0 0 0 0;
        1 1 1 1 1 1 1 1];
    drvtvCon = [0 1 2 3 4 5 6 7;
        0 0 2 6 12 20 30 42;
        0 0 0 6 24 60 120 210];    
    continuousCon = [0 1 2 3 4 5 6 7 0 -1 0 0 0 0 0 0;
        0 0 2 6 12 20 30 42 0 0 -2 0 0 0 0 0;
        0 0 0 6 24 60 120 210 0 0 0 -6 0 0 0 0;
        0 0 0 0 24 120 360 840 0 0 0 0 -24 0 0 0;
        0 0 0 0 0 120 720 2520 0 0 0 0 0 -120 0 0;
        0 0 0 0 0 0 720 5040 0 0 0 0 0 0 -720 0];
    % Waypoint constraints in A
    rowEnd = 0;
    for iPos = 1:n
        colStart = (iPos-1)*multFact + 1;
        colEnd   = colStart + multFact - 1;
        rowStart = rowEnd + 1;
        rowEnd   = rowStart + 1;
        A(rowStart:rowEnd,colStart:colEnd) = posCon;
    end
    nxtRow = rowEnd;
    for iDrv = 1:3
        nxtRow = nxtRow + 1;
        A(nxtRow, iDrv+1) = 1;
    end
    
    startRow = (nxtRow + 1);
    endRow   = startRow + size(drvtvCon,1) - 1;
    startCol = (n-1) * multFact + 1;
    endCol   = size(A,2);
    A(startRow:endRow,startCol:endCol) = drvtvCon;
    
    colEnd = 0;
    for iCon = 1:n-1
        startRow = endRow + 1;
        endRow   = startRow + size(continuousCon,1) - 1;
        colStart = colEnd + 1;
        colEnd   = colStart + size(continuousCon,2) - 1;
        A(startRow:endRow,colStart:colEnd) = continuousCon;
        colEnd = multFact * iCon;
    end    
    allCoeffs = A\b;    
else
    if(t > traj_time(end))
        desired_state.pos = waypoints0(end,:)';
        desired_state.vel = zeros(3,1);
        desired_state.acc = zeros(3,1);
    else
        t_index = find(traj_time <= t,1,'last');
        
        if(t == 0)
            desired_state.pos = waypoints0(1,:)';
            desired_state.vel = zeros(3,1);
            desired_state.acc = zeros(3,1);
        else
            if(t_index > 1)
                t = t - traj_time(t_index);
            end
            S = t/d0(t_index);
            coeffs = allCoeffs(((t_index-1)*multFact + 1): (t_index*8), :);
            desired_state.pos = ([1, S, S^2, S^3, S^4, S^5, S^6, S^7] * coeffs)';
            desired_state.vel = ([0, 1, 2*S, 3*S^2, 4*S^3, 5*S^4, 6*S^5, 7*S^6] * coeffs)';
            desired_state.acc = ([0, 0, 2, 6*S, 12*S^2, 20*S^3, 30*S^4, 42*S^5] * coeffs)';
        end
    end
    desired_state.yaw = 0;
    desired_state.yawdot = 0;
end
end
%Controller Function
function [F, M] = PIDController(t, state, des_state, params)
pDes = 0; %Desired roll velocity
qDes = 0; %Desired pitch velocity
Kp_X = 100;
Kd_X = 10; 
Kp_Y = Kp_X
Kd_Y = Kd_X;
Kp_Z = 100; 
Kd_Z = 10; 
Kp_phi = 100;
Kd_phi = 5;
Kp_theta = Kp_phi;
Kd_theta = Kd_phi;
Kp_psi = 100;
Kd_psi = 5;
%% Acceleration Commands
r1CommAcc = des_state.acc(1) + Kd_X*(des_state.vel(1) - state.vel(1)) + Kp_X*(des_state.pos(1) - state.pos(1)); 
r2CommAcc = des_state.acc(2) + Kd_Y*(des_state.vel(2) - state.vel(2)) + Kp_Y*(des_state.pos(2) - state.pos(2));
r3CommAcc = des_state.acc(3) + Kd_Z*(des_state.vel(3) - state.vel(3)) + Kp_Z*(des_state.pos(3) - state.pos(3));
phiDes   = 1/params.gravity * (r1CommAcc * sin(des_state.yaw) - r2CommAcc * cos(des_state.yaw));
thetaDes = 1/params.gravity * (r1CommAcc * cos(des_state.yaw) + r2CommAcc * sin(des_state.yaw));
%% Thrust 
F = params.mass * (r3CommAcc + params.gravity);
%% Moment
u2Phi   = Kp_phi * (phiDes - state.rot(1)) + Kd_phi * (pDes - state.omega(1));
u2Theta = Kp_theta * (thetaDes - state.rot(2)) + Kd_theta * (qDes - state.omega(2));
u2Psi   = Kp_psi * (des_state.yaw - state.rot(3)) + Kd_psi * (des_state.yawdot - state.omega(3));
M = [u2Phi;
    u2Theta;
    u2Psi];
end

%Simulator function
function [t_out, s_out, QP] = Simulator(Path, Stabilizer)
addpath('utils');
real_time = true;
max_time = 50;
params = sys_params;
disp('Initializing figures...');
h_fig = figure;
h_3d = gca;
axis equal
grid on
view(48.8, 25.8);
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]')
AUVcolors = lines(256);
set(gcf,'Renderer','OpenGL')
disp('Setting initial conditions...');
tstep    = 0.01;
cstep    = 0.02;
max_iter = max_time/cstep;
nstep    = cstep/tstep;
time     = 0;
err = []; 
% Start and stop position
des_start = Path(0, []);
des_stop  = Path(inf, []);
stop_pos  = des_stop.pos;
x0    = init_state(des_start.pos, 0);
xtraj = zeros(max_iter*nstep, length(x0));
ttraj = zeros(max_iter*nstep, 1);
x       = x0;
pos_tol = 0.01;
vel_tol = 0.01;
disp('Simulation Running....');

for iter = 1:max_iter

    timeint = time:tstep:time+cstep;

    tic;

    % Initialize AUV plot
    if iter == 1
        QP = AUVPlot(1, x0, 0.1, 0.04, AUVcolors(1,:), max_iter, h_3d);
        current_state = stateToQd(x);
        desired_state = Path(time, current_state);
        QP.UpdateAUVPlot(x, [desired_state.pos; desired_state.vel], time);
        h_title = title(sprintf('iteration: %d, time: %4.2f', iter, time));
    end

    % Run simulation
    [tsave, xsave] = ode45(@(t,s) AUVEOM(t, s, Stabilizer, Path, params), timeint, x);
    x    = xsave(end, :)';

    % Save to traj
    xtraj((iter-1)*nstep+1:iter*nstep,:) = xsave(1:end-1,:);
    ttraj((iter-1)*nstep+1:iter*nstep) = tsave(1:end-1);

    % Update AUV plot
    current_state = stateToQd(x);
    desired_state = Path(time + cstep, current_state);
    QP.UpdateAUVPlot(x, [desired_state.pos; desired_state.vel], time + cstep);
    set(h_title, 'String', sprintf('iteration: %d, time: %4.2f', iter, time + cstep))

    time = time + cstep; % Update simulation time
    t = toc;
    % Check to make sure ode45 is not timing out
    if(t> cstep*50)
        err = 'Ode45 Unstable';
        break;
    end
    % Pause to make real-time
    if real_time && (t < cstep)
        pause(cstep - t);
    end
    % Check termination criteria
    if terminate_check(x, time, stop_pos, pos_tol, vel_tol, max_time)
        break
    end
end
% Truncate xtraj and ttraj
xtraj = xtraj(1:iter*nstep,:);
ttraj = ttraj(1:iter*nstep);

% Truncate saved variables
QP.TruncateHist();
% Plot position
h_pos = figure('Name', ['AUV position']);
plot_state(h_pos, QP.state_hist(1:3,:), QP.time_hist, 'pos', 'vic');
plot_state(h_pos, QP.state_des_hist(1:3,:), QP.time_hist, 'pos', 'des');
legend('Actual Position','Desired Position')

% Assuming QP.state_hist and QP.state_des_hist are matrices of size (3 x num_points)
% and QP.time_hist is a vector containing corresponding timestamps

% Calculate position error between the two sets of data
position_error = QP.state_hist(1:3,:) - QP.state_des_hist(1:3,:);

% Create a figure for the position error plot
h_pos_error = figure('Name', 'Position Error');

% Plot the x-axis position error
subplot(3, 1, 1);
plot(QP.time_hist, position_error(1,:), 'r', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Error [m]');
title('X-axis Position Error');
grid on;

% Plot the y-axis position error
subplot(3, 1, 2);
plot(QP.time_hist, position_error(2,:), 'g', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Error [m]');
title('Y-axis Position Error');
grid on;

% Plot the z-axis position error
subplot(3, 1, 3);
plot(QP.time_hist, position_error(3,:), 'b', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Error [m]');
title('Z-axis Position Error');
grid on;





% Plot velocity
h_vel = figure('Name', ['AUV velocity']);
plot_state(h_vel, QP.state_hist(4:6,:), QP.time_hist, 'vel', 'vic');
plot_state(h_vel, QP.state_des_hist(4:6,:), QP.time_hist, 'vel', 'des');
legend('Actual Velocity','Desired Velocity')


% Assuming QP.state_hist and QP.state_des_hist are matrices of size (3 x num_points)
% and QP.time_hist is a vector containing corresponding timestamps

% Calculate velocity error between the two sets of data
velocity_error = QP.state_hist(4:6,:) - QP.state_des_hist(4:6,:);

% Create a figure for the velocity error plot
h_vel_error = figure('Name', 'Velocity Error');

% Plot the x-axis velocity error
subplot(3, 1, 1);
plot(QP.time_hist, velocity_error(1,:), 'r', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Error [m/s]');
title('X-axis Velocity Error');
grid on;

% Plot the y-axis velocity error
subplot(3, 1, 2);
plot(QP.time_hist, velocity_error(2,:), 'g', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Error [m/s]');
title('Y-axis Velocity Error');
grid on;

% Plot the z-axis velocity error
subplot(3, 1, 3);
plot(QP.time_hist, velocity_error(3,:), 'b', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Error [m/s]');
title('Z-axis Velocity Error');
grid on;





if(~isempty(err))
end
disp('finished.')
t_out = ttraj;
s_out = xtraj;
end