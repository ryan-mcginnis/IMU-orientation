function qinf = get_orientation_optim_quaternion(time, a, w, ind)
%Function to determine orientation of IMU.  Does so in two steps: 1) Define 
%initial orientation of device based on direction of gravity and 2) Define
%orientation thereafter by fusing acceleration and angular velocity
%estimates via optimization.  Method assumes that the device is initially 
%at rest.
%Written by Ryan S. McGinnis - ryan.mcginnis14@gmail.com - June 25, 2016
%
% Copyright (C) 2016  Ryan S. McGinnis
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% Method introduced and validated in:
% McGinnis, R. S., et al. "Accuracy of Femur Angles Estimated by IMUs 
% During Clinical Procedures Used to Diagnose Femoroacetabular 
% Impingement." IEEE transactions on bio-medical engineering 62.6 (2015): 
% 1503-1513.
%
%Inputs:
%1. time - time (s, nx1) 
%2. a - acceleration (g, nx3) 
%3. w - angular velocity (deg/s, nx3) 
%4. ind - index array specifying the samples where the device is at rest 
%         at the beginning of the trial (binary, nx1)
%
%Outputs:
%1. qinf - quaternion describing IMU orientation (nx4)


% Define signal characteristics based on initial still section
amag = sqrt(sum(a.^2,2));
wmag = sqrt(sum(w.^2,2));
g_mag = mean(amag(ind));
w0_mag = mean(wmag(ind)); %deg/s
a_noise = std(amag(ind));
w_noise = std(wmag(ind)); %deg/s
a_lb = -2.5*a_noise; %-2.5 * sd of noise
a_ub = -a_lb; %2.5 * sd of noise
w_ub = w0_mag + 5 * w_noise; %5 * sd above mean


% Calculate gravity direction from acceleration at still instances
indG = amag < a_ub+g_mag & amag > a_lb+g_mag & wmag < w_ub;
aG = a(indG,:); 
aGnorm = sqrt(sum(aG.^2,2)); 
aG = [aG(:,1)./aGnorm, aG(:,2)./aGnorm, aG(:,3)./aGnorm];


% Define parameters to pass into optimization
params.indG = indG;
params.ind = ind;
params.direction = [0, 0, 1]; %direction of gravity in world frame
params.a = a;


% Solve for gyro scale factors and biases that min. sum of squared error
options = optimset('maxiter',300,'tolfun',1e-4,...
                   'tolx',1e-6,'MaxFunEvals',5000,...
                   'algorithm',{'levenberg-marquardt',.01},...
                   'plotfcns',{@optimplotx,@optimplotresnorm}); 

f = @(x)SF_calc(x,aG,time,w,params);
x0 = [0, 0, 0, 0, 0, 0]; %[sfx, sfy, sfy, bx, by, bz]
SF = lsqnonlin(f,x0,[],[],options);


%Calculate new orientation
wnew = w * diag([1+SF(1), 1+SF(2), 1+SF(3)]) - (w(:,1).^0)*[SF(4), SF(5), SF(6)];
qinf = get_orientation_quaternion(time, a, wnew, params.ind);

end

%Function to define objective function for SF and bias optimization-------%
function e = SF_calc(x,aG,t,w,params)

%Correct angular velocity with latest guess
wnew = w * diag([1+x(1), 1+x(2), 1+x(3)]) - (w(:,1).^0)*[x(4), x(5), x(6)];

%Define orientation as func of time and use to estimate direction of grav
qUF = get_orientation_quaternion(t, params.a, wnew, params.ind);
wG = quaternRot(quaternConj(qUF(params.indG,:)), params.direction); %gravity direction in body frame

%Define error vector
e = real(acos(dot(wG,aG,2))); %Minimize angle between predictions

%Weight errors by time so errors at end of trial weighted more heavily
e = e .* t(params.indG);

end