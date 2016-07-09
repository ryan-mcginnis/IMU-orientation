function Rinf = get_orientation_compfilter(time, a, w, ind)
%Function to determine orientation of IMU.  Does so in two steps: 1) Define 
%initial orientation of device based on direction of gravity and 2) Define
%orientation thereafter by fusing acceleration and angular velocity
%estimates via complementary filtering.  Method assumes that the device is 
%initially at rest.
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
% Method validated in:
% McGinnis RS, Cain SM, Davidson SP, Vitali RV, McLean SG, Perkins NC. 
% Validation of Complementary Filter Based IMU Data Fusion for Tracking 
% Torso Angle and Rifle Orientation. In ASME 2014 International Mechanical 
% Engineering Congress and Exposition 2014 Nov 14 (pp. V003T03A052-V003T03A052). 
% American Society of Mechanical Engineers.
%
%Inputs:
%1. time - time (s, nx1) 
%2. a - acceleration (g, nx3) 
%3. w - angular velocity (deg/s, nx3) 
%4. ind - index array specifying the samples where the device is at rest 
%         at the beginning of the trial (binary, nx1)
%
%Outputs:
%1. Rinf - DCM describing IMU orientation (3x3xn)

%Define signal characteristics based on still section
amag = sqrt(sum(a.^2,2));
wmag = sqrt(sum(w.^2,2));
g_mag = mean(amag(ind));
w0_mag = mean(wmag(ind)); %deg/s
a_noise = std(amag(ind));
w_noise = std(wmag(ind)); %deg/s
a_lb = -3*a_noise; %-3 * sd of noise
a_ub = -a_lb; %3 * sd of noise
w_ub = (w0_mag + 3 * w_noise)*pi/180; %3 * sd above mean

%Define initial Z direction (gravity)
Z = mean(a(ind,:)) ./ norm(mean(a(ind,:)));

%Define remaining axes 
%X in terms of (i,j,k) according to X = j x Z
X = cross([0,0,1],Z); X = X./norm(X);

%Re-define Y in terms of (i,j,k) to ensure orthogonality (Y = Z x X)
Y = cross(Z,X); Y = Y./norm(Y);

%Define initial DCM    
R = [X; Y; Z];

%Define nominal filter gains
Kp0 = 2;
Ki0 = 1;

%Initialize variables
Rinf = zeros(3,3,length(time)); Rinf(:,:,1) = R;
dc = zeros(size(time,1),3);

%Convert angular velocity to rad/s and remove bias
w = (w - (w(:,1).^0)*mean(w(ind,:))) * pi/180;
wh = w;

for t=2:length(time)
    %Propegate orientation
    dt = time(t) - time(t-1);
    tt = 0.5 * dt * (w(t-1,:) + w(t,:));
    Rt = Rinf(:,:,t-1)*(eye(3)+(2/(1+0.5*(tt*tt.')))*...
        (0.5*skew(tt) + 0.25*skew(tt)^2));
    
    %Assign Gain Values
    amag = norm(a(t,:));
    wmag = norm(wh(t,:));
    if amag < a_ub+g_mag && amag > a_lb+g_mag && wmag < w_ub
        Ki = Ki0;
        Kp = Kp0;
    elseif amag < a_ub+g_mag && amag > a_lb+g_mag && wmag < 5*w_ub
        Ki = Ki0/10;
        Kp = Kp0/10;
    else
        Ki = 0;
        Kp = 0;
    end
    
    %Estimate error relative to gravity direction
    at = a(t,:)./norm(a(t,:));
    we = cross(at, Rt(3,:));
    
    %Adjust angular velocity
    dc(t,:) = dc(t-1, :) + dt * (-Ki * we);
    wh(t,:) = w(t,:) - dc(t,:) + Kp * we;
    
    %Calculate corrected orientation
    theta = 0.5 * dt * (wh(t-1,:) + wh(t,:));
    Rinf(:,:,t) = Rinf(:,:,t-1)*(eye(3)+(2/(1+0.5*(theta*theta.')))*...
                      (0.5*skew(theta) + 0.25*skew(theta)^2));
end
end