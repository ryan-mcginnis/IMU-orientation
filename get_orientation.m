function Rinf = get_orientation(time, a, w, ind)
%Function to determine orientation of IMU.  Does so in two steps: 1) Define 
%initial orientation of device based on direction of gravity and 2) Define
%orientation thereafter by fusing acceleration and angular velocity
%estimates.  Method assumes that the device is initially at rest.
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
% Method introduced in:
% McGinnis, Ryan S., and Noel C. Perkins. "A Highly Miniaturized, Wireless 
% Inertial Measurement Unit for Characterizing the Dynamics of Pitched 
% Baseballs and Softballs." Sensors 12 (2012): 11933-11945.
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

%Define initial Z direction (gravity)
Z = mean(a(ind,:)) ./ norm(mean(a(ind,:)));

%Define remaining axes 
%X in terms of (i,j,k) according to X = j x Z
X = cross([0,0,1],Z); X = X./norm(X);

%Re-define Y in terms of (i,j,k) to ensure orthogonality (Y = Z x X)
Y = cross(Z,X); Y = Y./norm(Y);

%Define initial DCM    
R = [X; Y; Z];

%Initialize variables
Rinf = zeros(3,3,length(time)); Rinf(:,:,1) = R;

%Convert angular velocity to rad/s and remove bias
w = (w - (w(:,1).^0)*mean(w(ind,:))) * pi/180;

for t=2:length(time)
    dt = time(t) - time(t-1);
    tt = 0.5 * dt * (w(t-1,:) + w(t,:));
    Rinf(:,:,t) = Rinf(:,:,t-1)*(eye(3)+(2/(1+0.5*(tt*tt.')))*...
        (0.5*skew(tt) + 0.25*skew(tt)^2));
end
end