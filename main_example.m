%% IMU Orientation Main Example   
% This script provides example implementations of IMU orientaiton functions  
% Written by Ryan S. McGinnis - ryan.mcginis14@gmail.com - July 9, 2016

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

%% Load example IMU data
load('example_data.mat');
time = data.t;
a = data.a; %g
w = data.w; %deg/s
ind = data.ind_still; %id time period at start of trial wher subject still


%% Calculate orientation using DCM 
Rinf = get_orientation(time, a, w, ind);
a_dcm = dcmRot(Rinf,a);
w_dcm = dcmRot(Rinf,w);


%% Calculate orientation using quaternion 
qinf = get_orientation_quaternion(time, a, w, ind);
a_quatern = quaternRot(qinf,a);
w_quatern = quaternRot(qinf,w);


%% Calculate orientation using DCM and complementary filtering approach
Rinf_comp = get_orientation_compfilter(time, a, w, ind);
a_dcm_comp = dcmRot(Rinf_comp,a);
w_dcm_comp = dcmRot(Rinf_comp,w);


%% Calculate orientation using quaternion and complementary filtering approach
qinf_comp = get_orientation_compfilter_quaternion(time, a, w, ind);
a_quatern_comp = quaternRot(qinf_comp,a);
w_quatern_comp = quaternRot(qinf_comp,w);


%% Plot the results
% Plots show vertical component of world-fixed accelerometer measurment
figure; 
set(gcf,'name','dcm vs quaternion (red) implementations');
hold on;
plot(time,a_dcm(:,3));
plot(time,a_quatern(:,3),'r');
xlabel('Time (s)');
ylabel('Acceleration (g)');

figure; 
set(gcf,'name','dcm vs quaternion (red) implementations of complementary filter');
hold on;
plot(time,a_dcm_comp(:,3));
plot(time,a_quatern_comp(:,3),'r');
xlabel('Time (s)');
ylabel('Acceleration (g)');

figure; 
set(gcf,'name','dcm raw vs complementary filter (red)');
hold on;
plot(time,a_dcm(:,3));
plot(time,a_dcm_comp(:,3),'r');
xlabel('Time (s)');
ylabel('Acceleration (g)');

