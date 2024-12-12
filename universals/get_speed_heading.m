
function [ speed , heading ] = get_speed_heading(x, y, Hz)

%% get_speed_heading.m 
% Produce speed and heading time series from time series of position (x,y). 

% Written by Greg Dachner for VENLAB 2017
% gdachner@gmail.com
% Changed from Kevin Rio's code
% Updated 10/19/2018 to add Hz

vx = diff(x);
vy = diff(y);

speed  = sqrt( vx.^2 + vy.^2 );
heading = cart2pol(diff(x),diff(y));

speed = vertcat(speed(1,:),speed);
heading = vertcat(heading(1,:),heading);

% Convert to proper units
speed = speed * Hz; 
heading = -rad2deg(heading-deg2rad(90));

% Rotate polar to cartesian coordinats so they equal straight head / upward
% direction as 0 degrees
for j = 1:size(heading,2)
    for i = 1:size(heading,1)
         if heading(i,j) > 180
             heading(i,j) = heading(i,j) - 360;
         end
     end
end

