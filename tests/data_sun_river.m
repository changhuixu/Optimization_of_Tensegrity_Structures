% Sun River tensegrity structure data
global tube_pts cable_pts fixed_pts tube_lengths tube_masses cable_len;
global e3 g Es;

% This matrix contains indexes of points that tubes join.
% The start index of the i'th tube is tube_pts(i,1), 
% 				 and the end point is tube_pts(i,2)
tube_pts = [1,2; 
			3,4; 
			5,6; 
			7,8];

% This matrix contains indexes of points that cables join.
% Same format as tube_pts.
cable_pts = [1,7;
			 1,8;
			 2,3;
			 2,4;
			 2,5;
			 2,6;
			 3,7;
			 3,8;
			 4,7;
			 4,8;
			 5,7;
			 5,8;
			 6,7;
			 6,8];

% The points with these indexes are fixed
fixed_pts = [1;3;5];

% These are the tube lengths (in meters)
tube_lengths = [3.0000;
				2.2000;
				2.2000;
				2.4000];

% The masses/length are in kg per meter.
% tube_mass_per_length =

%    0.4800
%    0.2800
%    0.2800
%    0.2800

% Combining the above information, we get:
tube_masses = [1.4400;
			   0.6160;
			   0.6160;
			   0.6720];

% Point i has position x(:,i) = the i'th column (coordinates in meters)
x = [ -1.2000,1.0981,    0,    -0.5694,  0,    -0.5694, -1.0000, -1.0000;
        0,      0,    -0.7000, -0.7000, 0.7000, 0.7000,  1.2000, -1.2000;
        0,    2.5981,    0,     2.1250,  0,     2.1250,  0.7000,  0.7000];

% Cable lengths in meters
cable_len =[1.4036;
			1.4036;
			2.9062;
			1.8693;
			2.9062;
			1.8693;
			2.2583;
			1.3191;
			2.4137;
			1.5704;
			1.3191;
			2.2583;
			1.5704;
			2.4137];

g = 9.82;
Es = 2.5e6;
e3 = [0;0;1];


