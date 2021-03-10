function [ dt , tf , T , g , fs , wf , qb , n ] = parameter ()

dt = 0.1 ; tf = 10 ; T = 0 : dt : tf ;
g = 9.81 ;
wf = 0.15 * pi ;  % fundamental pulsation [ 0.15 ~ 0.05 ]
fs = [ 5 5 ] ;  % Fourier series

qb = [ 90*pi/180  -90*pi/180  5  7.0
       90*pi/180  -90*pi/180  5  7.0 ] ;  % qb = [ q_max  q_min  qd1_max  qd2_max ]

% qb = [ 80*pi/180  -80*pi/180  1.5  2
%        80*pi/180  -80*pi/180  1.5  2 ] ;  % qb = [ q_max  q_min  qd1_max  qd2_max ]

  
n = size( qb , 1 ) ; % number of axes

