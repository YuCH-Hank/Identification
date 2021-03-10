clc ; clear ; close all ;

%=============================== Parameters ===============================
[ dt , tf , T , g , fs , wf , qb , n ] = parameter () ;

% fs Fourier series
% n  axis
Acol = sum (2 * fs + 1) ;

%=============================== A * x <= B ===============================
A = []; B = [];
for t = 0 : dt : tf
    A3 = [] ; B1 = [];
    for i = 1 : n
        A2 = [] ;
        for l = 1 : fs(i)
            A1 = [  sin(wf*l*t)/(wf*l)   -cos(wf*l*t)/(wf*l)
                   -sin(wf*l*t)/(wf*l)    cos(wf*l*t)/(wf*l)   
                       cos(wf*l*t)           sin(wf*l*t)
                      -cos(wf*l*t)          -sin(wf*l*t)
                    -wf*l*sin(wf*l*t)      wf*l*cos(wf*l*t)
                     wf*l*sin(wf*l*t)     -wf*l*cos(wf*l*t)  ] ;
            A2 = [A2, A1];
        end
        Aq0 = [ 1 -1  0  0  0  0 ]' ;
        A2  = [A2, Aq0];
        A3  = [A3, zeros(size(A3,1), size(A2,2)); zeros(size(A2,1), size(A3,2)), A2 ] ;
        B1  = [B1; qb(i,1) ; -qb(i,2) ; qb(i,3) ; qb(i,3) ; qb(i,4) ; qb(i,4) ] ;
    end
    A = [A; A3 ];
    B = [B; B1 ];
end
clear A1 A2 A3 B1 Aq0 t i l

%============================= Aeq * x = Beq ==============================  
Aeq = [] ;
for i = 1 : n
    Aeq2 = [];
    for l = 1 : fs(i)
        Aeq1 = [  sin(wf*l*tf)/(wf*l)  -(cos(wf*l*tf)-1)/(wf*l)
                    cos(wf*l*tf)-1           sin(wf*l*tf)
                  -wf*l*sin(wf*l*tf)     wf*l*(cos(wf*l*tf)-1)
                        0                   -1/(wf*l)
                        1                       0
                        0                      wf*l           ] ;
        Aeq2 = [Aeq2, Aeq1] ;
    end
    Aeqq0 = [ 0 0 0 1 0 0]' ;
    Aeq2 = [Aeq2, Aeqq0];
    Aeq  = [Aeq, zeros(size(Aeq,1), size(Aeq2,2)); zeros(size(Aeq2,1), size(Aeq,2)), Aeq2 ] ;
end
Beq = zeros( 6 * n , 1 ) ;
clear Aeq1 Aeq2 Aeqq0 i l

%=============================== fmincon ================================== 
X0 = 1000000 ; % X0大小與圖形變化成正比,太大將造成圖形多為弦波(挑選時由小往上加)
x0 = rand( size ( A , 2 ) , 1 ) * X0 - ( X0 / 2 ) ; % -X0/2 ~ X0/2

options = optimoptions('fmincon','Algorithm','sqp','Display','iter') ;
[ x , cond_A ] = fmincon( @optimfun , x0 , A , B , Aeq , Beq , [] , [] , []      , options ) ;
%                fmincon( FUN,        X  , A,  B , Aeq , Beq , LB , UB , NONLCON , options , varargin)
%   FMINCON attempts to solve problems of the form:
%    min F(X)  subject to:  A*X  <= B, Aeq*X  = Beq (linear constraints)

clear A B Aeq Beq x0
Q = xfun(x) ;

%% ================================ Desire ==================================
period = 10;
% 用 0.1 算 , 再以 0.001內插
Tinterp = 0 : dt / 100 : tf ;

%--------------------------------- PosDes ---------------------------------
PosDesx1 = Q ( 1 : n : size( Q , 1 ) , 1 ) ;
PosDesy1 = Q ( 2 : n : size( Q , 1 ) , 1 ) ;

PosDesx2 = interp1( T , PosDesx1 , Tinterp , 'spline' ) ;
PosDesy2 = interp1( T , PosDesy1 , Tinterp , 'spline' ) ;

%--------------------------------- VelDes ---------------------------------
VelDesx1 = Q ( 1 : n : size( Q , 1 ) , 2 ) ;
VelDesy1 = Q ( 2 : n : size( Q , 1 ) , 2 ) ;

VelDesx2 = interp1( T , VelDesx1 , Tinterp , 'spline' ) ;
VelDesy2 = interp1( T , VelDesy1 , Tinterp , 'spline' ) ;

%--------------------------------- AccCmd ---------------------------------
AccDesx1 = Q ( 1 : n : size( Q , 1 ) , 3 ) ;
AccDesy1 = Q ( 2 : n : size( Q , 1 ) , 3 ) ;

AccDesx2 = interp1( T , AccDesx1 , Tinterp , 'spline' ) ;
AccDesy2 = interp1( T , AccDesy1 , Tinterp , 'spline' ) ;

clear PosDesx1 PosDesy1 VelDesx1 VelDesy1 AccDesx1 AccDesy1

num = size( PosDesx2 , 2 ) ;
%--------------------------------- 10 個週期 ------------------------------
for i = 1 : period
    PosDesX( num * ( i - 1 ) + 1 : num * i , 1 ) = PosDesx2 ;
    PosDesY( num * ( i - 1 ) + 1 : num * i , 1 ) = PosDesy2 ;
    VelDesX( num * ( i - 1 ) + 1 : num * i , 1 ) = VelDesx2 ;
    VelDesY( num * ( i - 1 ) + 1 : num * i , 1 ) = VelDesy2 ;
    AccDesX( num * ( i - 1 ) + 1 : num * i , 1 ) = AccDesx2 ;
    AccDesY( num * ( i - 1 ) + 1 : num * i , 1 ) = AccDesy2 ;
end

clear PosDesx2 PosDesy2 VelDesx2 VelDesy2 AccDesx2 AccDesy2 period

time = linspace( 0, length(PosDesX(:,1)) - 1,  length(PosDesX(:,1)))' * 0.001 ;

%--------------------------------- txt ------------------------------------
DesireData = [ time PosDesX  PosDesY  VelDesX  VelDesY  AccDesX  AccDesY ] ;

Generate_txt( DesireData, 'Trajectory');
Generate_txt( DesireData(1:size( PosDesX , 1 )/5, : ), 'Trajectory_Test');

%============================= plot results ===============================
figure(1)
for i = 1 : n
    subplot( n , 1 , i )
    plot( T , Q ( i : n : size( Q , 1 ) , 1 ) , '-' , ...
          T , Q ( i : n : size( Q , 1 ) , 2 ) , '--' , ...
          T , Q ( i : n : size ( Q , 1 ) , 3 ) , '-.' ) ;
    title( [ 'Axis ' , num2str(i) ] , 'FontWeight' , 'bold' , 'FontSize' , 12 ) ;
    xlabel('Time (s)') ; ylabel('Pos (rad) , Vel (rad/s) , Acc (rad/s^2)', 'FontSize' , 10 ) ;
    legend( [ 'Pos [' , num2str(max(Q(i:n:size(Q,1),1))) , '  ' , num2str(min(Q(i:n:size(Q,1),1))) , ']' ] , ...
            [ 'Vel [' , num2str(max(Q(i:n:size(Q,1),2))) , '  ' , num2str(min(Q(i:n:size(Q,1),2))) , ']' ] , ...
            [ 'Acc [' , num2str(max(Q(i:n:size(Q,1),3))) , '  ' , num2str(min(Q(i:n:size(Q,1),3))) , ']' ] ) ;
    grid on ;
end

J1rpm = max( abs ( Q( 1 : n : size( Q , 1 ) , 2 ) ) ) * ( 60 / ( 2 * pi ) ) * 50
J2rpm = max( abs ( Q( 2 : n : size( Q , 1 ) , 2 ) ) ) * ( 60 / ( 2 * pi ) ) * 50

figure(2)
plot( T , Q ( 1 : n : size( Q , 1 ) , 1 ) , 'k-' ) ;
title( 'Axis 1' , 'FontWeight' , 'bold' , 'FontSize' , 12 ) ;
xlabel('Time (s)') ; ylabel('Pos (rad)') ;
axis( [ 0 , tf , qb(1,2) , qb(1,1) ] );
grid on ;

figure(3)
plot( T , Q ( 2 : n : size( Q , 1 ) , 1 ) , 'k-' ) ;
title( 'Axis 2' , 'FontWeight' , 'bold' , 'FontSize' , 12 ) ;
xlabel('Time (s)') ; ylabel('Pos (rad)') ;
axis( [  0 , tf , qb(2,2) , qb(2,1) ] );
grid on ;



        
