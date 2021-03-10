function [ cond_A ] = optimfun ( x )

[ dt , tf , T , g , fs , wf , qb , n ] = parameter () ;

d = 1 ;
for t = 0 : dt : tf
    c = 0 ;
	q = zeros ( n , 3 ) ;
    for i = 1 : n
        for l = 1 : fs(i)
            a = ( 2 * l - 1 ) + c ; b = ( 2 * l ) + c ;
            pos = ( x(a) / ( wf * l ) ) * sin( wf * l * t ) - ( x(b) / ( wf * l ) ) * cos( wf * l * t ) ;
            vel = x(a) * cos( wf * l * t ) + x(b) * sin( wf * l * t ) ;
            acc = - x(a) * wf * l * sin( wf * l * t ) + x(b) * wf * l * cos( wf * l * t ) ;
            q(i,1) = q(i,1) + pos ;
            q(i,2) = q(i,2) + vel ;
            q(i,3) = q(i,3) + acc ;
        end
        c = c + ( 2 * fs(i) + 1 ) ;
        q(i,1) = q(i,1) + x(c) ;
    end
    
    % Row: number of axes, Column: 1 Pos, 2 Vel, 3 Acc

%     b1 = 0.24 ;
%     fi = [ q(1,3)  q(1,3)+q(2,3)  b1*((2*q(1,3)+q(2,3))*cos(q(2,1))-q(2,2)*(2*q(1,2)+q(2,2))*sin(q(2,1)))  b1*(-(2*q(1,3)+q(2,3))*sin(q(2,1))-q(2,2)*(2*q(1,2)+q(2,2))*cos(q(2,1)))   q(1,2)      0     sign(q(1,2))       0     ;
%               0    q(1,3)+q(2,3)           b1*((q(1,2)^2)*sin(q(2,1))+q(1,3)*cos(q(2,1)))                             b1*((q(1,2)^2)*cos(q(2,1))-q(1,3)*sin(q(2,1)))                    0     q(2,2)        0        sign(q(2,2)) ] ;
   
    fi = [ q(1,3)  q(1,3)+q(2,3)  (2*q(1,3)+q(2,3))*cos(q(2,1))-q(2,2)*(2*q(1,2)+q(2,2))*sin(q(2,1))   q(1,2)      0     sign(q(1,2))       0          ;
              0    q(1,3)+q(2,3)           (q(1,2)^2)*sin(q(2,1))+q(1,3)*cos(q(2,1))                     0       q(2,2)        0        sign(q(2,2)) ] ;
      


    A( d : d + ( n - 1 ) , : ) = fi ;
    d = d + n ;
end

for i = 1 : size ( A , 2 )
	norm_A( : , i ) = A( : , i ) ./ norm( A( : , i ) ) ;
end
cond_A = cond( norm_A ) ;


