function [ Q ] = xfun ( x )

[ dt , tf , T , g , fs , wf , qb , n ] = parameter () ;

d = 1 ;
for t = 0 : dt : tf
    c = 0 ;
	q = zeros ( n , 3 ) ;
    for i = 1 : n
        for l = 1 : fs(i)
            a = ( 2 * l - 1 ) + c ; b = ( 2 * l ) + c ;
            pos = ( x(a) / ( wf * l ) ) * sin( wf * l * t ) - ( x(b) / ( wf * l ) ) * cos( wf * l * t ) ;
            vel =   x(a) * cos( wf * l * t )                +   x(b) * sin( wf * l * t ) ;
            acc = - x(a) * wf * l * sin( wf * l * t )       +   x(b) * wf * l * cos( wf * l * t ) ;
            q(i,1) = q(i,1) + pos ;
            q(i,2) = q(i,2) + vel ;
            q(i,3) = q(i,3) + acc ;
        end
        c = c + ( 2 * fs(i) + 1 ) ;
        q(i,1) = q(i,1) + x(c) ;
    end
    
    Q( d : d + ( n - 1 ) , : ) = q ;
    d = d + n ;
end





