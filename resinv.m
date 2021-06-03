function [x1,x2,l,m,hist] = resinv(W, l0, m0, v1, v2, x10, x20, TOL,f)
    % An implementation of the generalization of residual inverse
    % iteration proposed by Plestenjak in
    % B. Plestenjak: Numerical methods for nonlinear two-parameter eigenvalue
    % problems, BIT Numer. Math. 56 (2016) 241-262
    % Note A{1} is -A{1} in §3 in plestenjaks paper
   t0=tic();
   T1=@(l,m) -W{1,1}+W{1,2}*l+W{1,3}*m;
   T2=@(l,m) -W{2,1}+W{2,2}*l+W{2,3}*m;
   T1x=@(l,m,x) -W{1,1}*x+l*(W{1,2}*x)+m*(W{1,3}*x);
   T2x=@(l,m,x) -W{2,1}*x+l*(W{2,2}*x)+m*(W{2,3}*x);
   
   m_prealloc=100;
   hist.resnormnl = NaN(m_prealloc,1);
   hist.resnorm=NaN(m_prealloc,1);
   hist.resnormB=NaN(m_prealloc,1);
   hist.time_count=NaN(m_prealloc,1);
   hist.equalnorm=NaN(m_prealloc,1);
   
            
   dT1 = decomposition(T1(l0,m0));
   dT2 = decomposition(T2(l0,m0));
   T1solve=@(x) dT1\x;
   T2solve=@(x) dT2\x;

   x1=x10/(v1'*x10);
   x2=x20/(v2'*x20);
   l=l0;
   m=m0;
   k=1;
   

   while (k==1 || (hist.resnormnl(k-1)> TOL && k<200))
       a1=T1solve(W{1,2}*x1);
       a2=T2solve(W{2,2}*x2);
       b1=T1solve(W{1,3}*x1);
       b2=T2solve(W{2,3}*x2);

       g1=v1'*(T1solve(T1x(l,m,x1)));
       g2=v2'*(T2solve(T2x(l,m,x2)));
       AA=[v1'*a1 v1'*b1 ; v2'*a2 v2'*b2 ];

       dd=-AA\[g1;g2];
       dl=dd(1);
       dm=dd(2);
       l=l+dl;
       m=m+dm;
       
       x1=x1-T1solve(T1x(l,m,x1));
       x2=x2-T2solve(T2x(l,m,x2));
       x1=x1/(v1'*x1);
       x2=x2/(v2'*x2);
        
     
       hist.resnormnl(k) = norm(T1x(l,f(x1),x1))/norm(x1);
       hist.resnorm(k)=norm(T1x(l,m,x1))/norm(x1);
       hist.resnormB(k)=norm(T2x(l,m,x2))/norm(x1);
       hist.time_count(k)=toc(t0);
       hist.equalnorm(k)=norm(x1-x2)/norm(x1)/norm(x2);
       
       k=k+1;
   end
end
