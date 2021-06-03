function [x1,x2,x3,l,m1,m2,hist] = resinv3(W,l0,m10,m20,v1,v2,v3,x10,x20,x30,TOL,f1,f2)
    % An implementation of the generalization of residual inverse
    % iteration proposed by Plestenjak in
    % B. Plestenjak: Numerical methods for nonlinear two-parameter eigenvalue
    % problems, BIT Numer. Math. 56 (2016) 241-262
    % Note A{1} is -A{1} in §3 in plestenjaks paper
   t0=tic();
   T1=@(l,m1,m2) -W{1,1}+W{1,2}*l+W{1,3}*m1 + W{1,4}*m2;
   T2=@(l,m1,m2) -W{2,1}+W{2,2}*l+W{2,3}*m1 + W{2,4}*m2;
   T3=@(l,m1,m2) -W{3,1}+W{3,2}*l+W{3,3}*m1 + W{3,4}*m2;
   T1x=@(l,m1,m2,x) -W{1,1}*x+l*(W{1,2}*x)+m1*(W{1,3}*x)+m2*(W{1,4}*x);
   T2x=@(l,m1,m2,x) -W{2,1}*x+l*(W{2,2}*x)+m1*(W{2,3}*x)+m2*(W{2,4}*x);
   T3x=@(l,m1,m2,x) -W{3,1}*x+l*(W{3,2}*x)+m1*(W{3,3}*x)+m2*(W{3,4}*x);
   
   m_prealloc=100;
   hist.resnorm1=NaN(m_prealloc,1);
   hist.resnorm2=NaN(m_prealloc,1);
   hist.resnorm3=NaN(m_prealloc,1);
   hist.time_count=NaN(m_prealloc,1);
   hist.diff1=NaN(m_prealloc,1);
   hist.diff2=NaN(m_prealloc,1);
   hist.diff3=NaN(m_prealloc,1);
   
            
   dT1 = decomposition(T1(l0,m10,m20));
   dT2 = decomposition(T2(l0,m10,m20));
   dT3 = decomposition(T3(l0,m10,m20));
   T1solve=@(x) dT1\x;
   T2solve=@(x) dT2\x;
   T3solve=@(x) dT3\x;

   x1=x10/(v1'*x10);
   x2=x20/(v2'*x20);
   x3=x30/(v3'*x30);
   l=l0;
   m1=m10;
   m2=m20;
   k=1;
   

   while (k==1 || (hist.resnormnl(k-1)> TOL && k<200))
       
       a1=T1solve(W{1,2}*x1);
       a2=T2solve(W{2,2}*x2);
       a3=T3solve(W{3,2}*x3);
       b1=T1solve(W{1,3}*x1);
       b2=T2solve(W{2,3}*x2);
       b3=T3solve(W{3,3}*x3);
       c1=T1solve(W{1,4}*x1);
       c2=T2solve(W{2,4}*x2);
       c3=T3solve(W{3,4}*x3);

       g1=v1'*(T1solve(T1x(l,m1,m2,x1)));
       g2=v2'*(T2solve(T2x(l,m1,m2,x2)));
       g3=v3'*(T3solve(T3x(l,m1,m2,x3)));
       
       AA=[v1'*a1 v1'*b1 v1'*c1; 
           v2'*a2 v2'*b2 v2'*c2;
           v3'*a3 v3'*b3 v3'*c3];

       dd=-AA\[g1;g2;g3];
       dl=dd(1);
       dm1=dd(2);
       dm2=dd(3);
       l=l+dl;
       m1=m1+dm1;
       m2=m2+dm2;
       
       x1=x1-T1solve(T1x(l,m1,m2,x1));
       x2=x2-T2solve(T2x(l,m1,m2,x2));
       x3=x3-T3solve(T3x(l,m1,m2,x3));
       x1=x1/(v1'*x1);
       x2=x2/(v2'*x2);
       x3=x3/(v3'*x3);
       if k==1000
           sdl = 5;
       end
       if k==1200
           sdl = 5;
       end
       if k==1500
           sdl = 5;
       end
       
       hist.diff1(k) = norm(x1-x2);
       hist.diff2(k) = norm(x1-x3);
       hist.diff3(k) = norm(x2-x3);
       hist.resnormnl(k) = norm(T1x(l,f1(x1),f2(x1),x1))/norm(x1);
       hist.resnorm1(k)=norm(T1x(l,m1,m2,x1))/norm(x1);
       hist.resnorm2(k)=norm(T2x(l,m1,m2,x2))/norm(x2);
       hist.resnorm3(k)=norm(T3x(l,m1,m2,x3))/norm(x3);
       hist.time_count(k)=toc(t0);
       
       k=k+1;
   end
end
