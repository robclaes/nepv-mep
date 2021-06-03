function [x1,l,m,hist]=resinv_symm(W, l0, m0, v1, x10, TOL,f,max_iter)
   if nargin<8
       max_iter = 200
   end
   t0=tic(); 
   T1=@(l,m) -W{1,1}+W{1,2}*l+W{1,3}*m;
   T2=@(l,m) -W{2,1}+W{2,2}*l+W{2,3}*m;
   T1x=@(l,m,x) -W{1,1}*x+l*(W{1,2}*x)+m*(W{1,3}*x);
   T2x=@(l,m,x) -W{2,1}*x+l*(W{2,2}*x)+m*(W{2,3}*x);
   
   hist.resnormnl = NaN(max_iter, 1);
   hist.resnorm=NaN(max_iter,1);
   hist.resnormB=NaN(max_iter,1);
   hist.time_count=NaN(max_iter,1);
   
   dT1 = decomposition(T1(l0,m0));
   dT2 = decomposition(T2(l0,m0));
   T1solve=@(x) dT1\x;
   vT1 = v1'/dT1;
   vT2 = v1'/dT2;
   u1 = vT1*W{1,2};
   u2 = vT2*W{2,2};
   w1 = vT1*W{1,3};
   w2 = vT2*W{2,3};
   
   
   x1=x10/(v1'*x10);
   l=l0;
   m=m0;
   k=1;
   
   while (k==1 || (hist.resnormnl(k-1)> TOL && k<=max_iter))
       g1 = vT1*T1x(l,m,x1);
       g2 = vT2*T2x(l,m,x1);
       AA = [u1*x1 w1*x1;u2*x1 w2*x1];

       dd=-AA\[g1;g2];
       dl=dd(1);
       dm=dd(2);

       l=l+dl;
       m=m+dm;
       x1 = x1-T1solve(T1x(l,m,x1));
       x1=x1/(v1'*x1);
       
       hist.resnormnl(k) = norm(T1x(l,f(x1),x1))/norm(x1);
       hist.resnorm(k)=norm(T1x(l,m,x1))/norm(x1);
       hist.resnormB(k)=norm(T2x(l,m,x1))/norm(x1);
       hist.time_count(k)=toc(t0);
       k=k+1;
   end

end
