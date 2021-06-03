function [x1,l,m1,m2,hist]=resinv_symm3(W, l0, m10, m20, v1, x10, TOL,f1,f2)
   t0=tic();
   T1=@(l,m1,m2) -W{1,1}+W{1,2}*l+W{1,3}*m1 + W{1,4}*m2;
   T2=@(l,m1,m2) -W{2,1}+W{2,2}*l+W{2,3}*m1 + W{2,4}*m2;
   T3=@(l,m1,m2) -W{3,1}+W{3,2}*l+W{3,3}*m1 + W{3,4}*m2;
   T1x=@(l,m1,m2,x) -W{1,1}*x+l*(W{1,2}*x)+m1*(W{1,3}*x)+m2*(W{1,4}*x);
   T2x=@(l,m1,m2,x) -W{2,1}*x+l*(W{2,2}*x)+m1*(W{2,3}*x)+m2*(W{2,4}*x);
   T3x=@(l,m1,m2,x) -W{3,1}*x+l*(W{3,2}*x)+m1*(W{3,3}*x)+m2*(W{3,4}*x);
   
   m_prealloc=100;
   hist.resnormnl = NaN(m_prealloc, 1);
   hist.resnorm1=NaN(m_prealloc,1);
   hist.resnorm2=NaN(m_prealloc,1);
   hist.resnorm3=NaN(m_prealloc,1);
   hist.time_count=NaN(m_prealloc,1);
   
   dT1 = decomposition(T1(l0,m10,m20));
   dT2 = decomposition(T2(l0,m10,m20));
   dT3 = decomposition(T3(l0,m10,m20));
   T1solve=@(x) dT1\x;
   
   vT1 = v1'/dT1;
   vT2 = v1'/dT2;
   vT3 = v1'/dT3;
   u1 = vT1*W{1,2};
   u2 = vT2*W{2,2};
   u3 = vT3*W{3,2};
   w1 = vT1*W{1,3};
   w2 = vT2*W{2,3};
   w3 = vT3*W{3,3};
   z1 = vT1*W{1,4};
   z2 = vT2*W{2,4};
   z3 = vT3*W{3,4};
   
   x1=x10/(v1'*x10);
   l=l0;
   m1=m10;
   m2=m20;
   k=1;
   
   while (k==1 || (hist.resnormnl(k-1)> TOL && k<200))
       g1 = vT1*T1x(l,m1,m2,x1);
       g2 = vT2*T2x(l,m1,m2,x1);
       g3 = vT3*T3x(l,m1,m2,x1);
       AA = [u1*x1 w1*x1 z1*x1;
             u2*x1 w2*x1 z2*x1;
             u3*x1 w3*x1 z3*x1];

       dd=-AA\[g1;g2;g3];
       dl=dd(1);
       dm1=dd(2);
       dm2=dd(3);

       l=l+dl;
       m1=m1+dm1;
       m2=m2+dm2;
       x1 = x1-T1solve(T1x(l,m1,m2,x1));
       x1=x1/(v1'*x1);
       
       hist.resnormnl(k)=norm(T1x(l,f1(x1),f2(x1),x1))/norm(x1);
       hist.resnorm1(k)=norm(T1x(l,m1,m2,x1))/norm(x1);
       hist.resnorm2(k)=norm(T2x(l,m1,m2,x1))/norm(x1);
       hist.resnorm3(k)=norm(T3x(l,m1,m2,x1))/norm(x1);
       hist.time_count(k)=toc(t0);
       k=k+1;
   end

end
