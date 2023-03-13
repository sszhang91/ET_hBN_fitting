w=transpose(500:0.1:8000);
method='makima';

wb=W;
eb1=E1;
eb2=E2;


wb=nonzeros(wb);
indexb=find(wb);

eps_b1=interp1(wb,eb1(indexb),w,method,'extrap');
eps_b2=interp1(wb,eb2(indexb),w,method,'extrap');


figure; hold;
plot(w,eps_b1);
plot(W, E1,'o');
plot(w,eps_b2);
plot(W,E2,'x');

