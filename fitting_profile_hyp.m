clear all
profile_exp = importdata('profile_hyp2.dat');
xdata = profile_exp(:,1)+0.25; ydata = profile_exp(:,2); %center data

plot(xdata,ydata,'Ok','Markersize', 9); hold on

% fitting equation 
% here i fixed some of the parans to reduce the param space
F = @(x,t)0.8636 + x(1).*exp(-2*t./1.6).*sin(pi.*(t-...
    x(2))./0.415*2)./t.^(1/2) + x(3).*exp(-t./1.6).*sin(pi.*(t-...
    x(4))./0.415)./t;

% fitting parameters 
%         A1      xc1        A2       xc2
x0 = [0.4020   0.0823   0.0828    0.0000]; x = x0; % seed
lb = [0, 0, 0, -10, 0, 0, 0, -inf, 0, 0]; % lower bound
ub = [20, 20, 2, 20, 1, 20, 3, 20, 2, 4]; % upper bound

%% Fitting options
curvefitoptions = optimset( 'Display', 'iter' );
[x,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,xdata,ydata,lb,ub,curvefitoptions);

%plot the fitting result
xp = linspace(0.86,xdata(end),1000); yp = F(real(x),xp);
plot(xp,yp, '-b', 'linewidth', 2)
xlabel('L (\mum)')
ylabel('s_3 (a.u.)')