%% General Preamble + loading from files
% for sanghyuk


% the files as
% freqs is first row of S3mat.csv
% the y data is row 2:301 of S3mat.csv
% x data is in xs


%S3_profiles=importdata('S3_mat.csv',',',1).Data;
imported=importdata('S3_mat.csv',','); %rn the data has been saved in a csv
S3_profiles=imported(2:301,:); %the data is in rows 2:301
freqs_fit=imported(1,:);%the freqs are the header in row 1
xs=importdata('xs.csv',',');
%% the manual parts that will need to be changed for every line and iteration
%here, select the frequency to fit, plot, and then manually select
freq_num=17; %just which column it is in
%first plot it to decide which parts to use
linefull=S3_profiles(:,freq_num);
freq=freqs_fit(freq_num)
figure();
plot(xs,linefull,'Ok','Markersize',9); 

%% now to get the part to fit
% for column 10, freq: 1397, use the part between 2 and 4
ind_start=80;
ind_end=203;

line2fit=flip(linefull(ind_start:ind_end));
xs2fit=xs(ind_start:ind_end)-(xs(ind_start)-0.25); %consistently start at 0.25

figure();
plot(xs2fit,line2fit,'Ok','Markersize',9); hold on

%% guess some starting parameters
%peak to peak of biggest one is ~0.743 to 1.261, so set offset
%as 0.2591+0.743=1.0021
%the fitting function is:


Lp=1.5;
lambdafit=0.2231;
offset=0.755;
slope=0;
a=1;
st_ind=10; %ignore the first wiggly bit
end_ind=90;

% F_new=@(x,t)offset+x(1).*exp(-2*t./Lp)...
%     .*sin(2*pi.*(t-x(2))./lambdafit)...
%     ./t.^(1/2) + ...
%     x(3).*exp(-t./Lp)...
%     .*sin(pi.*(t-x(4))./lambdafit)./t;

% % varying all the parameters other than offset
F_new=@(x,t)offset+x(1).*exp(-2*t./x(5))...
    .*sin(2*pi.*(t-x(2))./lambdafit)...
    ./t.^(1/2) + ...
    x(3).*exp(-t./x(5))...
    .*sin(pi.*(t-x(4))./lambdafit)./(t); ...
    %+slope.*t; %maybe a slope???

%try adding an 'a' term as a decay
%force slope and lambda fit
% F_new=@(x,t)offset+x(1).*exp(-2*t./x(5))...
%     .*sin(2*pi.*(t-x(2))./lambdafit)...
%     ./t.^(1/2) + ...
%     x(3).*exp(-t./x(5))...
%     .*sin(pi.*(t-x(4))./lambdafit)./t...
%     +slope.*t; %maybe a slope???



%% fitting parameters 

%      A1  xc1   A2   xc2    Lp lambdafit offset a slope
x0 = [1, 1, 1, 1, Lp, lambdafit, offset, a, slope]; x = x0; % seed
lb = [0, 0,0, 0, 0, 0,0,0, 0,-1]; % lower bound
ub = [3, 15, 5, 5, 3, 2,2, 5,1 ]; % upper bound


%% Fitting options
curvefitoptions = optimset( 'Display', 'final' ,'MaxFunEvals',2e5,'TolFun',1e-8);
[x,resnorm,~,exitflag,output] = lsqcurvefit(F_new,x0,xs2fit(st_ind:end_ind),line2fit(st_ind:end_ind),lb,ub,curvefitoptions);


%plot the fitting result
xp = linspace(xs2fit(st_ind),xs2fit(end_ind),1000); yp = F_new(real(x),xp); plot(xp,yp, '-', 'linewidth', 2)
xlabel('L (\mum)')
ylabel('s_3 (a.u.)')  


