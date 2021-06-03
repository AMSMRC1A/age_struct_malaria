%% Example to test nonlinear fitting of the force of mortality
a = 0:1:60; % ages
b1 = 0;
b2 = 0.05;
b3 = 0.505;
b4 = 0.01;
b5 = 0.055;
y = b1 + b2*exp(-b3*a) + b4*exp(b5*a); % generate synthetic data based on our initial rough fitting
modelfun = @(b,x) b(1) + b(2).*exp(-b(3)*x) + ...
    b(4).*exp(x*b(5)); % define the functional form to fit
beta0 = [0.025 0.1 0.05 0.1 0.1]; % initial guess at the coefficients
X = a;

opts = statset('Display','iter','TolFun',1e-10,'MaxIter',500); % set some options for the fitting
mdl = fitnlm(X,y, modelfun, beta0,'Options',opts); % perform the nonlinear fitting
b_est = mdl.Coefficients{:,1}; % extract the fitted coefficients

% plot the initial data versus the fitted line
plot(a,modelfun([b1 b2 b3 b4 b5]',a),'Linewidth',2);
hold on;
plot(a,modelfun(b_est,a),'Linewidth',2);
legend('true','fitted');







