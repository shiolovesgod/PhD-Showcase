
function [curve, gof] = circFit(xz, sNo, weights)

opts=fitoptions('method','nonlinear');
opts.algorithm='trust-region';
opts.maxiter=5000;
opts.maxfunevals=5000;
opts.upper=[inf, 30, inf]; %R, x0, y0,

x= xz(:,1); z = xz(:,2);

%R, x0, y0, z0 --> order of the variables for setting limits
switch sNo
    case {1, 2, 3} %anterior cornea and lens
        %Anterior Surface
        
        opts.lower=[-inf, -10, -inf];
        opts.startpoint=[7, mean(x), mean(z)]; %7 lens; starting point R, 7; 1mm z0, 1.8mm (posterior cornea)
       fitEqn = '(z0+R)-sqrt(R^2-(x-x0)^2)';
       fitFcn = @(b,x) (b(3)+b(1))-sqrt(b(1)^2-(x-b(2)).^2);
    case 4 %lens only
        %Posterior Surface
                         %R,   x0,  y0
        opts.lower=[-inf, -10, -inf];
        opts.startpoint=[7, mean(x), mean(z)];
        fitEqn = '(z0+R)+sqrt(R^2-(x-x0)^2)';
        fitFcn = @(b,x) (b(3)+b(1))+sqrt(b(1)^2-(x-b(2)).^2);
end


if nargin < 3
    ftype=fittype(fitEqn, 'Dependent', 'z','Independent','x', 'Options', opts); 
    [curve, gof] = fit(x,z, ftype);
else
    curve = fitnlm(x, z,fitFcn,opts.startpoint,'Weight', weights,...
        'CoefficientNames',{'R','x0','z0'});
    gof = [];
end
