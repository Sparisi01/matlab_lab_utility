%% Funzione che minimizza il chi-quadro e calcola gli errori sui parametri per un modello dato
function [par, errpar, yfin, chi2norm, dof, p_value] = regressione_non_lin(x, y, dy, model, parameters, ub, lb)
arguments
    x(1,:) % x data
    y (1,:) % y data 
    dy (1,:) % y uncertainties
    model (1,1) % Model function (parameters, x) => y
    parameters (1,:) % Initial model parameters
    ub (1, :) {mustBeNonempty} = ones(size(parameters)) * inf; % Upper bound parameters
    lb (1, :) {mustBeNonempty} = ones(size(parameters)) * -inf; % Lower bound parameters
end

% Check consistent array dimension ----------------------------------------

if size(y) ~= size(x)
    error("x and y must be of the same size");
end

if (length(ub) ~= length(parameters)) || (length(lb) ~= length(parameters))
    error("ub, lb and parameters must be of the same size");
end

% -------------------------------------------------------------------------

options=optimset('lsqnonlin');

% Definizione funzione scarti a partire dal modello
scarti = @(par,xd,yd,ed) (model(par,xd)-yd)./ed;

[par,resnorm,~,~,~,~,jacobian] = lsqnonlin(scarti,parameters,lb,ub,options,x,y,dy);

%Covariance Matrix
covar = inv(jacobian' * jacobian);
%Variance
var = diag(covar);
%sigma (+- sigma = intervallo di confidenza al 68.3 % considerando un solo par per volta)
sigma=sqrt(var);
sigmaf=full(sigma);
dof = (length(x)-length(par));
chi2norm=resnorm/dof;

errpar=sigmaf*sqrt(chi2norm);
% Data with the new parameters
yfin = model(par,x);
p_value = chi2cdf(resnorm,dof);
end
