close all

% Frequenza di campionamento
Ws = 20;
T = 2*pi / Ws;
% Numero di dati raccolti
n_samples = 20;
% Funzione da discretizzare
func = @(t) sin(4*t);

t = -T*n_samples:0.01:T*n_samples;
data = zeros(size(t));

figure();
hold on
grid on
for ii = -n_samples:n_samples
    % sinc(t) = sinc(pi*t) in matlab
    data = data + (func(ii*T)*sinc(pi*(t - ii*T)/(T*pi)));   
    %fplot(@(t) func(ii*T)*sinc(pi*(t - ii*T)/T),"Color",[0.6,0.6,0.6])
end
s = scatter((-n_samples:n_samples)*T,func((-n_samples:n_samples)*T),"MarkerEdgeColor",[0,0,1]);
p = plot(t,data,"Color",[1,0,0],"LineWidth",3);
f = fplot(func,"Color",[0,0,1]);
legend([s,p,f],["samplings","ricostruzione","funzione originale"]);

