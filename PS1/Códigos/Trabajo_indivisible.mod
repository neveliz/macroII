close all; 
clc;

var c w r y h k i lambda $\lambda$ (long_name='TFP')
    productividad ${\frac{y}{h}}$ (long_name='productividad');
    
varexo eps_a;

parameters beta $\beta$ (long_name='discount factor')
    delta $\delta$ (long_name='depreciation rate')
    theta $\theta$ (long_name='capital share')
    gamma $\gamma$ (long_name='AR coefficient TFP')
    A $A$ (long_name='labor disutility parameter')
    h_0 ${h_0}$ (long_name='full time workers in steady state')
    sigma_eps $\sigma_e$ (long_name='TFP shock volatility')
    B $B$ (long_name='composite labor disutility parameter')
    ;

//Calibración

beta = 0.99;
delta = 0.025;
theta = 0.36;
gamma = 0.95;
A = 2;
sigma_eps=0.007;
h_0=0.53;

model;
//1. Euler Equation
1/c = beta*((1/c(+1))*(r(+1) +(1-delta)));

//2. Labor FOC
(1-theta)*(y/h) = B*c;

//3. Resource constraint
c = y +(1-delta)*k(-1) - k;

//4. LOM capital
k = (1-delta)*k(-1) + i;

//5. Production function
y = lambda*k(-1)^(theta)*h^(1-theta);

//6. Real wage
r = theta*(y/k(-1));

//7. Real interest rate
w = (1-theta)*(y/h);

//8. LOM TFP
log(lambda)=gamma*log(lambda(-1))+eps_a;

//9. productividad
productividad= y/h;

end;

steady_state_model;
B=-A*(log(1-h_0))/h_0;
lambda = 1;
h = (1-theta)*(1/beta -(1-delta))/(B*(1/beta -(1-delta)-theta*delta));
k = h*((1/beta -(1-delta))/(theta*lambda))^(1/(theta-1));
i = delta*k;
y = lambda*k^(theta)*h^(1-theta);
c = y-delta*k;
r =  1/beta - (1-delta);
w = (1-theta)*(y/h);
productividad = y/h;
end;

steady;

shocks;
var eps_a; stderr sigma_eps;
end;

check;
steady;
stoch_simul(order=1,irf=20,loglinear,hp_filter=1600, periods=200, simul_replic = 10000);

% Recolecto todos los resultados del Dynare
% Recopilamos los resultados
results = get_simul_replications(M_,options_);

% Determinamos el N° de Simulaciones a descartar
N = 0;

% Aplicando el Filtro de Hodrick - Prescott
simulated_series_filtered = NaN(size(results));
for ii=1:options_.simul_replic
    [trend, cycle]=sample_hp_filter(results(:,:,ii)',1600);
    simulated_series_filtered(:,:,ii)=cycle';
end

% Sacamos la Posición de las Variables
c_pos=strmatch('c',M_.endo_names,'exact');
w_pos = strmatch('w',M_.endo_names,'exact');
r_pos=strmatch('r',M_.endo_names,'exact');
y_pos = strmatch('y',M_.endo_names,'exact');
h_pos=strmatch('h',M_.endo_names,'exact');
k_pos=strmatch('k',M_.endo_names,'exact');
i_pos = strmatch('i',M_.endo_names,'exact');
z_pos = strmatch('z',M_.endo_names,'exact');
g_pos = strmatch('g',M_.endo_names,'exact');
prod_pos=strmatch('productividad',M_.endo_names,'exact');

%% Estadísticas Descriptivas de las Variables
%  Desviación y Correlación Promedio - 10.000 Simulaciones

% Definimos las posiciones de las variables
var_positions = [c_pos; w_pos; r_pos; y_pos; h_pos; k_pos; i_pos; z_pos; g_pos; prod_pos];

% Nombres de las Variables
var_names = M_.endo_names_long(var_positions,:);

% Calculamos la Desviación Estándar
std_mat = std(simulated_series_filtered(var_positions,:,N+1:end),0,2)*100;

% Almacenamos todos los resultados
corr_mat = zeros(6,options_.simul_replic - N);
stats_model = zeros(6,options_.simul_replic - N);

% Calculo Correlaciones
for ii=1:options_.simul_replic - N
   corr_mat(1,ii)=corr(results(y_pos,:,ii)',results(y_pos,:,ii)');
   corr_mat(2,ii)=corr(results(y_pos,:,ii)',results(c_pos,:,ii)');
   corr_mat(3,ii)=corr(results(y_pos,:,ii)',results(i_pos,:,ii)');
   corr_mat(4,ii)=corr(results(y_pos,:,ii)',results(k_pos,:,ii)');
   corr_mat(5,ii)=corr(results(y_pos,:,ii)',results(h_pos,:,ii)');
   corr_mat(6,ii)=corr(results(y_pos,:,ii)',results(prod_pos,:,ii)');
end

% Calculo Desviaciones Relativas
for jj = 1:options_.simul_replic - N
    stats_model(1,jj) = std_mat(c_pos,:,jj)/std_mat(y_pos,:,jj);
    stats_model(2,jj) = std_mat(i_pos,:,jj)/std_mat(y_pos,:,jj);
    stats_model(3,jj) = std_mat(h_pos,:,jj)/std_mat(y_pos,:,jj);
    stats_model(4,jj) = std_mat(w_pos,:,jj)/std_mat(y_pos,:,jj);
    stats_model(5,jj) = std_mat(h_pos,:,jj)/std_mat(w_pos,:,jj);
    stats_model(6,jj) = corr(results(c_pos,:,jj)',results(y_pos,:,jj)');
end

% Tabla 1: Estadísticas de las Variables
fprintf('\n');
fprintf('----------------------------------------------------- \n');
fprintf('%-40s \n',"ESTADÍSTICAS DE LAS VARIABLES*");
fprintf('----------------------------------------------------- \n');
fprintf('%-20s \t %11s \t %11s \n','','std(x)','corr(y,x)')
for ii=1:size(corr_mat,1)
    fprintf('%-20s \t %3.2f (%3.2f) \t %3.2f (%3.2f) \n',var_names{ii,:},mean(std_mat(ii,:,:),3),std(std_mat(ii,:,:),0,3),mean(corr_mat(ii,:),2),std(corr_mat(ii,:),0,2))
end
fprintf('----------------------------------------------------- \n');

% Tabla 2: Desviaciones Relativas Variables
fprintf('\n');
fprintf('------------------------- \n');
fprintf('%-40s \n','DESVIACIONES RELATIVAS*');
fprintf('------------------------- \n');
fprintf('std(c)/std(y) \t %3.2f \n', mean(stats_model(1,:),'all'));
fprintf('std(I)/std(y) \t %3.2f \n', mean(stats_model(2,:),'all'));
fprintf('std(h)/std(y) \t %3.2f \n', mean(stats_model(3,:),'all'));
fprintf('std(w)/std(y) \t %3.2f \n', mean(stats_model(4,:),'all'));
fprintf('std(h)/std(w) \t %3.2f \n', mean(stats_model(5,:),'all'));
fprintf('corr(h,w)     \t %3.2f \n', mean(stats_model(6,:),'all'));
fprintf('------------------------- \n');

fprintf('\n');
fprintf('%-40s \n','* Nota:');
fprintf('----------------------------------------------------- \n');
fprintf('- N° de Simulaciones consideradas: %d de %d \n',options_.simul_replic-N,options_.simul_replic);
fprintf('- N° de Periodos: %d\n', options_.periods);
fprintf('----------------------------------------------------- \n');

%%% Histograma de las Estadísticas
%% Según Paper Hansen and Wright (1992)

% Crear una figura
fig_histogram = figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);

% Crear subplots para cada histograma
subplot(2, 4, 1);
histogram(std_mat(4,:,:), 'FaceColor', 'b', 'FaceAlpha', 0.3);
hold on;
xline(mean(std_mat(4,:,:)), 'LineWidth', 2, 'Color', 'r', 'LineStyle', '-');
title('Desviación del Producto');
xlabel('std(y)');
ylabel('Frecuencia');

subplot(2, 4, 2);
histogram(stats_model(1,:), 'FaceColor', 'b', 'FaceAlpha', 0.3);
hold on;
xline(mean(stats_model(1,:)), 'LineWidth', 2, 'Color', 'r', 'LineStyle', '-');
title('Desviación Relativa Consumo/Producto');
xlabel('std(c)/std(y)');
ylabel('Frecuencia');

subplot(2, 4, 3);
histogram(stats_model(2,:), 'FaceColor', 'b', 'FaceAlpha', 0.3);
hold on;
xline(mean(stats_model(2,:)), 'LineWidth', 2, 'Color', 'r', 'LineStyle', '-');
title('Desviación Relativa Inversión/Producto');
xlabel('std(I)/std(y)');
ylabel('Frecuencia');

subplot(2, 4, 4);
histogram(stats_model(3,:), 'FaceColor', 'b', 'FaceAlpha', 0.3);
hold on;
xline(mean(stats_model(3,:)), 'LineWidth', 2, 'Color', 'r', 'LineStyle', '-');
title('Desviación Relativa Horas/Producto');
xlabel('std(h)/std(y)');
ylabel('Frecuencia');

subplot(2, 4, 5);
histogram(stats_model(4,:), 'FaceColor', 'b', 'FaceAlpha', 0.3);
hold on;
xline(mean(stats_model(4,:)), 'LineWidth', 2, 'Color', 'r', 'LineStyle', '-');
title('Desviación Relativa Salario/Producto');
xlabel('std(w)/std(y)');
ylabel('Frecuencia');

subplot(2, 4, 6);
histogram(stats_model(5,:), 'FaceColor', 'b', 'FaceAlpha', 0.3);
hold on;
xline(mean(stats_model(5,:)), 'LineWidth', 2, 'Color', 'r', 'LineStyle', '-');
title('Desviación Relativa Horas/Salario');
xlabel('std(h)/std(w)');
ylabel('Frecuencia');

subplot(2, 4, 7);
histogram(stats_model(6,:), 'FaceColor', 'b', 'FaceAlpha', 0.3);
hold on;
xline(mean(stats_model(6,:)), 'LineWidth', 2, 'Color', 'r', 'LineStyle', '-');
title('Correlación Horas/Salario');
xlabel('corr(h,w)');
ylabel('Frecuencia');

% Añadir un título general a la figura
sgtitle('Estadísticas: Modelo - Trabajo Indivisible');

% Ajustar el tamaño de la figura y el papel
set(fig_histogram, 'PaperPositionMode', 'auto');
set(fig_histogram, 'PaperOrientation', 'landscape');
set(fig_histogram, 'PaperUnits', 'normalized');
set(fig_histogram, 'PaperPosition', [0 0 1 1]);

% Guardar la figura en formato PNG
exportgraphics(fig_histogram, 'stats_Trabajo_indivisible.png', 'Resolution', 300);