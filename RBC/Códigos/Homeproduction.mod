var k km kh l ym yh cm ch hm hh rm rh wm wh zm zh; 
varexo eps_m eps_h;

parameters beta delta theta rho A a eta e;
beta = 0.99;
delta = 0.025;
theta = 0.36;
rho = 0.95;
A = 1;
a = 0.3;
eta = 0.08;
e = 0.8;

model;
//1 Euler of Market Consumption
(a*cm^(e-1))/(a*cm^(e)+(1-a)*ch^(e)) = beta*((a*cm(+1)^(e-1))/(a*cm(+1)^(e)+(1-a)*ch(+1)^(e)))*(rm(+1) + (1-delta));

//2 Euler of Home Consumption
(a*cm^(e-1))/(a*cm^(e)+(1-a)*ch^(e)) = beta*((a*cm(+1)^(e-1))/(a*cm(+1)^(e)+(1-a)*ch(+1)^(e)))*(1-delta) + beta*((1-a)*ch(+1)^(e-1)/(a*cm(+1)^(e)+(1-a)*ch(+1)^(e)))*rh(+1);

//3 Market Labor Supply
A/l = (a*cm^(e-1))*wm/(a*cm^(e)+(1-a)*ch^(e));

//4 Home Labor Supply
A/l = ((1-a)*ch^(e-1))*wh/(a*cm^(e)+(1-a)*ch^(e));

//5 Market Resource Constraint
cm + k - (1-delta)*k(-1) = ym;

//6 Home Resource Constraint
ch = yh;

//7 Labor Constraint
l = 1 - hh - hm;

//8 Capital Allocation Constraint
k = km + kh;

//9 Market Production Function
ym = zm*km^(theta)*hm^(1-theta);

//10 Home Production Function
yh = zh*kh^(eta)*hh^(1-eta);

//11 Market Capital Return
rm = theta*(ym/km);

//12 Home Capital Return
rh = eta*(yh/kh);

//13 Market Real Wage
wm = (1-theta)*(ym/hm);

//14 Home Real Wage
wh = (1-eta)*(yh/hh);

//15 Market Shock Process
log(zm) = rho*log(zm(-1)) + eps_m;

//16 Home Shock Process
log(zh) = rho*log(zh(-1)) + eps_h;
end;

initval;
k = 1; km = 0.5; kh = 0.5; l = 0.33; ym = 1; yh = 1; cm = 0.8; ch = 0.8;
hm = 0.33; hh = 0.33; rm = 0.1; rh = 0.1; wm = 1; wh = 1; zm = 1; zh = 1;
end;

steady;

shocks;
var eps_m = 0.007^2;
var eps_h = 0.007^2;
end;

check;

stoch_simul(order=1,irf=20,loglinear,hp_filter=1600, periods=179) cm k km kh hh;
