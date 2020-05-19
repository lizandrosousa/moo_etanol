C=[]
t=[]

// função que descreve o modelo do fermentador
function g=fermenter(V)
    to=0
    Co=[5;0;V(2);0] 

function [f]=fun(t,C,V)
 
    Ks=23.587; Ki=103.502; Pmax=201.443; ms=0.272; kd=0.054;betamix=2.201; n=0.365; a=2.6; b=-0.282

    mumax=2.305*exp(-61.786/V(1))+-15.326*exp(-152.713/V(1))
    Yxs=-0.187*exp(-83.083/V(1))+0.05*exp(-3.288/V(1))
    Ypx=41.086*exp(-52.601/V(1))+43.330*exp(-52.015/V(1))
    Imix=150
       
      f(1)= (mumax*C(3)/(Ks+C(3)+(C(3))^2/Ki)*(1-C(4)/Pmax))*C(1)-kd*C(1)
      f(2)=kd*C(1)
      f(3)= -((betamix*(1-a*exp(b*C(1)))*C(3)*Imix^n)*((mumax*C(3)/(Ks+C(3)+(C(3))^2/Ki)*(1-C(4)/Pmax))/Yxs +ms))/((betamix*(1-a*exp(b*C(1)))*C(3)*Imix^n)+(((mumax*C(3)/(Ks+C(3)+(C(3))^2/Ki)*(1-C(4)/Pmax))/Yxs +ms)))*C(1)
      f(4)= Ypx*(mumax*C(3)/(Ks+C(3)+(C(3))^2/Ki)*(1-C(4)/Pmax))*C(1)

endfunction 
    t = V(3)
    dlist=list(fun,V)
    C=ode(Co,to,V(3),dlist)
    dev= 32- C(4)
//w1=0.3; w2=0.6; w3=0.1
    g(1)=-(C(4)/t)/9.3 +1*(max(0,dev))//  Produtividade
    g(2)=-(V(2)-C(3))/V(2) +1*(max(0,dev))//  Conversão
    g(3)=-(C(4)/(V(2)-C(3)))/0.511 +1*(max(0,dev))//  Eficiencia
    
    //x=-(w1*(C(4)/t)/9.3 +w2*(V(2)-C(3))/V(2) + w3/0.511*(C(4)/(V(2)-C(3))))+1000*max(0,dev)

          
endfunction


PopSize     = 10;
Proba_cross = 0.5;
Proba_mut   = 0.3;
NbGen       = 10;
NbCouples   = 110;
Log         = %T;
nb_disp     = 10; // Nb point to display from the optimal populationop
pressure    = 0.1;

ga_params = init_param();
ga_params = add_param(ga_params,'dimension',3);
ga_params = add_param(ga_params,'minbound',[26,70,0.8]);
ga_params = add_param(ga_params,'maxbound',[42,170,16]);

[pop_opt, fobj_pop_opt, pop_init, fobj_pop_init] = optim_nsga2(fermenter, PopSize,NbGen, Proba_mut, Proba_cross, Log, ga_params)

//pause
////
//plot3d(fobj_pop_opt(:,3),fobj_pop_opt(:,1),fobj_pop_opt(:,2),alpha=35,theta=45,leg='Yield@Productivity@Conversion')
