function [Tabla61] = Errores_estandar_Hansen(P,Q,R,S,gamma)
%--------------------------------------------------------------------------------
% Propósito :  Calcula los errores estandar del Modelo
%--------------------------------------------------------------------------------
% Inputs    : P     : MxM Coeficientes de ley de movimiento lineal recursivo para las 
%                      variables estado
%             Q     : Mx1 Coeficientes de ley de movimiento lineal recursivo
%             R     : NxN Coeficientes de ley de movimiento lineal recursivo
%             S     : Nx1 Coeficientes de ley de movimiento lineal recursivo
%             gamma : 1x1 coeficiente del proceso estocastico AR(1)
%--------------------------------------------------------------------------------
% Output    : Tabla : Mx2 Errores estandar de las variables del modelo
%                         Porcentaje del error estandar de las variables 
%                         endogenas con respecto al PIB
%--------------------------------------------------------------------------------

a     = P(1,1);
b     = Q(1,1);
tol   = 0.0000001;

N   = size(R,1);
V_V = NaN(N,1);

for i=1:N
c        = R(i,1);
d        = S(i,1);
V_V(i,1) =coeficienteV_Y(a,b,c,d,gamma,tol);
end

Tabla_fila1 = sqrt(V_V);
Tabla_fila2 = NaN(N,1);
for i=1:N
    Tabla_fila2(i,1) = (Tabla_fila1(i,1)/Tabla_fila1(1,1))*100;
end 

Tabla61 = [Tabla_fila1 Tabla_fila2];

end 
