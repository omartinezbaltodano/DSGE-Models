function [P,Q] = solucionT1(F,G,H,L,M,N)
%--------------------------------------------------------------------------------
% Propósito :  Esta función utiliza el resultado del teorema 1 de Uhlig
%              para encontrar las matrices P y Q de la solucion 
%                 x_t = P x_{t-1} + Q z_t
%              del modelo escrito de la forma
%                 0 = E[F x_{t+1}+ G x_t + H x_{t-1}+ L z_{t+1} + M z_{t}]
%                 z_{t+1} = N z_t + mu_{t+1}
%--------------------------------------------------------------------------------
% Inputs    : F  : NxN Matriz de coeficientes de las variables endogenas en t+1
%             G  : NxN Matriz de coeficientes de las variables endogenas en t
%             H  : NxN Matriz de coeficientes de las variables endogenas en t-1
%             L  : Nx1 Matriz de coeficientes de las variables estocasticas en t+1
%             M  : Nx1 Matriz de coeficientes de las variables estocasticas en t
%             N  : kxk Matriz de coeficientes de la ley de movimiento de las 
%                      variables estocasticas
%--------------------------------------------------------------------------------
% Output    : P  : NxN Coeficientes de ley de movimiento lineal recursivo
%             Q  : Nx1 Coeficientes de ley de movimiento lineal recursivo
%--------------------------------------------------------------------------------

% para entontrar P
auxB = zeros(size(G,1),size(G,2));
P    = fsolve(@(P)F*(P*P)+G*P + H, zeros(size(H,1),size(H,2)));

% para encontrar Q
k    = size(N,1);  
V    = kron(N',F)+kron(ones(k,k),F*P+G); 

A    = L*N+M;
vecA = reshape(A, [], 1); 
Q    = inv(V)*(-vecA);

end 
