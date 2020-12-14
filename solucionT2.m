function [P,Q,R,S] = solucionT2(A,B,C,D,F,G,H,J,K,L,M,N)
%--------------------------------------------------------------------------------
% Propósito :  Esta función utiliza el resultaddo del teorema 2 de Uhlig
%              para encontrar las matrices P,Q, R y S de la solucion 
%                 x_t = P x_{t-1} + Q z_t
%                 y_t = R x_{t-1} + S z_t
%              del modelo escrito de la forma
%                 0 = A x_{t} + B x_{t-1} + C y_{t} + D z_{t}
%                 0 = E[F x_{t+1}+ G x_t + H x_{t-1}+J y_{t+1} + K y_{t} + L z_{t+1} + M z_{t}]
%                 z_{t+1} = N z_t + mu_{t+1}
%--------------------------------------------------------------------------------
% Inputs    : A  : NxM Matriz de coeficientes de las variables estado en t
%             B  : NxM Matriz de coeficientes de las variables estado en t-1
%             C  : NxN Matriz de coeficientes de las variables control en t
%             D  : Nxk Matriz de coeficientes de las variables estocasticas en t
%             F  : lxM Matriz de coeficientes de las variables estado en t+1
%             G  : lxM Matriz de coeficientes de las variables estado en t
%             H  : lxN Matriz de coeficientes de las variables estado en t-1
%             J  : lxM Matriz de coeficientes de las variables endogenas en t+1
%             K  : lxN Matriz de coeficientes de las variables endogenas en t
%             L  : lxk Matriz de coeficientes de las variables estocasticas en t+1
%             M  : lxk Matriz de coeficientes de las variables estocasticas en t
%             N  : kxk Matriz de coeficientes de la ley de movimiento de las 
%--------------------------------------------------------------------------------
% Output    : P  : MxM Coeficientes de ley de movimiento lineal recursivo para las 
%                      variables estado
%             Q  : Mx1 Coeficientes de ley de movimiento lineal recursivo
%             R  : NxN Coeficientes de ley de movimiento lineal recursivo
%             S  : Nx1 Coeficientes de ley de movimiento lineal recursivo
%--------------------------------------------------------------------------------

% encontrando P
auxA = F-J*inv(C)*A;
auxB = zeros(size(auxA,1),size(auxA,2));
auxC = - (J*inv(C)*B-G+K*inv(C)*A);
auxD = -K*inv(C)*B+H;

P    = fsolve(@(P)auxA*(P*P)+auxC*P + auxD, zeros(1,1));

% Encontrando Q, S, R

k        = size(N,1);   
R        = -inv(C)*(A*P+B);
Aux1     =  (J*inv(C)*D-L)*N+K*inv(C)*D-M;
vecAux1  = reshape(Aux1, [], 1);  % Vectorizacion por columnas
Aux2     = kron(N',F-J*inv(C)*A) + kron(ones(k,k),J*R+F*P+G-K*inv(C)*A);
vecQ     = inv(Aux2)*Aux1;
Q        = vecQ;                  % solo en este caso porque Q es 1x1 
S        = -inv(C)*(A*Q+D);  

end 
