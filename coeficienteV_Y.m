function VY=coeficienteV_Y(a,b,c,d,gamma,tol)
%--------------------------------------------------------------------------------
% Propósito :  Esta función calcula el valor del coeficiente que acomana a
%              la varianza del termino de error en el calculo de la
%              varianza de la desviacion del producto 
%                        V[Y] = VY*V[epsilon] 
%--------------------------------------------------------------------------------
% Inputs    : a     : 1x1 coeficiente de la solucion de k_{t+1}
%             b     : 1x1 coeficiente de la solucion de k_{t+1}
%             c     : 1x1 coeficiente de la solucion de Y_{t}
%             d     : 1x1 coeficiente de la solucion de Y_{t}
%             gamma : 1x1 coeficiente del AR(1) del proceso estocastico
%--------------------------------------------------------------------------------
% Output    : VY    : 1x1 Coeficientes VY en la formula de la V[Y]
%--------------------------------------------------------------------------------
VY=0;
for i=1:1000
    shortsum=0;
    if i==1
        shortsum=1;
    else
        for j=1:(i+1)
            shortsum=shortsum+a^(j-1)*gamma^(i+1-j);
        end
    end
    increment=c*b*shortsum+d*gamma^i;
    VY=VY+increment*increment;
    if abs(increment)<tol
        'tol achieved'
        increment
        i
        break
    end
end
VY=VY+d*d;
    