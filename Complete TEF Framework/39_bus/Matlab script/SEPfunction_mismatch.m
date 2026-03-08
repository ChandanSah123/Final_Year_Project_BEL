function [c,ceq]=SEPfunction(x,Pm,E,C,D,H,Y)
c=[];
Pcoi=0;
for i=1:numel(Pm)
    for j=i+1:numel(Pm)
        Pcoi=Pcoi+D(i,j)*cos(x(i)-x(j));
    end
end
Pcoi=Pcoi*2;
g=10;
ceq=zeros(numel(Pm),1);
for i=1:numel(Pm)
    for k=1:g
        if x(k)>pi
            x(k)=-(2*pi-x(k));
        end
        if x(k)<-pi
            x(k)=(2*pi+x(k));
        end
    end
    for j=1:numel(Pm)        
        if i~=j
            ceq(i)=ceq(i)-(C(i,j)*sin(x(i)-x(j))+D(i,j)*cos(x(i)-x(j)));
        else
            ceq(i)=ceq(i)+Pm(i)-E(i)^2*real(Y(i,i));
        end
    end
    ceq(i)=ceq(i)-(H(i)/sum(H))*sum(Pm-E.^2.*real(diag(Y))')+(H(i)/sum(H))*Pcoi;
end
 
end