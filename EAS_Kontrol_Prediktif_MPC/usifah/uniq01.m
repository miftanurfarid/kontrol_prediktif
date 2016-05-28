function gam=uniq01(x,T)
nc=2;
R = [1.4311 0.92];
Q = [1.4320 1.4];
U(1,2) = 54.3372;
U(2,1) = 47.1062;
for i=1 : nc 
    sum=0;
    summ=0;
    for k=1:nc
        sum=sum+R(k)*x(k);
        summ=summ+Q(k)*x(k);
    end
    Fi(i) = R(i)*x(i)/sum;
    Tetha(i) = Q(i)*x(i)/summ;
    El(i)=5*(R(i)-Q(i))-(R(i)-1);
end
for i=1:nc
    for j=1:nc
        if i==j 
            tho(j,i)=1;
        end
        tho(j,i)=exp(-U(j,i)/1.987/T);
    end
end
for i=1:nc
    sigma=0;
    sigma1=0;
    sigma2=0;
    for j=1:nc
        sigma3=0;
        for k=1:nc
            sigma3=sigma3+Tetha(k)*tho(k,j);
        end
        sigma1=sigma1+Tetha(j)*tho(j,i);
        sigma2=sigma2+Tetha(j)*tho(i,j)/sigma3;
        sigma=sigma+x(j)*El(j);
    end
    s(i)=5*Q(i)*log(Tetha(i)/Fi(i));
    lngammaC(i)=log(Fi(i)/x(i))+s(i)+El(i)-Fi(i)*sigma/x(i);
    lngammaR(i)=Q(i)*(1-log(sigma1)-sigma2);
    lngamma(i)=lngammaC(i)+lngammaR(i);
    gam(i)=exp(lngamma(i));
end

        
    

    