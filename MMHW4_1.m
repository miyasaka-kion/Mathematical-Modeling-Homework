clear
syms  t f(t) g(t)

f(t) = (0.3E1+exp(1).^((1/30).*t)).^(-2).*((-0.405E1)+(-0.845E1).*exp(1).^((1/15).*t)+exp(1).^((1/30).*t).*(0.253E2+(-0.8E0).*t));
g(t) = diff(f,t);
lt = 0;
for i = 1:200
   newt = double(lt - f(lt)/g(lt));
   if abs(double(lt-newt)) < 0.00001
       break;
   end
   lt = newt;
end

