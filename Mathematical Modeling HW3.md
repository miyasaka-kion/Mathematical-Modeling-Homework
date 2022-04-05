# Mathematical Modeling HW 3

Miyasaka Kion

## Newton Iteration

The objective function is:
$$
f(x)=(0.65-0.01 x)\left(200 e^{c x}\right)-0.45 x
$$
and we want to find the root of $f'(x)$, where
$$
f'(x) = F(\text{x})=(-2 c x+130 c-2) \exp (c x)-0.45
$$

```mathematica
Clear["Global`*"]
L := {0.0230, 0.0235, 0.0240, 0.0245, 0.0250, 0.0255, 0.0260, 0.0265, 
  0.0270}
(*c= 0.028*)

f[x_] := (130 c \[Minus] 2 \[Minus] 2*c*x)*Exp[c*x] \[Minus] 0.45;

ans := Table[FixedPoint[(# - f[#]/f'[#]) &, 1, 10], {c, L}]
ttl = Transpose[{L, ans}] // TableForm
```

And the output is 
$$
c\qquad\quad x\\
\left(
\begin{array}{cc}
 0.023 & 14.5159 \\
 0.0235 & 15.8497 \\
 0.024 & 17.1166 \\
 0.0245 & 18.3213 \\
 0.025 & 19.4682 \\
 0.0255 & 20.5611 \\
 0.026 & 21.6037 \\
 0.0265 & 22.5992 \\
 0.027 & 23.5507 \\
\end{array}
\right)
$$
, where column 1 reprsents $c$ and column 2 represents $x$.

(才发现好像要求要用MATLAB)

So we implement the algorithm in matlab.

```matlab
syms  x f(x) g(x) 
L = [0.0230, 0.0235, 0.0240, 0.0245, 0.0250, 0.0255, 0.0260, 0.0265, 0.0270];
ans = zeros(1,length(L));
for j = 1:length(L)
    c = L(j);
    f(x) = (130*c - 2 - 2*c*x)*exp(c*x) - 0.45;
    for i = 1:20
        g(x) = diff(f,x);
       newx = double(lx - f(lx)/g(lx));
       if abs(lx-newx) < 0.00001
           break
       end
       lx = newx;
    end
    ans(j) = lx;
end

ans

res = [L; ans].'
```

and the output is

```
ans =

   14.5159   15.8497   17.1166   18.3213   19.4682   20.5611   21.6037   22.5992   23.5507


res =

    0.0230   14.5159
    0.0235   15.8497
    0.0240   17.1166
    0.0245   18.3213
    0.0250   19.4682
    0.0255   20.5611
    0.0260   21.6037
    0.0265   22.5992
    0.0270   23.5507
```

which comes the same result.

## Random Search

The objective function is defined as
$$
\begin{aligned}
&z=3.2+1.7\left[6 \sqrt{(x-1)^{2}+(y-5)^{2}}^{ 0.91}+8 \sqrt{(x-3)^{2}+(y-5)^{2}}^{0.91}\right.\\
&+8{\sqrt{(x-5)^{2}+(y-5)^{2}}}^{0.91}+21{\sqrt{(x-1)^{2}+(y-3)^{2}}}^{0.91}\\
&+6{\sqrt{(x-3)^{2}+(y-3)^{2}}}^{0.91}+3{\sqrt{(x-5)^{2}+(y-3)^{2}}}^{0.91}\\
&+18{\sqrt{(x-1)^{2}+(y-1)^{2}}}^{0.91}+8{\sqrt{(x-3)^{2}+(y-1)^{2}}}^{0.91}\\
&\left.+6 \sqrt{(x-5)^{2}+(y-1)^{2}}^{0.91}\right] / 84
\end{aligned}
$$
We can first use mathematica to earn a relatively accurate solution.

```mathematica
f[x_, y_, a_, b_] := Sqrt[((x - a)^2 + (y - b)^2)]^0.91
z[x_, y_] := 
 3.2 + 1.7*(6*f[x, y, 1, 5] + 8*f[x, y, 3, 5] + 8*f[x, y, 5, 5] + 
      21*f[x, y, 1, 3] + 6*f[x, y, 3, 3] + 3*f[x, y, 5, 3] + 
      18*f[x, y, 1, 1] + 8*f[x, y, 3, 1] + 6*f[x, y, 5, 1])/84
Minimize[z[x, y], {x, y}]
```

and the output is

```mathematica
{6.46298, {x -> 1.61649, y -> 2.76591}}
```

The code below is random search implemented in MATLAB.

```matlab
clear 
syms x y f(x,y,a,b) z(x,y)
f(x,y,a,b) = sqrt((x-a)^2+(y-b)^2)^0.91;
z(x,y) = 3.2+1.7*( 6*f(x,y,1,5) + 8*f(x,y,3,5) + 8*f(x,y,5,5) + 21*f(x,y,1,3) + 6*f(x,y,3,3) + 3*f(x,y,5,3) +18*f(x,y,1,1) + 8*f(x,y,3,1) + 6*f(x,y,5,1))/84;
% z(x,y) = simplify(z(x,y));
range_y = [0 6];
range_x = range_y;
INF = 0x7fffffff;
N = 0;
sz = [1 2];
zmin = INF;
minx = -INF;
miny = -INF;
std_ans = 6.46298;
A = [];

while N <= 1000
    N = N + 50;
    
    tic
    for i = 1:N
       cx = (range_x(2) - range_x(1)) * rand + range_x(1);
       cy = (range_y(2) - range_y(1)) * rand + range_y(1);
       [cx,cy];
       cz = double(z(cx,cy));
       if cz < zmin 
           zmin = cz;
           xmin = cx;
           ymin = cy;
       end      
    end
    time_cost = toc;
    bs = abs(std_ans - zmin)/std_ans;
    kkk = ['N = ', num2str(N) ,', minimum value it found is ', num2str(zmin),', the relative bias is ', num2str(bs),' and costs ', num2str(time_cost) ,'s.'];
    disp(kkk); 
    A = [A, [N; zmin; bs; time_cost]];
    
end


plot(A(1,:),A(2,:));
exportgraphics(gcf,'MMHW321.png','Resolution',300)
plot(A(1,:),A(3,:));
exportgraphics(gcf,'MMHW322.png','Resolution',300)
plot(A(1,:),A(4,:));
exportgraphics(gcf,'MMHW323.png','Resolution',300)


```

and the out put is 

```matlab
>> MMHW3_2
N = 1, minimum value it found is 8.5979, the relative bias is 0.33033 and costs 0.011482s.
N = 2, minimum value it found is 8.5979, the relative bias is 0.33033 and costs 0.026942s.
N = 3, minimum value it found is 7.1652, the relative bias is 0.10866 and costs 0.037644s.
N = 4, minimum value it found is 7.157, the relative bias is 0.10738 and costs 0.044422s.
N = 5, minimum value it found is 6.6327, the relative bias is 0.026267 and costs 0.064061s.
N = 6, minimum value it found is 6.6097, the relative bias is 0.022699 and costs 0.085334s.
N = 7, minimum value it found is 6.6097, the relative bias is 0.022699 and costs 0.072112s.
N = 8, minimum value it found is 6.6097, the relative bias is 0.022699 and costs 0.092397s.
N = 9, minimum value it found is 6.6097, the relative bias is 0.022699 and costs 0.093514s.
N = 10, minimum value it found is 6.5823, the relative bias is 0.018461 and costs 0.11738s.
N = 11, minimum value it found is 6.4734, the relative bias is 0.0016143 and costs 0.11806s.
>> MMHW3_2
N = 50, minimum value it found is 6.489, the relative bias is 0.0040194 and costs 0.55517s.
N = 100, minimum value it found is 6.4783, the relative bias is 0.0023668 and costs 1.1277s.
N = 150, minimum value it found is 6.4717, the relative bias is 0.0013508 and costs 1.7138s.
N = 200, minimum value it found is 6.4669, the relative bias is 0.00061321 and costs 2.1679s.
N = 250, minimum value it found is 6.4669, the relative bias is 0.00061321 and costs 2.7764s.
N = 300, minimum value it found is 6.4669, the relative bias is 0.00061321 and costs 3.1765s.
N = 350, minimum value it found is 6.4669, the relative bias is 0.00061321 and costs 3.5677s.
N = 400, minimum value it found is 6.4669, the relative bias is 0.00061321 and costs 4.1528s.
N = 450, minimum value it found is 6.4637, the relative bias is 0.00011295 and costs 4.6667s.
N = 500, minimum value it found is 6.4633, the relative bias is 4.4306e-05 and costs 4.9786s.
N = 550, minimum value it found is 6.4633, the relative bias is 4.4306e-05 and costs 5.588s.
N = 600, minimum value it found is 6.4633, the relative bias is 4.4306e-05 and costs 6.8878s.
N = 650, minimum value it found is 6.4633, the relative bias is 4.4306e-05 and costs 7.1201s.
N = 700, minimum value it found is 6.4633, the relative bias is 4.4306e-05 and costs 8.1593s.
N = 750, minimum value it found is 6.4633, the relative bias is 4.4306e-05 and costs 8.7114s.
N = 800, minimum value it found is 6.4633, the relative bias is 4.4306e-05 and costs 8.893s.
N = 850, minimum value it found is 6.4633, the relative bias is 4.4306e-05 and costs 9.8643s.
N = 900, minimum value it found is 6.4633, the relative bias is 4.4306e-05 and costs 9.5887s.
N = 950, minimum value it found is 6.4633, the relative bias is 4.4306e-05 and costs 9.6306s.
N = 1000, minimum value it found is 6.4633, the relative bias is 4.4306e-05 and costs 10.2337s.
N = 1050, minimum value it found is 6.4633, the relative bias is 4.4306e-05 and costs 11.1621s.
```



### Visualize

the minimum value w.r.t. N:

![MMHW321](https://tva1.sinaimg.cn/large/e6c9d24egy1h0qktzreiaj213h0u00tt.jpg)

The relative bias w.r.t. N:

![MMHW322](https://tva1.sinaimg.cn/large/e6c9d24egy1h0qkv0wrsyj211l0u0ab3.jpg)

Time cost w.r.t. N:

![MMHW323](https://tva1.sinaimg.cn/large/e6c9d24egy1h0qkvlhq9nj212a0u0my8.jpg)

It can be seen that the larger N is, the more accurate the result is (but the accuracy remains unchanged after a certain level), and the operation time as a whole has been increasing within a certain range.



## Multivariable Newton’s method

Objective funcion:
$$
\begin{aligned}
z=& x\left(10+31 x^{-0.5}+1.3 y^{-0.2}\right)-18 x \\
&+y\left(5+15 y^{-0.4}+0.8 x^{-0.08}\right)-10 y
\end{aligned}
$$
To find global optimum, we need $\nabla z  = \bold 0$, i.e.
$$
\left\{-\frac{0.064 y}{x^{1.08}}+\frac{15.5}{x^{0.5}}+\frac{1.3}{y^{0.2}}-8,\frac{0.8}{x^{0.08}}-\frac{0.26 x}{y^{1.2}}-\frac{6.}{y^{0.4}}+\frac{15}{y^{0.4}}-5\right\} = \bold 0
$$
Let $\bold F : \mathbb R ^2\rightarrow \mathbb R ^2 $ be a vector field.
$$
\bold F = \left\{-\frac{0.064 y}{x^{1.08}}+\frac{15.5}{x^{0.5}}+\frac{1.3}{y^{0.2}}-8,\frac{0.8}{x^{0.08}}-\frac{0.26 x}{y^{1.2}}-\frac{6.}{y^{0.4}}+\frac{15}{y^{0.4}}-5\right\}
$$
then we can apply multivariable Newton's method to solve $\bold F = \bold 0$.



we can use some trick to avoid inputing the tedious function.

```mathematica
Clear[z, x, y]
z [x_, y_] := 
 x*(10 + 31 x^(-0.5) + 1.3 y^(-0.2)) - 18 x + 
  y (5 + 15 y^(-0.4) + 0.8 x^(-0.08)) - 10 y
Grad[z[x, y], {x, y}] // ToMatlab
```

output:

```mathematica
[(-8)+0.155E2.*x.^(-0.5E0)+0.13E1.*y.^(-0.2E0)+(-0.64E-1).*x.^( ...
  \
-0.108E1).*y,(-5)+0.8E0.*x.^(-0.8E-1)+(-0.26E0).*x.*y.^(-0.12E1)+ ...
\
  15.*y.^(-0.4E0)+(-0.6E1).*y.^(-0.4E0)];
```



### Multivariable Newton’s method, MATLB:

the iteratioin step can be written as:
$$
\bold x_{\text{new}} = \bold x_{\text{pre}} - \left( \operatorname{Jacobian}(\bold F(\bold x_{\text{pre}})) \right)^{-1} \bold F(\bold x_{\text{pre}})
$$


```matlab
clear 
syms z(x,y) x y F(x,y) CA(x,y) ICA(x,y)
z(x,y) = x*(10 + 31*x^(-0.5) + 1.3*y^(-0.2)) - 18*x + y *(5 + 15* y^(-0.4) + 0.8*x^(-0.08)) - 10 *y;
F(x,y) = [(-8)+0.155E2.*x.^(-0.5E0)+0.13E1.*y.^(-0.2E0)+(-0.64E-1).*x.^( ...
  -0.108E1).*y, (-5)+0.8E0.*x.^(-0.8E-1)+(-0.26E0).*x.*y.^(-0.12E1)+ ...
  15.*y.^(-0.4E0)+(-0.6E1).*y.^(-0.4E0)];
CA(x,y) = jacobian(F(x,y),[x y]);
ICA(x,y) = simplify((CA(x,y))^-1);

initial_val = [1;1];
cx = initial_val;
epsi = 0.00001;
INF = 0x7fffffff;
cnt = 0;
while 1==1
    newx = double(cx - double(ICA(cx(1),cx(2)) ) * double(F(cx(1),cx(2))).' );
    cnt = cnt + 1;
    if norm(newx - cx) <= epsi
        break;
    end
    cx = newx;
end



rng = [1,8];
for i = rng(1):rng(2)
    for j = rng(1):rng(2)
        cnt = 0;
        cx = [i; j];
        while 1==1
            newx = double(cx - double(ICA(cx(1),cx(2)) ) * double(F(cx(1),cx(2))).' );
            cnt = cnt + 1;
            if norm(newx - cx) <= epsi
                break;
            end
            cx = newx;
        end
        ANS(i-rng(1)+1,j-rng(1)+1) = cnt;
    end
end
newx
double(z(newx(1),newx(2)))
ANS

```

Out put

```matlab

newx =

    4.6896
    5.8520


ans =

   52.0727


ANS =

     7     7     7     7     7     7     7     7
     7     6     6     6     6     6     6     6
     7     6     5     5     5     5     5     5
     6     6     5     5     4     4     4     5
     6     6     5     5     4     4     4     5
     6     6     5     5     5     5     5     5
     6     6     5     5     5     5     5     5
     7     6     6     6     6     6     6     6
```



Hence the optimum value is $x^* = (4.6896,5.8520)$ and the optimum value is $52.0727$.

The meaning of `ANS` matrix is, if the initiative value is `cx = (i,j)`, the number of iterations required is `ANS(i,j)` (in the precision of $\epsilon = 10 ^{-5}$), i.e.
$$
\text{ANS}\\
\begin{array}{lllllllll}
    & x_1 &x_2 &x_3 &x_4 &x_5 &x_6 &x_7 &x_8 \\
y_1 & 7 & 7 & 7 & 7 & 7 & 7 & 7 & 7 \\
y_2 & 7 & 6 & 6 & 6 & 6 & 6 & 6 & 6 \\
y_3 & 7 & 6 & 5 & 5 & 5 & 5 & 5 & 5 \\
y_4 & 6 & 6 & 5 & 5 & 4 & 4 & 4 & 5 \\
y_5 & 6 & 6 & 5 & 5 & 4 & 4 & 4 & 5 \\
y_6 & 6 & 6 & 5 & 5 & 5 & 5 & 5 & 5 \\
y_7 & 6 & 6 & 5 & 5 & 5 & 5 & 5 & 5 \\
y_8 & 7 & 6 & 6 & 6 & 6 & 6 & 6 & 6
\end{array}
$$
