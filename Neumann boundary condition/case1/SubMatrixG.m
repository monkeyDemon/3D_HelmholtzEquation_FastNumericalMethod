function  [BMRE2, BMIM] = SubMatrixG(n, wavek, a, b)

h = (b - a) / (n + 1) ;
p = a : h : b ;
col = p(2 : n + 1) ;
BMRE=zeros(1,n);
BMIM=zeros(1,n);
%      2.1 Compute BMRE(1,j), except for the singular part 
%          i.e., -2/(pi*x)+x*log(x)/pi.
for j = 1 : n 
    m1 = 1 ;
    subh = 2 * h ;
	sum2 = 0 ;
    m1 = m1 * 2 ;
	subh = subh * 0.5 ;
    sum1 = sum2 ;
    sum2 = 0 ;
	for k = 1 : 2 : m1-1
        xx = p(j) + 2 * k * h / m1 ;
        fi = phi(j, n, p, xx) ;
    	fenmu = abs(xx - p(2)) ;
        if fenmu > 1e-12
            xxx = wavek * fenmu ;
            if j == 1
                bj11 = Bessel_mdf1Y1(xxx) ;
            else
                bj11 = Bessel_mdf2Y1(xxx) ;
            end
            wjmtemp = bj11 / fenmu ;
        else
            if j == 1
                wjmtemp = wavek * (-0.196057111305026) ;
            else
                wjmtemp = 0 ;
            end
        end
        sum2 = sum2 + wjmtemp * fi ;
    end
    sum2 = 0.5 * sum1 + subh * sum2 ;
    while abs(sum1-sum2) > h^3
        m1 = m1 * 2 ;
        subh = subh * 0.5 ;
        sum1 = sum2 ;
        sum2 = 0 ;
    	for k = 1 : 2 : m1-1
            xx = p(j) + 2 * k * h / m1 ;
            fi = phi(j, n, p, xx) ;
        	fenmu = abs(xx - p(2)) ;
            if fenmu > 1e-12
                xxx = wavek * fenmu ;
                if j == 1
                    bj11 = Bessel_mdf1Y1(xxx) ;
                else
                    bj11 = Bessel_mdf2Y1(xxx) ;
                end
                wjmtemp = bj11 / fenmu ;
            else
                if j == 1
                    wjmtemp = wavek * (-0.196057111305026) ;
                else
                    wjmtemp = 0 ;
                end
            end
            sum2 = sum2 + wjmtemp * fi ;
        end
        sum2 = 0.5 * sum1 + subh * sum2 ;
    end
    BMRE(1,j) = sum2 * 0.5 * wavek ;
end
BMRE(1,1) = BMRE(1,1) + wavek^2 * h * (log(wavek * h)-1.5)/(2 * pi) ;

%       2.3 Compute BMIM(1,j)	
for	j = 1 : n
	m1 = 1 ;
    subh = 2 * h ;
	sum2 = 0 ;
    m1 = m1 * 2 ;
	subh = subh * 0.5 ;
    sum1 = sum2 ;
    sum2 = 0 ;
    for k = 1 : 2 : m1-1
        xx = p(j) + 2 * k * h / m1 ;
        fi = phi(j, n, p, xx) ;
        fenmu = abs(xx - p(2)) ;
        if fenmu > 1e-12
            xxx = wavek * fenmu ;
            bj11 = Bessel_J1(xxx) ;
            wjmtemp = bj11 / fenmu ;
        else
            wjmtemp = 0.5 * wavek ;
        end
        sum2 = sum2 + wjmtemp * fi ;
    end
    sum2 = 0.5 * sum1 + subh * sum2 ;
	while  abs(sum1 - sum2) > h^3
            m1 = m1 * 2 ;
        	subh = subh * 0.5 ;
            sum1 = sum2 ;
            sum2 = 0 ;
            for k = 1 : 2 : m1-1
                xx = p(j) + 2 * k * h / m1 ;
                fi = phi(j, n, p, xx) ;
                fenmu = abs(xx - p(2)) ;
                if fenmu > 1e-12
                    xxx = wavek * fenmu ;
                    bj11 = Bessel_J1(xxx) ;
                    wjmtemp = bj11 / fenmu ;
                else
                    wjmtemp = 0.5 * wavek ;
                end
                sum2 = sum2 + wjmtemp * fi ;
            end
            sum2 = 0.5 * sum1 + subh * sum2 ;
    end
    BMIM(1,j) = -sum2 * 0.5 * wavek ;
end

%  2.5 Compute the hypersingular part
BMRE2 = BI1_Algorithm_III(n, p, col, BMRE);

% BMRE2 = BI1_Algorithm_II(n, p, BMRE);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function BMRE2 = BI1_Algorithm_III(n, p, col, BMRE)
h = p(3) - p(2) ;
work=zeros(1,n);
BMRE2=zeros(1,n);
work(1) = 2 / (pi * h) ;
work(2) = -(1 - log(2)) / (pi * h) ;
if n >2
    for i = 3 : n
        work(i) = -log((i - 1)^2 / ((i - 1)^2 - 1)) / (pi * h) ;
    end
end
% for	i = 1 : n
%     for j= 1 : n
%         ii = i - j ;
%         if ii < 0 
%             ii = -ii ;
%         end
%         BMRE(i, j) = BMRE(i, j) + work(ii + 1) ;
%     end
% end
BMRE2 = BMRE + work ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function BMRE2 = BI1_Algorithm_II(n, p, BMRE)
h = p(3) - p(2);
a = p(1);
b = p(n+1);
np = n+1;
for i = 1: 1: n
    sum1 = 0.0;
    sum2 = 0.0;
    for j = 1:1:n
        if i ~= j
            sum1 = sum1 + 1.0/((j-i)*(j-i)*h);
            sum2 = sum2 + 1.0/(j-i);
            BMRE(i,j) = BMRE(i,j) - 1.0/((j-i)*(j-i)*h)/pi;
        end

    end
    sum1 = sum1 + 0.5/((0-i)*(0-i)*h)+0.5/((np-i)*(np-i)*h);
    sum2 = sum2+ 0.5/(0-i)+0.5/(np-i);
    BMRE(i,i) = BMRE(i,i)+(sum1+1.0/h+(b-a)/(b-p(i))/(p(i)-a))/pi;
    if i > 1
        BMRE(i,i-1)=BMRE(i,i-1)-(1+sum2 -log((b-p(i))/(p(i)-a)))/(2*pi*h);
    end
	if  i < n
      BMRE(i,i+1)=BMRE(i,i+1)-(1-sum2 + log((b-p(i))/(p(i)-a)))/(2*pi*h);
    end
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function  by1 = Bessel_mdf1Y1(x)
r = [-0.4900604943e13, 0.1275274390e13, -0.5153438139e11, 0.7349264551e9, -0.4237922726e7, 0.8511937935e4] ;
s = [0.2499580570e14, 0.4244419664e12, 0.3733650367e10, 0.2245904002e8, 0.1020426050e6, 0.3549632885e3, 1.0] ;
p = [1.0, 0.183105e-2, -0.3516396496e-4, 0.2457520174e-5, -0.240337019e-6] ;
q =	[0.04687499995, -0.2002690873e-3, 0.8449199096e-5, -0.88228987e-6, 0.105787412e-6] ;
if x < 8
    bj11 = Bessel_J1(x) ;
	y = x^2 ;
    t = x * (r(1) + y * (r(2) + y * ( r(3) + y * (r(4) + y * (r(5) + y * r (6)))))) ;
	by1 = t / (s(1) + y * (s(2) + y * (s(3) + y * (s(4) + y * (s(5) + y * (s(6) + s(7) * y)))))) ;
    if x >1e-12
        by1 = by1 + 0.636619772367581 * (log(x) * (bj11 - 0.5 * x)) ;
    end
else
    z = 8 / x ;
    y = z^2 ;
    xx = x - 2.35619449019234 ;
    t1 = sqrt(0.636619772367581 / x) * sin(xx)*(p(1)+y*(p(2)+y*(p(3)+y*(p(4)+y*p(5))))) ;
    t2 = z*cos(xx)*(q(1)+y*(q(2)+y*(q(3)+y*(q(4)+y*q(5))))) ;
    by1 = t1 + t2 +(2/x-x*log(x))/pi ;
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function by1 = Bessel_mdf2Y1(x) 
r = [-0.4900604943e13,0.1275274390e13,-0.5153438139e11,0.7349264551e9,-0.4237922726e7,0.8511937935e4] ;,	
s = [0.2499580570e14,0.4244419664e12,0.3733650367e10,0.2245904002e8,0.1020426050e6, 0.3549632885e3,1.0] ;
p = [1.0, 0.183105e-2, -0.3516396496e-4, 0.2457520174e-5, -0.240337019e-6] ;
q =	[0.04687499995, -0.2002690873e-3, 0.8449199096e-5, -0.88228987e-6, 0.105787412e-6] ;
if x < 8
    bj11 = Bessel_J1(x) ;
	y=x^2 ;
    t = x*(r(1)+y*(r(2)+y*(r(3)+y*(r(4)+y*(r(5)+y*r(6)))))) ;
    by1=t/(s(1)+y*(s(2)+y*(s(3)+y*(s(4)+y*(s(5)+y*(s(6)+s(7)*y)))))) ;
    if x>1e-12
        by1 = by1 + 0.636619772367581*(log(x)*(bj11)) ;
    end
else
    z = 8 / x ;
	y = z^2 ;
    xx = x - 2.35619449019234 ;
	t1 = sqrt(0.636619772367581/x) * sin(xx) * (p(1)+y*(p(2)+y*(p(3)+y*(p(4)+y*p(5))))) ;
    t2 = z*cos(xx)*(q(1)+y*(q(2)+y*(q(3)+y*(q(4)+y*q(5))))) ;
    by1 = t1 + t2 + 0.636619772367581/x ;
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function by1 = Bessel_Y1(x) 
r = [-0.4900604943e13,0.1275274390e13,-0.5153438139e11,0.7349264551e9,-0.4237922726e7,0.8511937935e4];
s = [0.2499580570e14,0.4244419664e12,0.3733650367e10,0.2245904002e8,0.1020426050e6,0.3549632885e3,1.0] ;
p = [1.0,0.183105e-2,-0.3516396496e-4,0.2457520174e-5,-0.240337019e-6] ;
q = [0.04687499995,-0.2002690873-3,	0.8449199096e-5,-0.88228987e-6,0.105787412e-6] ;
bj1 = Bessel_j1(x) ;
if x < 8
    y=x^2 ;
    t = x*(r(1)+y*(r(2)+y*(r(3)+y*(r(4)+y*(r(5)+y*r(6)))))) ;
    by1 = t / (s(1)+y*(s(2)+y*(s(3)+y*(s(4)+y*(s(5)+y*(s(6)+s(7)*y))))))+2/pi*(log(x)*bj1-1/x) ;
else
    z = 8 / x ;
	y = z^2 ;
    xx = x - 0.75 * pi ;
    t = sin(xx)*(p(1)+y*(p(2)+y*(p(3)+y*(p(4)+y*p(5)))))+z*cos(xx)*(q(1)+y*(q(2)+y*(q(3)+y*(q(4)+y*q(5))))) ;
	by1 = sqrt(2/x/pi) * t ;
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function bj1 = Bessel_J1(x) 
r = [72362614232.0, -7895059235.0, 242396853.1,-2972611.439,15704.48260,-30.16036606] ;	
s = [144725228442.0, 2300535178.0, 18583304.74,99447.43394,376.9991397,1.0] ;
p = [1.0,0.183105e-2,-0.3516396496e-4, 0.2457520174e-5, -0.240337019-6] ;
q =	[0.04687499995, -0.2002690873e-3, 0.8449199096e-5,-0.88228987e-6,0.105787412e-6] ;
ax = abs(x) ;
if ax<8
    y = x^2 ;
    t =	x*(r(1)+y*(r(2)+y*(r(3)+y*(r(4)+y*(r(5)+y*r(6)))))) ;
	bj1 = t / (s(1)+y*(s(2)+y*(s(3)+y*(s(4)+y*(s(5)+y*s(6)))))) ;
else
    z = 8 / x ;
    y = z^2 ;
    xx = ax - 0.75 * pi ;
	t1 = sqrt(2/ax/pi) ;
    t2 = cos(xx)*(p(1)+y*(p(2)+y*(p(3)+y*(p(4)+y*p(5)))))- z*sin(xx)*(q(1)+y*(q(2)+y*(q(3)+y*(q(4)+y*q(5))))) ;
	bj1 = t1 * t2 * sign(x) ;
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%    Return the piece linear basis function phi_m(x) for x in [a, b].
function fi = phi(m,n,p,x) 
if x > p(m) && x < p(m+1)
    fi=(x-p(m))/(p(m+1)-p(m)) ;
else
    if x >= p(m+1) && x < p(m+2)
        fi = (p(m+2)-x)/(p(m+2)-p(m+1)) ;	    
	else
		fi = 0 ;
    end
end