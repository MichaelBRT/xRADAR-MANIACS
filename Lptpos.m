%% Lagrange point location

function [x,y] = Lptpos(mu, Lpt)
% function [x,y] = Lptpos(mu, Lpt)
%
% Calculate location of libration point Lpt
%
% see gammaL() 

mu2=1-mu;

if Lpt < 4
	LPTS = [(mu2-gammaL(mu,1)) (mu2+gammaL(mu,2)) (-mu-gammaL(mu,3))];
	x=LPTS(Lpt);
	y=0;
else
	x=.5-mu;
	if Lpt == 4
		y=.5*sqrt(3);
	elseif Lpt==5
		y=-.5*sqrt(3);
	end
end

end

%%%%%%%%%%%%%


%% gamma_L, the distance from the Lagrange point Lpt to the nearest primary mass

function [gamma]=gammaL(mu, Lpt)
% gamma=gammaL(mu, Lpt);
%
% Calculate ratio of libration point distance from closest primary to distance
% between two primaries  (example gammaL1 = (E-L1)/AU)
%
%----------------------------------------------------------------------
% CR3BP (Circular Restricted Three-Body [Gravitational] Problem)
% with the LARGER MASS, M1 (i.e., Sun) to the left of the origin at (-mu,0)
% and the smaller mass, M2, or the planet (i.e., Earth), is at (1-mu, 0)
%
%      (rotating frame)
%
%                 L4
%
% -L3------M1--+-----L1--M2--L2-
%
%                 L5
%			 |<-->|
%			  gamma
%
% Shane Ross (revised 7.8.97)

mu2 = 1 - mu;

poly1 = [1  -1*(3-mu)  (3-2*mu)  -mu   2*mu  -mu ];
poly2 = [1     (3-mu)  (3-2*mu)  -mu  -2*mu  -mu ];
poly3 = [1     (2+mu)  (1+2*mu)  -mu2 -2*mu2 -mu2];

rt1 = roots(poly1); rt2 = roots(poly2); rt3 = roots(poly3);

for k=1:5
        if isreal(rt1(k)) GAMMAS(1)=rt1(k); end
        if isreal(rt2(k)) GAMMAS(2)=rt2(k); end
        if isreal(rt3(k)) GAMMAS(3)=rt3(k); end
end
gamma=GAMMAS(Lpt);

end