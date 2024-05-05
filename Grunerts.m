function [T, r] = Grunerts(u_hat, P)
%{
This function uses Grunert's method of PnP to find camera rotation matrix
and position
INPUTS:
- u = matrix of 4 rows cam LOS vectors. x, y, z as columns
- P = matrix of 4 rows world points. x, y, z as columns
OUTPUTS:
- T = attitude transformation matrix
- r = position of camera in world frame
%}
arguments
    u_hat (4,3) double {mustBeFinite,mustBeReal}
    P (4,3) double {mustBeFinite,mustBeReal}
end
P1 = P(1,:)';
P2 = P(2,:)';
P3 = P(3,:)';
u1 = u_hat(1,:)';
u2 = u_hat(2,:)';
u3 = u_hat(3,:)';

a = norm(P2-P3);
b = norm(P3-P1);
c = norm(P1-P2);
d = (a^2-c^2)/b^2;

alpha = acos(u2'*u3);
beta = acos(u1'*u3);
gamma = acos(u1'*u2);

A4 = (d-1)^2 - 4*(c^2/b^2)*(cos(alpha))^2;
A3 = 4*(d*(1-d)*cos(beta)-(1-(a^2+c^2)/b^2)*cos(alpha)*cos(gamma) + ...
    2*(c^2/b^2)*(cos(alpha)^2)*cos(beta));
A2 = 2*(d^2-1+2*d^2*cos(beta)^2 + 2*((b^2-c^2)/(b^2))*cos(alpha)^2 - ...
    4*((a^2+c^2)/(b^2))*cos(alpha)*cos(beta)*cos(gamma) + ...
    2*((b^2-a^2)/(b^2))*cos(gamma)^2);
A1 = 4*(-d*(1+d)*cos(beta) + 2*((a^2)/(b^2))*cos(gamma)^2*cos(beta) ...
    -(1-(a^2+c^2)/b^2)*cos(alpha)*cos(gamma));
A0 = (1+d)^2 - 4*((a^2)/(b^2))*cos(gamma)^2;

v = roots([A4,A3,A2,A1,A0]);
v = real(v(imag(v) == 0));
v = v(v>=0);

u = ((-1+d)*v.^2 - 2*d*cos(beta)*v + 1 + d).*(2*(cos(gamma)-v*cos(alpha))).^-1;

s1 = sqrt((a^2)*((u.^2+v.^2-2*u.*v*cos(alpha)).^-1));
s2 = u.*s1;
s3 = v.*s2;

minErr = 10^100; %Initializing minimum error
iMinErr = 0; %Index of v with minimum error
rMinErr = zeros(3,1);
TMinErr = zeros(3,3);

%Making E
e12MCMF = (P2-P1)/norm(P2-P1);
e23MCMF = (P3-P2)/norm(P3-P2);
e31MCMF = (P1-P3)/norm(P1-P3);
E = [e12MCMF, e23MCMF,e31MCMF];

for i = 1:length(v)
    %Making A
    e12cam = (s2(i)*u2-s1(i)*u1)/norm(s2(i)*u2-s1(i)*u1);
    e23cam = (s3(i)*u3-s2(i)*u2)/norm(s3(i)*u3-s2(i)*u2);
    e31cam = (s1(i)*u1-s3(i)*u3)/norm(s1(i)*u1-s3(i)*u3);
    A = [e12cam, e23cam, e31cam];
    B = A*E';%B to use SVD on
    [svdU,~,svdV]=svd(B);
    M = eye(3);
    M(3,3) = det(svdU)*det(svdV);
    TMCMF2cam = svdU*M*svdV';
    r = P1-TMCMF2cam'*s1(i)*u1;
    err = 0;
    for j = 1:length(u_hat)
        ui_exp = (TMCMF2cam*(P(j,:)'-r))/norm(TMCMF2cam*(P(j,:)'-r));
        err = err+(acos(ui_exp' * u_hat(j,:)'))^2;
    end
    if err <= minErr
        minErr = err; %Initializing minimum error
        iMinErr = i; %Index of v with minimum error
        rMinErr = r;
        TMinErr = TMCMF2cam;
    end
end

T = TMinErr;
r = rMinErr;
end