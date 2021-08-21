## Copyright (C) 2021 Jérémy
## 
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {} {@var{retval} =} Implicite (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Jérémy <Jérémy@LAPTOP-67JK3OLE>
## Created: 2021-05-14

function T2 = Implicite ()

clear all
close all

%Param�tres physiques
alpha=0.1;
beta=0;
gamma=0;
L=2;
tau=(4*L**2)/(pi**2*alpha**2);

%Discr�tisation spatiale
N=10;   %nombre de points de discr�tisation
dx=L/(N-1);   %dx^2 > 2*alpha^2*dt

%Discr�tisation temporelle
tf=5*tau;   %temps final (on prend t=5tau)
dt=4;
Nt=floor(tf/dt);  %nombre de pas de temps

%Calcul de r
r=(alpha^2)*dt/(dx^2);

%Calcul de T^n pour n=0 � partir de CI
x=0:dx:L

T2=zeros(N-1,Nt+1);   %solution num�rique implicitee
Tex=zeros(N,Nt+1);    %solution analytique

T0=cos(pi*x/(2*L));   %condition initiale

T2(:,1)=T0(1:N-1);
Tex(:,1)=cos(pi*x/(2*L))';

%Construction de A et B

v0=ones(1,N-1);
v1=ones(1,N-2);

B=zeros(N-1,1);
B(1)=-2*r*beta*dx;
B(N-1)=r*gamma;

%R�solution avce le sch�ma implicite
A=diag((1+2*r)*v0)+diag(-r*v1,1)+diag(-r*v1,-1);
A(1,2)=-2*r;

%Boucle de temps
for i=1:Nt
  T2(:,i+1)=A\(T2(:,i)+B);
  Tex(:,i+1)=exp(-alpha^2*pi^2*i*dt/(4*L^2))*cos(pi*x/(2*L));
end

T2=[T2;gamma*ones(1,Nt+1)];%concat�nation de la derni�re ligne

figure;
plot(x,Tex(:,floor(tau/dt)),x,T2(:,floor(tau/dt)));
%on regarde le comportement � t=tau, la fin du r�gime transitoire
legend('Solution exacte','Solution implicite');
xlabel('x');
ylabel('T');


endfunction