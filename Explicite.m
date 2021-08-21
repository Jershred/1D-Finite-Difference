## Copyright (C) 2021 JÃ©rÃ©my
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
## @deftypefn {} {@var{retval} =} Explicite (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: JÃ©rÃ©my <JÃ©rÃ©my@LAPTOP-67JK3OLE>
## Created: 2021-05-14

function [T1,Tex,x] = Explicite (alpha,beta,gamma,L,tau,N,dx,tf,dt,Nt,r)
clear all
close all

%Paramètres physiques
alpha=0.1;
beta=0;
gamma=0;
L=2;
tau=(4*L**2)/(pi**2*alpha**2);

%Discrétisation spatiale
N=10;   %nombre de points de discrétisation
dx=L/(N-1);   %dx^2 > 2*alpha^2*dt

%Discrétisation temporelle
tf=5*tau;   %temps final (on prend t=5tau)
dt=4;
Nt=floor(tf/dt);  %nombre de pas de temps

%Calcul de r
r=(alpha^2)*dt/(dx^2);

%Calcul de T^n pour n=0 à partir de CI
x=0:dx:L

T1=zeros(N-1,Nt+1);   %solution numérique explicite

T0=cos(pi*x/(2*L));   %condition initiale

T1(:,1)=T0(1:N-1);
Tex(:,1)=cos(pi*x/(2*L))';

%Construction de A et B

v0=ones(1,N-1);
v1=ones(1,N-2);

B=zeros(N-1,1);
B(1)=-2*r*beta*dx;
B(N-1)=r*gamma;

%Résolution avce le schéma explicite
A=diag((1-2*r)*v0)+diag(r*v1,1)+diag(r*v1,-1);
A(1,2)=2*r;

%Boucle de temps
for i=1:Nt
  T1(:,i+1)=A*T1(:,i)+B;
  Tex(:,i+1)=cos(pi*x/(2*L))*exp(-alpha^2*pi^2*i*dt/(4*L^2));
end

T1=[T1;gamma*ones(1,Nt+1)];   %concaténation de la dernière ligne

figure;
plot(x,Tex(:,floor(tau/dt)),x,T1(:,floor(tau/dt)));
%on regarde le comportement à t=tau, la fin du régime transitoire
legend('Solution exacte','Solution explicite','Solution implicite', 'solution de Crank-Nicholson');
xlabel('x');
ylabel('T');

endfunction