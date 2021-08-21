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
## @deftypefn {} {@var{retval} =} Application_cas_2 (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Jérémy <Jérémy@LAPTOP-67JK3OLE>
## Created: 2021-06-14

function retval = Application_cas_2 (input1, input2)
clear all
clc

L=2; %Longueur de la barre
N=100; %Nombre de points
a=0.1; %alpha
b=0; %beta (condition de neumann)
Tc=100;
Tf=50;
dx=L/(N-1); %Pas d'espace
dt=0.01; %Pas de temps,
%Pas de temps choisi par rapport a la condition de stabilite
% sur le schema explicite.
r=(a^2*dt)/((dx)^2); 

x=0:dx:L; %Intervalle espace
t=0:dt:60; %Intervalle temps

T0=zeros(length(x),1); %Temperatures initiales dans les deux barres
T0(1:N/2,1)=Tc;
T0(N/2:end,1)=Tf;

%SCHEMA DE CRANK NICHOLSON

Ccn=zeros(N-1,1);
Ccn(1,1)=-2*r*b*dx; %Prise en compte de la condition de Neumann
Ccn(N-1,1)=2*r*b*dx; %prise en compte de la condition de Neumann

A1cn=diag(ones(N-1,1).*(1+r),0);
A2cn=diag(ones(N-2,1).*(-r/2),-1);
A3cn=diag(ones(N-2,1).*(-r/2),1);
Acn=A1cn+A2cn+A3cn;
Acn(1,2)=-r; %Prise en compte de la condition de Neumann
Acn(end,end-1)=-r; %Prise en compte de la condition de Neumann

B1cn=diag(ones(N-1,1).*(1-r),0);
B2cn=diag(ones(N-2,1).*(r/2),-1);
B3cn=diag(ones(N-2,1).*(r/2),1);
Bcn=B1cn+B2cn+B3cn;
Bcn(1,2)=r;%Prise en compte de la condition de Neumann
Bcn(end,end-1)=r; %Prise en compte de la condition de Neumann

Ainvcn=inv(Acn);

Tcn=zeros(N-1,length(t));
Tcn(:,1)=T0(1:N-1);

for j=2:length(t)
  Tcn(:,j)=Ainvcn*(Bcn*Tcn(:,j-1)+Ccn);
end


%GRAPHIQUES

xx=0:dx:L-dx; %Echelle d'espace pour les schemas de discretisation

figure 1
hold on
plot(xx,Tcn(:,1),'c')
legend('Crank-Nicholson');
title('t=0sec');
hold off

figure 2
hold on
plot(xx,Tcn(:,60/(3*dt)),'c')
legend('Crank-Nicholson');
title('t=20sec');
hold off

figure 3
hold on
plot(xx,Tcn(:,(2*60)/(3*dt)),'c')
legend('Crank-Nicholson');
title('t=40sec');
hold off

figure 4
hold on
plot(xx,Tcn(:,end),'c')
legend('Crank-Nicholson');
title('t=60sec');
hold off

endfunction
