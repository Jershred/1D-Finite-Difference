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
## @deftypefn {} {@var{retval} =} erreur_temporelle_cn (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Jérémy <Jérémy@LAPTOP-67JK3OLE>
## Created: 2021-06-07

function retval = erreur_temporelle_cn (input1, input2)
clear all
clc

L=2; %Longueur de la barre
a=0.1; %alpha
b=0; %beta
y=0; %gamma
N=100;%Nombre de points
dx=L/(N-1);
x=0:dx:L;
dt=[0.001,0.005,0.01,0.015];

Err_cn=zeros(1,length(N)); %Erreur de la methode C-N pour differents dt

for i=1:4
  %Chaque iteration correspond a un dt different
  
  r=(a^2*dt(1,i))/((dx)^2); 
  t=0:dt(1,i):60; %Intervalle temps
  
  %SOLUTION EXACTE
  Texacte_t=exp((-(a^2)*((pi^2)/(4*(L^2)))).*t);
  Texacte_x=cos((pi/(2*L)).*x);
  Texacte=Texacte_x(end-1)*Texacte_t(end); 
  % On prend l'avant derniere valeur de Texacte_x car 
  %elle correspond � la derniere valeur en x de Tcn
  
  %SCHEMA DE CRANK NICHOLSON
  T0cn=cos(pi*x/(2*L));
  T0cn=T0cn';

  Ccn=zeros(N-1,1);
  Ccn(1,1)=-2*r*b*dx; %Prise en compte de la condition de Neumann
  Ccn(N-1,1)=r*y; %prise en compte de la condition de Dirichlet

  A1cn=diag(ones(N-1,1).*(1+r),0);
  A2cn=diag(ones(N-2,1).*(-r/2),-1);
  A3cn=diag(ones(N-2,1).*(-r/2),1);
  Acn=A1cn+A2cn+A3cn;
  Acn(1,2)=-r; %Prise en compte de la condition de Neumann

  B1cn=diag(ones(N-1,1).*(1-r),0);
  B2cn=diag(ones(N-2,1).*(r/2),-1);
  B3cn=diag(ones(N-2,1).*(r/2),1);
  Bcn=B1cn+B2cn+B3cn;
  Bcn(1,2)=r;%Prise en compte de la condition de Neumann

  Ainvcn=inv(Acn);

  Tcn=zeros(N-1,length(t));
  Tcn(:,1)=T0cn(1:N-1);

  for j=2:length(t)
    Tcn(:,j)=Ainvcn*(Bcn*Tcn(:,j-1)+Ccn);
  endfor
  
  %Calcul de l'erreur de la methode de C-N pour le dt de l'iteration en cours
  Err_cn(1,i)=abs((Texacte-Tcn(end,end))/Texacte);

endfor

logdt=log10(dt); %Log du pas de temps
logErr_cn=log10(Err_cn); %log de l'erreur de la methode de C-N

%Approximation affine du log de l'erreur:
Affine_Err_cn=polyfit(logdt,logErr_cn,1)

hold on
plot(logdt,logErr_cn,'kx')
plot(logdt,polyval(Affine_Err_cn,logdt),'g')
legend('Crank-Nicholson','Approximation affine')
title('Schema Crank-Nicholson')
hold off
endfunction
