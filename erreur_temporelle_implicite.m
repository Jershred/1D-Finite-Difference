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
## @deftypefn {} {@var{retval} =} erreur_temporelle_implicite (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Jérémy <Jérémy@LAPTOP-67JK3OLE>
## Created: 2021-06-07

function retval = erreur_temporelle_implicite (input1, input2)
clear all
clc

L=2; %Longueur de la barre
a=0.1; %alpha
b=0; %beta
y=0; %gamma
N=100;%Nombre de points
dx=L/(N-1);
x=0:dx:L;
dt=[0.001,0.005,0.010,0.015];

Err_imp=zeros(1,length(N)); %Erreur de la methode implicite pour differents dt

for i=1:4 
  %Chaque iteration correspond � un dt different
  
  r=(a^2*dt(1,i))/((dx)^2); 
  t=0:dt(1,i):60; %Intervalle temps
  
  %SOLUTION EXACTE
  Texacte_t=exp((-(a^2)*((pi^2)/(4*(L^2)))).*t);
  Texacte_x=cos((pi/(2*L)).*x);
  Texacte=Texacte_x(end-1)*Texacte_t(end); 
  %on prend les valeurs finales de Texacte pour le calcul de l'erreur.
  % On prend l'avant derniere valeur de Texacte_x car la valeur de x 
  %correspond � celle de la derniere valeur de Timplicit
  
  %SCHEMA IMPLICITE
  T0implicit=cos(pi*x/(2*L));
  T0implicit=T0implicit';

  Bimplicit=zeros(N-1,1);
  Bimplicit(1,1)=-2*r*b*dx; %Prise en compte de la condition de Neumann
  Bimplicit(N-1,1)=r*y; %Prise en compte de la condition de Dirichlet

  A1implicit=diag(ones(N-1,1).*(1+2*r),0);
  A2implicit=diag(ones(N-2,1).*(-r),-1);
  A3implicit=diag(ones(N-2,1).*(-r),1);
  Aimplicit=A1implicit+A2implicit+A3implicit;
  Aimplicit(1,2)=-2*r;%Prise en compte de la condition de Neumann
  Ainv=inv(Aimplicit);

  Timplicit=zeros(N-1,length(t));
  Timplicit(:,1)=T0implicit(1:N-1);

  for j=2:length(t)
    Timplicit(:,j)=Ainv*(Timplicit(:,j-1)+Bimplicit);
  endfor
 
 %Calcul de l'erreur de la methode implicite pour le dt de l'iteration en cours
  Err_imp(1,i)=abs((Texacte-Timplicit(end,end))/Texacte);

endfor

logdt=log10(dt); %Log du pas de temps
logErr_imp=log10(Err_imp); %log de l'erreur de la methode implicite

%Approximation affine du log de l'erreur:
Affine_Err_imp=polyfit(logdt,logErr_imp,1)

hold on
plot(logdt,logErr_imp,'kx')
plot(logdt,polyval(Affine_Err_imp,logdt),'m')
legend('Implicite','Approximation affine')
title('Schema implicite')
hold off
endfunction
