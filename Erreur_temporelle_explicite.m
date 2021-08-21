## Copyright (C) 2021 J√©r√©my
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
## @deftypefn {} {@var{retval} =} Erreur_temporelle_explicite (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: J√©r√©my <J√©r√©my@LAPTOP-67JK3OLE>
## Created: 2021-06-07

function retval = Erreur_temporelle_explicite (input1, input2)
clear all
clc

L=2; %Longueur de la barre
a=0.1; %alpha
b=0; %beta
y=0; %gamma
N=100;%Nombre de points
dx=L/(N-1);
x=0:dx:L;
dt=[0.001,0.002,0.003,0.004];
% Pour respecter les criteres de stabilite du schema explicite pour N=100,
% on choisit des valeurs de dt < 0.005.

%Erreur de la mÈthode exponentielle pour diffÈrents dt
Err_exp=zeros(1,length(N));

for i=1:4
  %Chaque iteration correspond ‡ un dt different
  
  r=(a^2*dt(1,i))/((dx)^2); 
  t=0:dt(1,i):60; %Intervalle temps
  
  %SOLUTION EXACTE
  Texacte_t=exp((-(a^2)*((pi^2)/(4*(L^2)))).*t);
  Texacte_x=cos((pi/(2*L)).*x);
  Texacte=Texacte_x(end-1)*Texacte_t(end); 
  % On prend l'avant derniere valeur de Texacte_x car 
  %elle correspond ‡ la derniere valeur en x de Texplicit
  
  %SCHEMA EXPLICITE
  T0explicit=cos(pi*x/(2*L));
  T0explicit=T0explicit';

  Bexplicit=zeros(N-1,1);
  Bexplicit(1,1)=-2*r*b*dx; %Prise en compte de la condition de Neumann
  Bexplicit(N-1,1)=r*y; %Prise en compte de la condition de Dirichlet

  A1explicit=diag(ones(N-1,1).*(1-2*r),0);
  A2explicit=diag(ones(N-2,1).*r,-1);
  A3explicit=diag(ones(N-2,1).*r,1);
  Aexplicit=A1explicit+A2explicit+A3explicit;
  Aexplicit(1,2)=2*r; %Prise en compte de la condition de Neumann

  Texplicit=zeros(N-1,length(t));
  Texplicit(:,1)=T0explicit(1:N-1);

  for j=2:length(t)
    Texplicit(:,j)=Aexplicit*Texplicit(:,j-1)+Bexplicit;
  endfor
  
  %Calcul de l'erreur de la mÈthode explicite pour le dt de l'iteration en cours
  %Err_exp(1,i)=abs((Texacte(end-1,end)-Texplicit(end,end))/Texacte(end-1,end));
  Err_exp(1,i)=abs((Texacte-Texplicit(end,end))/Texacte);
endfor

logdt=log10(dt); %Log du pas de temps
logErr_exp=log10(Err_exp); %log de l'erreur de la methode explicite

%Approximation affine du log de l'erreur:
Affine_Err_exp=polyfit(logdt,logErr_exp,1)

hold on
plot(logdt,logErr_exp,'kx')
plot(logdt,polyval(Affine_Err_exp,logdt),'r')
legend('Explicite','Approximation affine')
title('Schema explicite')
hold off

endfunction
