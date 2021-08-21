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
## @deftypefn {} {@var{retval} =} erreur_spatiale_explicite (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: JÃ©rÃ©my <JÃ©rÃ©my@LAPTOP-67JK3OLE>
## Created: 2021-06-07

function retval = erreur_spatiale_explicite (input1, input2)
clear all
clc

L=2; %Longueur de la barre
a=0.1; %alpha
b=0; %beta
y=0; %gamma
dt=0.001; %Pas de temps
t=0:dt:30; %Intervalle temps

N=[20 50 100 150]';%Nombre de points
delta_x=zeros(1,length(N)); %dx pour le graphe final

%Erreur de la méthode exponentielle pour chaque série de point dans N.
Err_exp=zeros(1,length(N));

for i=1:4 
  %Chaque iteration correspond à un nombre différent 
  %de point pris en espace (vecteur N).
  Nx=N(i,1);
  dx=L/(Nx-1); %Pas d'espace
  delta_x(1,i)=dx; %A la fin, delta_x contiendra les 4 valeurs prisent par dx.
  r=(a^2*dt)/((dx)^2); 
  x=0:dx:L; %Intervalle espace
  
  %SOLUTION EXACTE
  Texacte_t=exp((-(a^2)*((pi^2)/(4*(L^2)))).*t);
  Texacte_x=cos((pi/(2*L)).*x);
  Texacte=Texacte_x(end-1)*Texacte_t(end); 
  %on prend les valeurs finales de Texacte pour le calcul de l'erreur.
  %On prend l'avant dernière valeur de Texacte_x car la valeur de x 
  %correspond à celle de la dernière valeur de Texplicit
  
  %SCHEMA EXPLICITE
  T0explicit=cos(pi*x/(2*L));
  T0explicit=T0explicit';

  Bexplicit=zeros(Nx-1,1);
  Bexplicit(1,1)=-2*r*b*dx; %Prise en compte de la condition de Neumann
  Bexplicit(Nx-1,1)=r*y; %Prise en compte de la condition de Dirichlet

  A1explicit=diag(ones(Nx-1,1).*(1-2*r),0);
  A2explicit=diag(ones(Nx-2,1).*r,-1);
  A3explicit=diag(ones(Nx-2,1).*r,1);
  A3explicit(1,2)=2*r; %Prise en compte de la condition de Neumann
  Aexplicit=A1explicit+A2explicit+A3explicit;

  Texplicit=zeros(Nx-1,length(t));
  Texplicit(:,1)=T0explicit(1:Nx-1);

  for j=2:length(t)
    Texplicit(:,j)=Aexplicit*Texplicit(:,j-1)+Bexplicit;
  endfor
  
  %Calcul de l'erreur de la methode explicite pour le dx de l'iteration en cours
  Err_exp(1,i)=abs((Texacte-Texplicit(end,end))/Texacte);

endfor

logdx=log10(delta_x); %Log du pas d'espace
logErr_exp=log10(Err_exp); %log de l'erreur de la methode explicite

%Approximation affine du log de l'erreur:
Affine_Err_exp=polyfit(logdx,logErr_exp,1)

hold on
plot(logdx,logErr_exp,'kx')
plot(logdx,polyval(Affine_Err_exp,logdx),'r')
legend('Explicite','Approximation affine')
title('Schema explicite')
hold off
endfunction
