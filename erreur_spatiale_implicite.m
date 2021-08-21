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
## @deftypefn {} {@var{retval} =} erreur_spatiale_implicite (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Jérémy <Jérémy@LAPTOP-67JK3OLE>
## Created: 2021-06-07

function retval = erreur_spatiale_implicite (input1, input2)
clear all
clc

L=2; %Longueur de la barre
a=0.1; %alpha
b=0; %beta
y=0; %gamma
dt=0.001; %Pas de temps,
t=0:dt:30; %Intervalle temps

N=[20 50 100 150]';%Nombre de points
delta_x=zeros(1,length(N)); %dx pour le graphe final

%Erreur de la methode impliite pour chaque serie de point dans N
Err_imp=zeros(1,length(N));

for i=1:4 
  %Chaque iteration correspond � un nombre different
  %de point pris en espace (vecteur N).
  Nx=N(i,1);
  dx=L/(Nx-1); %Pas d'espace
  delta_x(1,i)=dx; %A la fin, delta_x contiendra les 4 valeurs prisent par dx.
  r=(a^2*dt)/((dx)^2); 
  x=0:dx:L; %Intervalle espace
  
  %SOLUTION EXACTE
  Texacte_t=exp((-(a^2)*((pi^2)/(4*(L^2)))).*t);
  Texacte_x=cos((pi/(2*L)).*x);
  Texacte=zeros(length(x),length(t));
  for j=1:length(t)
    for k=1:length(x)
      Texacte(k,j)=Texacte_x(k)*Texacte_t(j);
    endfor
  endfor
  
  %SCHEMA IMPLICITE
  T0implicit=cos(pi*x/(2*L));
  T0implicit=T0implicit';

  Bimplicit=zeros(Nx-1,1);
  Bimplicit(1,1)=-2*r*b*dx; %Prise en compte de la condition de Neumann
  Bimplicit(Nx-1,1)=r*y; %Prise en compte de la condition de Dirichlet

  A1implicit=diag(ones(Nx-1,1).*(1+2*r),0);
  A2implicit=diag(ones(Nx-2,1).*(-r),-1);
  A3implicit=diag(ones(Nx-2,1).*(-r),1);
  Aimplicit=A1implicit+A2implicit+A3implicit;
  Aimplicit(1,2)=-2*r;%Prise en compte de la condition de Neumann
  Ainv=inv(Aimplicit);

  Timplicit=zeros(Nx-1,length(t));
  Timplicit(:,1)=T0implicit(1:Nx-1);

  for j=2:length(t)
    Timplicit(:,j)=Ainv*(Timplicit(:,j-1)+Bimplicit);
  endfor
 
  %Calcul de l'erreur de la methode implicite pour le dx de l'iteration en cours
  Err_imp(1,i)=abs((Texacte(end-1,end)-Timplicit(end,end))/Texacte(end-1,end));

endfor

logdx=log10(delta_x); %Log du pas d'espace
logErr_imp=log10(Err_imp); %log de l'erreur de la methode implicite

%Approximation affine du log de l'erreur:
Affine_Err_imp=polyfit(logdx,logErr_imp,1)

hold on
plot(logdx,logErr_imp,'kx')
plot(logdx,polyval(Affine_Err_imp,logdx),'m')
legend('Implicite','Approximation affine')
title('Schema implicite')
hold off

endfunction
