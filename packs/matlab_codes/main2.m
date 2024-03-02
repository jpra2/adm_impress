% Simulador para resolver a equacao de difusao em 2D 
% Desenvolvedor: Prof. Fernando R.L. Contreras
% 
%% Este codigo roda somente escoamento monofasico
clear all
clc
format long
global coord centelem elem esurn1 esurn2 nsurn1 nsurn2 bedge inedge ...
    normals esureface1 esureface2 esurefull1 esurefull2 elemarea dens ...
    visc satlimit pormap bcflag courant totaltime filepath foldername gravitational...
    benchmark pmetodo interpol iteration erromethod correction strategy...
    typecorrection jacob;
%%========================================================================%

[coord,centelem,elem,esurn1,esurn2,nsurn1,nsurn2,bedge,inedge,normals,...
    esureface1,esureface2,esurefull1,esurefull2,elemarea,dens,visc,...
    satlimit,pormap,bcflag,courant,totaltime,filepath,foldername,kmap,...
    wells] = preprocessor;

%% NOTAS E PENDENCIAS
% 1. Para o interpolador com correcao de pontos harmonicos precisa ainda
% implementar o caso artigo Zhang Kobaise figura 12. 
% 3. Deve-se investir nos precondicionadores
% 4. Deve-se incluir condicao de contorno de Neumann no eLPEW2
% 7. O termo gravitacional consistente necessita trabalhar
%% Habilite esta funcao para obter distocao de malhas estruturadas
%[auxcoord]=distortedramd;
%% funcao que modificacao de bedge
% esta funcao deve ser ativado para malhas patologicas
%[bedge]=modificationbedge(bedge);

%% calculo o flag do elemento que deseja
%   a=6287;
%   b=445;
%   c=5740;
%   d=0;
%   [elemento]=searchelement(a,b,c,d)
%% escolha o tipo de erro discreto que deseja usar
% erromethod1 ---> erro utilizado por Gao e Wu 2010;  amplamente utilizado
% erromethod2 --->  ''     ''     por Lipnikov et al 2010
% erromethod3 --->  ''     ''     por Eigestad et al 2005
% erromethod4 --->  ''     ''     por Shen e Yuan 2015
% erromethod6 --->  ''     ''     por M. Starnoni 2019, para o caso gravitacional
erromethod='erromethod1';
%% defina o tipo de metodo  
% tpfa      --> (TPFA)
% mpfad     --> (MPFA-D) 
% lfvLPEW   --> (MPFA-HD)ou (MPFA-QL)
% lfvHP     --> (MPFA-H)
% lfvEB     --> metodo baseado na face (MPFA-BE), ainda os testes nao foram feitos
% nlfvLPEW  --> (NL-TPFA)
% nlfvDMPSY --> (NL-DMP) (Gao e Wu, 2013) e (Sheng e Yuan, 20...)
% nlfvHP    --> (NL-TPFA-H) metodo nao linear baseado em pontos harmonicos
% nlfvPPS   --> 
% interpfree
pmetodo='nlfvLPEW';
%% metodo de interacao: picard, newton, broyden, secant,
% m�todo de iterecao proprio de m�todos n�o lineares iterfreejacobian,iterdiscretnewton, JFNK
% iteration='iterdiscretnewton';
 iteration='iterbroyden';
% iteration='JFNK';
% iteration='fullpicard';
% iteration='MPE'; 
% iteration='RRE'; % picard com acelerador rank reduced extrapolation
%  iteration='AA';  % picard com aceleracao de Anderson
%iteration='iterhybrid';
%iteration='fsolver';
%% Para metodo nao-linear e iteracao Broyden ou Newton
% escolha a estrategia de calculo do jacobiano
%jacob='classic';
jacob='nonclassic';
%% qual tipo de interpolacao deseja utilizar
interpol='LPEW2';
%interpol='LPEW1';
%interpol='LS';
%interpol='eLS';
%interpol='eLPEW2';
%% correcao dos pontos harmonicos
% digite 'yes' ou 'no'
correction='no';
% qual tipo de correcao deseja utilizar
typecorrection='firstcorrection'; % correcao utilizando express. simplif.
%typecorrection= 'secondcorrection'; % correcao utilizando ponto medio da face
%typecorrection='thirdcorrection'; % correcao utilizando metodo Kobaise
%% digite segundo o benchmark
% procure o teste que deseja rodar no arquivo "benchmarks.m"
benchmark='shenyuan16';
%benchmark='gaowu5'; 
%benchmark='starnonigrav1';
%% com termo gravitacional
% com termo gravitacional 'yes' ou 'no'
gravitational='no';
% quando pretende incluir termo gravitacional deve utilizar estrategia
% 'starnoni' ou 'inhouse' ou 'inhouse1'
strategy= 'inhouse';
%strategy='GravConsist';

%% adequacao das permeabilidades e otros parametros fisico-geometricos 
%segundo cada caso ou problema
[elem,kmap,normKmap,pressurexact,bedge,fonte,velexact,gravelem,gravno,...
    gravface,grav_elem_escalar]=benchmarks(kmap,elem,bedge);
%   mm=find(bedge(:,4)==202);% ???
%   bedge(mm',4)=201; %???
% F faces na vizinhanca de um elemento
% V 
% N
[F,V,N]=elementface;

%% pre-processador local
[pointarmonic,parameter,gamma,p_old,tol,nit,nflagface,nflagno,...
    weightDMP,Hesq,Kde,Kn,Kt,Ded,auxface,calnormface,gravresult,gravrate,weight,s,wg]=...
    preprocessorlocal(kmap,N,gravelem,gravface);
