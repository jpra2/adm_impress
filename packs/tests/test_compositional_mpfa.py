# TODO implementar o composicional mpfa

## Passo0
"""
1 - Definir a quantidade de componentes
2 - Definir a quantidade de fases
3 - Definir as propriedades iniciais do reservatório:
    composicao inicial
    permeabilidade
4 - Definir as condicoes de contorno:
    Pressao prescrita no contorno
    Fluxo prescrito no contorno
    composicao do contorno
5 - Definir os pocos se houverem:
    Pressao Prescrita
    Fluxo Prescrito
    composicao do poco caso o poco seja injetor
6 - Definir o tipo de permeabilidade relativa
7 - Definir a Temperatura
"""

## Passo 1
"""
1 - Calcular as propriedades das fases:
    Verificar a composição:
        Calculo de estabilidade
        Calculo de flash
    viscosidade
    densidade massica
    densidade molar
    permeabilidade relativa
2 - Calcular os pesos de cada no da malha
3 - Definir o valor do parametro do fluxo total e fracionario em cada edge
4 - Montar a matriz de transmissibilidade
"""