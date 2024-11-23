import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#Lendo os dados do arquivo .txt localizado na mesma pasta deste arquivo .py
tabela = pd.read_csv("gem_q00012.txt" , sep = "\t", header = None).to_numpy()
# sep : como os dados são separados. exp: ";" ou "\t" se for tabulação
# header: ignora cabeçario dos dados, 0 1 linha, None se não tiver.
#definindo as variáveis x e y que serão plotadas
#import pdb; pdb.set_trace()

x = tabela[:,0]
y = tabela[:,1]
#import pdb; pdb.set_trace()
#plot dos pontos
plt.scatter(x,y, color = 'black', label = 'Legenda')


#Títulos, legendas e salvando a figura
plt.legend(loc='best')

plt.xlabel('nome do eixo x', size = 10)

plt.ylabel('nome do eixo y', size = 10)

plt.title('Título do Gráfico', size = 15)

plt.savefig('nomedografico.png', dpi = 300)
plt.show()
