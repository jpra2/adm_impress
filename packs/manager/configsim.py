import matplotlib.pyplot as plt
from packs import defpaths
from pathlib import Path
from typing import List
import numpy as np
import yaml

DEFAULT = ''


class GlobalData:
    time_funcs = {}  # Dictionary to store timing information for functions

    
    def generate_time_plot(self, prefix_folder: str=''):
        filename = Path('barplot_functions_times.png')
        arq = str(Path(defpaths.flying) / prefix_folder / filename)

        plt.clf()
        plt.rcParams['text.usetex'] = True
        plt.rcParams['legend.fontsize'] = 14
        plt.rcParams['axes.labelsize'] = 16
        plt.rcParams["font.family"] = "Times New Roman"
        fig = plt.figure()
        ax = fig.add_subplot()

        text = list(self.time_funcs.keys())
        values = list(self.time_funcs.values())
        values = [round(v, 2) for v in values]

        bars = ax.bar(text, values, color='lightcoral')

        ax.set_ylabel('CPU Time (s)')
        ax.set_title('Function Execution Times')
        plt.xticks(rotation=45, ha='right')
        # for bar in bars:
        #     height = bar.get_height()
        #     ax.annotate(f'{height:.2f}',
        #                 xy=(bar.get_x() + bar.get_width() / 2, height),
        #                 xytext=(0, 3),  # 3 points vertical offset
        #                 textcoords="offset points",
        #                 ha='center', va='bottom')
        fig.tight_layout()
        fig.savefig(arq, dpi=500)
        plt.close()

    def generate_cumulative_plot(self, from_data:List[List[str]]=[], labels: List[str]=[], alias_from_data: List[str]=[], prefix_folder: str=''):
        filename = Path('barplot_functions_times_cumulative.png')
        arq = str(Path(defpaths.flying) / prefix_folder / filename)
        
        grupos = []
        for names in from_data:
            grupo = np.array([self.time_funcs.get(name, 0) for name in names])
            grupos.append(grupo)
        grupos = np.array(grupos)
        
        plt.clf()
        plt.rcParams['text.usetex'] = True
        plt.rcParams['legend.fontsize'] = 18
        plt.rcParams['axes.labelsize'] = 18
        plt.rcParams["font.family"] = "Times New Roman"
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot()
        
        soma = grupos[0].copy()
        soma[:] = 0
        for i, grupo in enumerate(grupos):
            ax.bar(alias_from_data, grupo, bottom=soma, label=labels[i])
            soma[:] += grupo
        
        # ax.set_title('Tempo', fontsize=14)
        plt.xticks(rotation=45, ha='right', fontsize=18)
        # ax.set_xlabel('Categorias')
        ax.set_ylabel('CPU Time (s)')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 1.0, box.height])
        # Put a legend to the right of the current axis
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        # plt.legend()
        # ax.margins(0.5)
        ax.set_ylim(0, soma.max() + 1)
        fig.tight_layout()
        fig.savefig(arq, dpi=500)
        plt.close()
        
    def file_to_export_times(self, path: str='') -> Path:
        # return str(Path(defpaths.flying) / 'functions_times.yaml')
        
        if path == '':
            return Path(defpaths.flying) / 'functions_times_generic.yaml'
        else:
            if path.endswith('.yaml') or path.endswith('.yml'):
                return Path(defpaths.flying) / path
            else:
                raise ValueError("The file extension must be .yaml or .yml")
        
    def load_times_from_file(self, path: str) -> dict:
        file_path = self.file_to_export_times(path)
        if file_path.is_file():
            with open(file_path, 'r') as file:
                self.time_funcs = yaml.safe_load(file)
    
    def export_times(self, path: str):
        file_path = self.file_to_export_times(path)
        with open(file_path, 'w') as file:
            yaml.dump(self.time_funcs, file, default_flow_style=False)
               
    def update_cumulative_times(self, funcname: str, path: str=''):
        
        elapsed_time = self.time_funcs.get(funcname, 0)
        time2 = self.time_funcs.get(funcname + '_cum', 0)
        self.time_funcs.update({funcname + '_cum': elapsed_time + time2})
        self.time_funcs.update({funcname: 0})
        self.export_times(path)
    
gdata = GlobalData()
