import yaml
from pathlib import Path
from packs import defpaths

def load_solver_times(path: str=''):
    ext = '.yaml'
    if path == '':
        path3 = Path(defpaths.remove_folder) / 'solver_times.yaml'
    else:
        path3 = Path(path)
        if path3.suffix.lower() != ext.lower():
            path3 = path3.with_suffix(ext)
    
    path3.parent.mkdir(exist_ok=True, parents=True)
    path3.touch(exist_ok=True)
    
    times = yaml.safe_load(path3.read_text(encoding='utf-8')) or {}
    return times

def save_solver_times(times: dict, path: str=''):
    ext = '.yaml'
    if path == '':
        path3 = Path(defpaths.remove_folder) / 'solver_times.yaml'
    else:
        path3 = Path(path)
        if path3.suffix.lower() != ext.lower():
            path3 = path3.with_suffix(ext)
    
    data = load_solver_times(path=path3)
    data.update(times)
    
    conteudo = yaml.dump(data, allow_unicode=True)
    path3.write_text(conteudo, encoding='utf-8')