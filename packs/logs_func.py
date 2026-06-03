import os
import shutil
from rich import print as rprint
import numpy as np
from loguru import logger
# from rich.pretty import d

from packs import defnames, defpaths
from packs.io_operations import save_solver_times

def preparar_ambiente_logs(caminho_pasta: str, **kwargs):
    if os.path.exists(caminho_pasta):
        shutil.rmtree(caminho_pasta)
    os.makedirs(caminho_pasta)

    logger.remove() # Silencia o terminal completamente
    rprint(f'[green]Ambiente de logs preparado em: {caminho_pasta}[/green]')

def create_generic_logger(logs_path: str, tag: str, logfile: str, **kwargs):
    nome_arquivo = os.path.join(logs_path, tag, logfile)
    caminho_pasta_tag = os.path.dirname(nome_arquivo)
    preparar_ambiente_logs(caminho_pasta_tag)
    
    logger.add(
        nome_arquivo,
        rotation="10 MB", # Tamanho máximo do arquivo
        retention=0,      # Apaga o log antigo ao rotacionar, mantendo apenas o atual
        enqueue=True,     # Performance assíncrona
        level="DEBUG",
        filter=lambda record: record["extra"].get("tag") == tag,
        format="{time:HH:mm:ss} | {level: <8} | {message}"
    )
    
    return logger.bind(tag=tag)

def criar_logger_tag(logs_path: str, tag: str, **kwargs):
    logfile = 'solver.log'
    return create_generic_logger(logs_path, tag, logfile)

def create_operator_logger(op_logs_path: str='', tag: str='', **kwargs):
    logfile = 'operator.log'
    if op_logs_path == '':
        op_logs_path = defpaths.DEFAULT_OPERATOR_LOG_PATH
    if tag == '':
        tag = 'OP'
    return create_generic_logger(op_logs_path, tag, logfile)

def check_estag_residuo(residuo: list, logger, it: int, **kwargs):
    taxa_convergencia = residuo[it-1] / residuo[it-2]
    if 0.99999 < taxa_convergencia <= 1.0:
        logger.warning(f"Estagnação detectada: taxa de {taxa_convergencia:.3f}, iteração {it}")
        return True
    return False

def check_divergencia(residuo: list, logger, it: int, divergence_tolerance: float, **kwargs):
    resp = False
    if residuo[it-1] > divergence_tolerance:
        resp = True
    elif np.isnan(residuo[it-1]):
        resp = True
    elif np.isinf(residuo[it-1]):
        resp = True
    if resp == True:
        logger.warning(f"Divergência detectada. Residuo final {residuo[it-1]:.3f}, iteração {it}.")
    return resp

def check_divergencia_and_estag(residuo: list, logger, it: int, divergence_tolerance: float, **kwargs):
    return  check_divergencia(residuo, logger, it, divergence_tolerance=divergence_tolerance, **kwargs) or check_estag_residuo(residuo, logger, it, **kwargs)

def check_stopping_criteria(residuo: list, logger, it: int, divergence_tolerance: float, iteration_check: int, **kwargs):
    if not (it & iteration_check-1):
        return check_divergencia_and_estag(residuo, logger, it, divergence_tolerance, **kwargs)
    else: 
        return False
    
def plot_convergence_check(logger, it: int, residuo: float, siga: bool, maxit: int, **kwargs):
    if it > maxit:
        msg = f'Atingiu Limite de iteracoes. It: {it}, Maxit: {maxit}, residuo: {residuo}'
        logger.warning(msg)
    if siga:
        msg = f'Convergencia atingida. It: {it}, residuo: {residuo}'
        logger.success(msg)
    
def init_log_local(funcname, logs_path=defpaths.DEFAULT_LOG_PATH, tag=defnames.DEFAULT_TAG, **kwargs):
    
    log_local = criar_logger_tag(logs_path, tag)
    log_local.info(f'Iniciando {funcname}')
    return log_local