import cProfile
import io
import pstats
from pstats import SortKey

def profile(func):
    def wrapper(*args, **kwargs):
        pr = cProfile.Profile()
        pr.enable()
        retval = func(*args, **kwargs)
        pr.disable()
        s = io.StringIO()
        sortby = SortKey.CUMULATIVE  # 'cumulative'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        print(s.getvalue())
        return retval

    return wrapper

## rodar o cprofile:
#ptyhon -m cProfile -o arquivo.profile arquivo.py

## visualizar o profile com snakeviz
#snakeviz arquivo.profile

## converter de profile para dot
#gprof2dot -f pstats arquivo.profile -o arquivo.dot

## converter de dot para svg
#dot -Tsvg arquivo.dot > arquivo.svg
#dot -Tpng arquivo.dot > arquivo.png
