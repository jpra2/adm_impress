import os
# parent_dir = os.path.dirname(__file__)

# os.system('sudo docker run -it -v ' + parent_dir + ':/pytest elliptic_scy:latest -c "cd/pytest; bash"')
# os.system('sudo docker run -it -v ' + parent_dir + ':/pytest elliptic_scy:latest bash -c "cd /pytest/preprocessor; python3 definicoes.py"')

# sudo docker run -it -v  $PWD:/pytest elliptic_scy:latest bash -c "cd /pytest/preprocessor; bash"
# os.system('sudo docker run -it -v  $PWD:/pytest desenvolvimento-scipy:latest bash -c "cd /pytest; bash"')
# os.system('sudo docker run -it -v  $PWD:/pytest gabrielmmats/impress:latest bash -c "cd /pytest; bash"')
os.system('sudo docker run -it -v  $PWD:/pytest facsa/impress-workspace:latest bash -c "cd /pytest; bash"')

# profiling
# os.system('python -m cProfile -s cumtime definicoes.py')
