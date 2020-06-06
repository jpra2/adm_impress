import time
import concurrent.futures

def do(a, b):
    t0=time.time()
    print('doing something')
    time.sleep(1)
    print(time.time()-t0)
    return a+b
t1=time.time()
aa=[]
with concurrent.futures.ProcessPoolExecutor() as executor:
    a=[1,2,5,6,7,8,9, 10,20,50,9,10,22,4,5,6,6,6,78,8]
    b=2
    results=[executor.submit(do,aa,b) for aa in a]
    for f in concurrent.futures.as_completed(results):
        aa.append(f.result())


print(time.time()-t1,'done all')
import pdb; pdb.set_trace()
