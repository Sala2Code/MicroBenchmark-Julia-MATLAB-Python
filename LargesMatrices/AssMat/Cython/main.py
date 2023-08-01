from ParamBeta import *

t=0
n=10
for i in range(n):
    t0 = time.time()
    main()
    t += time.time() - t0
t = t/n
print(t)
