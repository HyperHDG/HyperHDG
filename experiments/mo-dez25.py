import matplotlib.pyplot as plt
import pandas as pd
import sys

paths = sys.argv[1:]
print(paths)

for path in paths:
  dt = pd.read_csv(path)
  ppath = path.split("-")[1].strip("s")
  print(ppath)
  plt.plot(dt.it, dt.res, label=ppath)

plt.legend()
#plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.yscale("log")
plt.xlabel("iterations")
plt.ylabel("relative residual")
plt.savefig("mo-dez25.png")
plt.show()
