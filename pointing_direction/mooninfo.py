import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

a = [-5, 31, 50, 9, 73, 20]
b = [60, 31, 25, 4, -7, 3]

print(np.where(np.array(a) >= 30)[0])
print(np.where(np.array(b) >= 30)[0])
