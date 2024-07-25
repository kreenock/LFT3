import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_fwf("mooninfo_2024.txt")
print(df.columns)
plt.plot(df['Elon'], df['Elat'])
plt.show()