from bitmaplst import lst,lst2
import numpy as np
import matplotlib.pyplot as plt

data = np.array(lst, dtype=np.float32)
data2 = np.array(lst2, dtype= np.float32)

plt.imshow(
    data,
    cmap="gray",
    vmin=0.0,
    vmax=1.0,
    interpolation="nearest"
)

plt.colorbar()
plt.show()


plt.imshow(
    data2,
    cmap="gray",
    vmin=0.0,
    vmax=1.0,
    interpolation="nearest"
)

plt.show()

print("script done")
