

# import numpy as np
# import rcm


# a = rcm.B(0,1,2,3,0,7,8,6,8,9,15,5)
# print(a)


import tkinter as tk
from tkinter import filedialog

root = tk.Tk()
root.withdraw()

# file_path = filedialog.askopenfilename(
#   title="Select geometry xyz file.",
#   initialdir='.',
#   filetypes=(
#       # ("jpeg files", "*.jpg"),
#       # ("csv files", "*.csv"),
#       ("xyz files", "*.xyz"),
#       ("gif files", "*.gif*")
#       # ("png files", "*.png")
#     )
#   )

file_path = 'C:/Users/edu51/MJ_Python/rcm-python/cp6t6_0p_blyp35_geom.xyz'

from rdkit import Chem
a = Chem.rdmolfiles.MolFromXYZFile(file_path)


print(file_path)