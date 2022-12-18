

# import numpy as np
# import rcm


# a = rcm.B(0,1,2,3,0,7,8,6,8,9,15,5)
# print(a)


import tkinter as tk
from tkinter import filedialog

root = tk.Tk()
root.withdraw()

file_path = filedialog.askopenfilename(
  title="lel",
  filetypes=(
      ("jpeg files", "*.jpg"),
      ("csv files", "*.csv"),
      ("gif files", "*.gif*"),
      ("png files", "*.png")
    )
  )





print(file_path)