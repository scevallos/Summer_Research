from Tkinter import *
from ttk import *
import itertools
import tkFileDialog

def askdirectory():
  dirname = tkFileDialog.askdirectory()
  if dirname:
    var.set(dirname)

def UserFileInput(status,name):
  optionFrame = Frame(root)
  optionLabel = Label(optionFrame)
  optionLabel["text"] = name
  optionLabel.pack(side=LEFT)
  text = status
  var = StringVar(root)
  var.set(text)
  w = Entry(optionFrame, textvariable= var)
  w.pack(side = LEFT)
  optionFrame.pack()
  return w, var

def Print_entry():
  print var.get()

if __name__ == '__main__':
  root = Tk()

  dirBut = Button(root, text='askdirectory', command = askdirectory)
  dirBut.pack(side = RIGHT)
  getBut = Button(root, text='print entry text', command = Print_entry)
  getBut.pack(side = BOTTOM)

  w, var = UserFileInput("", "Directory")

  root.mainloop()