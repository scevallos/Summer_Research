from Tkinter import *
import tkFileDialog as fd

def GUI():
	window = Tk()
	B = Button(window, text="Hi", command= print_hi).grid(columnspan=2)

	a = "x"
	b = "z"
	c = "a"
	d = 'p'

	if a == "" or b == "" or c == "" or d == "":
		print 'true'

	window.mainloop()
def print_hi():
	print 'hi'

def main():
	GUI()

if __name__=="__main__":
    main()