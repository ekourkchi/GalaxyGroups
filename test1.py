#Import Tkinter
from Tkinter import *

#Main Frame
class Application(Frame):
    def __init__(self, master):  #initialize the grid and widgets
        Frame.__init__(self,master)
        self.grid()
        self.redFUN() #initialize the red frame's Function
        self.greenFUN() #initialize the green frame's Function
        self.widgets() #To show that you can still place non-Frame widgets 
    def widgets(self):
        self.mylabel = Label (self, text = "Hello World!")
        self.mylabel.grid()
    def redFUN(self): #The 'self' means that it is an instance of the main frame
        #Init the red frame
        self.redFrame = Frame(root, width = 100, height = 50,pady = 5,
                              bg = "red")
        self.redFrame.grid()



    def greenFUN(self): #Child of the mainframe
        self.greenFrame = Frame(root, width = 100, height = 50,pady = 5,
                          bg = "green") #it is green!
        self.greenFrame.grid()








#These lines of code are used for the grid
root = Tk()
root.title("Frame Example")
root.geometry("300x300")
app = Application(root)

root.mainloop()