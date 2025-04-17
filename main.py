import sys
print(f"Using Python version: {sys.version}")

import tkinter as tk
from GUI.app import App

def run_app():
    root = tk.Tk()
    app = App(root)
    root.mainloop()

if __name__ == "__main__":
    run_app()
