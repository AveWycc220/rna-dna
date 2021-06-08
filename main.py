import eel
import os
from middle.controler import *

""" CONST """
PATH = os.path.dirname(os.path.abspath(__file__)) + '\\'

if __name__ == '__main__':
    eel.init(PATH)
    eel.start(PATH + 'front\\index.html', mode='chrome', port=0, size=(1000, 720))
