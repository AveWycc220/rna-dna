import wx
import eel
import numpy as np
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.cross_decomposition import CCA
from sklearn.preprocessing import LabelEncoder
from scipy.spatial import distance


class Program:
    DNA_1 = None
    DNA_2 = None
    RNA_1 = None
    RNA_2 = None
    cca = CCA(n_components=1)
    le = LabelEncoder()

    @eel.expose
    def get_data(self, type, wildcard="*"):
        app = wx.App(None)
        style = wx.FD_OPEN | wx.FD_FILE_MUST_EXIST
        dialog = wx.FileDialog(None, 'Open', wildcard=wildcard, style=style)
        if dialog.ShowModal() == wx.ID_OK:
            path = dialog.GetPath()
        else:
            path = None
        dialog.Destroy()
        if path:
            with open(path, "r") as f:
                text = f.read()
                text = re.sub('\d', '', text)
                text = re.sub('^[CAGTcagt]+$', '', text)
                text = re.sub('[^A-Za-z0-9]', '', text)
                text = text.replace(' ', '')
                text = text.replace('\n', '')
                if type == 1:
                    self.DNA_1 = text
                    text = text.replace('T', 'U')
                    text = text.replace('t', 'u')
                    self.RNA_1 = text
                if type == 2:
                    self.DNA_2 = text
                    text = text.replace('T', 'U')
                    text = text.replace('t', 'u')
                    self.RNA_2 = text
        return path

    @eel.expose
    def run(self):
        if self.DNA_1 and self.DNA_2:
            temp_DNA_1 = self._convert_dna(self.DNA_1)
            temp_DNA_2 = self._convert_dna(self.DNA_2)
            shape = min(temp_DNA_1.shape[0], temp_DNA_2.shape[0])
            DNA_1_C, DNA_2_C = self.cca.fit_transform(temp_DNA_1[:shape].reshape(-1, 1), temp_DNA_2[:shape].reshape(-1, 1))
            cor = np.corrcoef(DNA_1_C, DNA_2_C, rowvar=False).diagonal(1)[0]
            return cor, distance.cosine(temp_DNA_1[:shape], temp_DNA_2[:shape])
        return 'Выберете файлы'

    def _convert_dna(self, l):
        return self.le.fit_transform([*l])

