import wx
import eel
import re
from sklearn.cross_decomposition import CCA
from scipy.stats import pearsonr
from scipy.spatial.distance import cosine


class Program:
    DNA_1 = None
    DNA_2 = None
    cca = CCA(n_components=1)

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
                text = text.lower()
                if 'U' in text or 'u' in text:
                    text = text.replace('u', 't')
                text = re.sub('^[CAGTcagt]+$', '', text)
                text = re.sub('[^A-Za-z0-9]', '', text)
                text = text.replace(' ', '')
                text = text.replace('\n', '')
                if type == 1:
                    self.DNA_1 = text
                if type == 2:
                    self.DNA_2 = text
        return path

    @eel.expose
    def run(self):
        if self.DNA_1 and self.DNA_2:
            length = min(len(self.DNA_1), len(self.DNA_2))
            temp_DNA_1 = self._convert(self.DNA_1)[:length]
            temp_DNA_2 = self._convert(self.DNA_2)[:length]
            p_real = pearsonr(list(map(lambda x: x.real, temp_DNA_1)), list(map(lambda x: x.real, temp_DNA_2)))
            p_imag = pearsonr(list(map(lambda x: x.imag, temp_DNA_1)), list(map(lambda x: x.imag, temp_DNA_2)))
            p = (p_real[0] + p_imag[0])/2
            if p < 0.0:
                p = abs(p)
            p_real = cosine(list(map(lambda x: x.real, temp_DNA_1)), list(map(lambda x: x.real, temp_DNA_2)))
            p_imag = cosine(list(map(lambda x: x.imag, temp_DNA_1)), list(map(lambda x: x.imag, temp_DNA_2)))
            c = (p_real + p_imag)/2
            if c > 1.0:
                c = 1 - (c - 1)
            return p, c
        return 'Выберете файлы'

    @staticmethod
    def _convert(l):
        def convert_letters(x):
            if x == 'a':
                return 0 + 1j
            if x == 'c':
                return 1 + 0j
            if x == 'g':
                return -1 + 0j
            if x == 't':
                return 0 - 1j
        return list(map(convert_letters, l))

