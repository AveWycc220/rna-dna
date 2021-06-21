import wx
import eel
import re
import math
import os
from skbio.alignment import StripedSmithWaterman, local_pairwise_align
from sklearn.cross_decomposition import CCA
from scipy.stats import pearsonr
import swalign
from scipy.spatial.distance import cosine, cityblock, minkowski

""" CONST """
PATH = os.path.dirname(os.path.abspath(__file__)) + '\\'


class Program:
    DNA_1 = None
    path_DNA_1 = None
    path_DNA_2 = None
    path_res = None
    DNA_2 = None
    cca = CCA(n_components=1)
    sw = swalign.LocalAlignment(swalign.NucleotideScoringMatrix(2, -1))

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
        if type == 1:
            self.path_DNA_1 = path
        if type == 2:
            self.path_DNA_2 = path
        if type == 3:
            self.path_res = path
        return path

    @eel.expose
    def run(self, length=None, threshold=None):
        if self.path_DNA_1 and self.path_DNA_2 and not length and not threshold:
            with open(self.path_DNA_1, "r") as f:
                text = f.read()
                text = re.sub('\d', '', text)
                text = text.lower()
                if 'U' in text or 'u' in text:
                    text = text.replace('u', 't')
                text = re.sub('[^cagt]', '', text)
                text = re.sub('[^A-Za-z0-9]', '', text)
                text = text.replace(' ', '')
                text = text.replace('\n', '')
                self.DNA_1 = text
            with open(self.path_DNA_2, "r") as f:
                text = f.read()
                text = re.sub('\d', '', text)
                text = text.lower()
                if 'U' in text or 'u' in text:
                    text = text.replace('u', 't')
                text = re.sub('[^cagt]', '', text)
                text = re.sub('[^A-Za-z0-9]', '', text)
                text = text.replace(' ', '')
                text = text.replace('\n', '')
                self.DNA_2 = text
            length_dna = min(len(self.DNA_1), len(self.DNA_2))
            temp_DNA_1 = self._convert(self.DNA_1)[:length_dna]
            temp_DNA_2 = self._convert(self.DNA_2)[:length_dna]
            if temp_DNA_1.count(temp_DNA_1[0]) == len(temp_DNA_1) or temp_DNA_2.count(temp_DNA_2[0]) == len(temp_DNA_2):
                cb_real = cityblock(list(map(lambda x: x.real, temp_DNA_1)), list(map(lambda x: x.real, temp_DNA_2)))
                cb_imag = cityblock(list(map(lambda x: x.imag, temp_DNA_1)), list(map(lambda x: x.imag, temp_DNA_2)))
                cb = round((cb_real + cb_imag)/2, 5)
                mk_real = minkowski(list(map(lambda x: x.real, temp_DNA_1)), list(map(lambda x: x.real, temp_DNA_2)))
                mk_imag = minkowski(list(map(lambda x: x.imag, temp_DNA_1)), list(map(lambda x: x.imag, temp_DNA_2)))
                mk = round((mk_real + mk_imag)/2, 5)
                return 'Одна из ДНК неизменяема. Очевидно корреляция - 0, коэффициент Отиаи - 1. Выведем другие данные', cb, mk
            p_real = pearsonr(list(map(lambda x: x.real, temp_DNA_1)), list(map(lambda x: x.real, temp_DNA_2)))
            p_imag = pearsonr(list(map(lambda x: x.imag, temp_DNA_1)), list(map(lambda x: x.imag, temp_DNA_2)))
            p = round((p_real[0] + p_imag[0])/2, 5)
            if p < 0.0:
                p = round(abs(p), 5)
            p_real = cosine(list(map(lambda x: x.real, temp_DNA_1)), list(map(lambda x: x.real, temp_DNA_2)))
            p_imag = cosine(list(map(lambda x: x.imag, temp_DNA_1)), list(map(lambda x: x.imag, temp_DNA_2)))
            c = round((p_real + p_imag)/2, 5)
            if c > 1.0:
                c = round(1 - (c - 1), 5)
            return p, c
        if self.path_DNA_1 and self.path_DNA_2 and self.path_res and length and threshold:
            with open(self.path_DNA_1, "r") as f:
                text = f.read()
                text = re.sub('\d', '', text)
                text = text.lower()
                if 'U' in text or 'u' in text:
                    text = text.replace('u', 't')
                text = re.sub('[^cagt]', '', text)
                text = re.sub('[^A-Za-z0-9]', '', text)
                text = text.replace(' ', '')
                text = text.replace('\n', '')
                self.DNA_1 = text
            with open(self.path_DNA_2, "r") as f:
                text = f.read()
                text = re.sub('\d', '', text)
                text = text.lower()
                if 'U' in text or 'u' in text:
                    text = text.replace('u', 't')
                text = re.sub('[^cagt]', '', text)
                text = re.sub('[^A-Za-z0-9]', '', text)
                text = text.replace(' ', '')
                text = text.replace('\n', '')
                self.DNA_2 = text
            length_dna = min(len(self.DNA_1), len(self.DNA_2))
            temp_DNA_1 = self._convert(self.DNA_1)[:length_dna]
            temp_DNA_2 = self._convert(self.DNA_2)[:length_dna]
            if length > length_dna:
                return 'Заданная длина больше, чем длина ДНК'
            DNA_1_list = []
            DNA_2_list = []
            correlation_list = []
            for i in range(len(temp_DNA_1)-length):
                if len(set(temp_DNA_1[i:i + length])) != 1 and len(set(temp_DNA_2[i:i + length])) != 1:
                    window_p_real = pearsonr(list(map(lambda x: x.real, temp_DNA_1[i:i + length])),
                                             list(map(lambda x: x.real, temp_DNA_2[i:i + length])))
                    window_p_imag = pearsonr(list(map(lambda x: x.imag, temp_DNA_1[i:i + length])),
                                             list(map(lambda x: x.imag, temp_DNA_2[i:i + length])))
                    window_p = round((window_p_real[0] + window_p_imag[0]) / 2, 5)
                    if window_p < 0.0:
                        window_p = round(abs(window_p), 5)
                    if window_p > threshold:
                        DNA_1_list.append(self.DNA_1[i:i + length])
                        DNA_2_list.append(self.DNA_2[i:i + length + math.ceil(length * 0.15)])
                        correlation_list.append(window_p)
            res = []
            str_res = ''
            for i in range(len(DNA_1_list)):
                a = StripedSmithWaterman(DNA_1_list[i].upper())
                res.append(a(DNA_2_list[0].upper()).aligned_query_sequence)
                print(a(DNA_2_list[0].upper()).aligned_query_sequence)
            for i in range(len(res)):
                str_res = str_res + f'Корреляционный коэффициент: {correlation_list[i]}\nDNA_1: {DNA_1_list[i].upper()}\nDNA_2: {DNA_2_list[i].upper()}\n---------\nRES: {res[i].replace("-","")}\n\n'
            if len(str_res) == 0:
                str_res = 'Нет подходящих данных.'
            with open(self.path_res, "w") as f:
                f.write(str_res)
            os.startfile(self.path_res)
            return True
        return 'Выберете файлы'

    @staticmethod
    def _convert(l):
        def convert_letters(x):
            if x == 'a':
                return 0 + 1j
            elif x == 'c':
                return 1 + 0j
            elif x == 'g':
                return -1 + 0j
            elif x == 't':
                return 0 - 1j
        return list(map(convert_letters, l))

    def open_file(self):
        os.startfile(self.path_res)
