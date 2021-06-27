import wx
import eel
import re
import math
import os
import numpy as np
from decimal import Decimal
import matplotlib.pyplot as plt
from math import factorial
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
        if type == 1:
            app = wx.App(None)
            style = wx.FD_OPEN | wx.FD_FILE_MUST_EXIST
            dialog = wx.FileDialog(None, 'Open', wildcard=wildcard, style=style)
            if dialog.ShowModal() == wx.ID_OK:
                path = dialog.GetPath()
            else:
                path = None
            dialog.Destroy()
            self.path_DNA_1 = path
            with open(path, "r") as f:
                text = f.read()
                text = re.sub('\d', '', text)
                text = text.lower()
                if 'U' in text or 'u' in text:
                    text = text.replace('u', 't')
                text = re.sub('[^cagt]', '', text)
                text = re.sub('[^A-Za-z0-9]', '', text)
                text = text.replace(' ', '')
                text = text.replace('\n', '')
                return path, len(text)
        if type == 2:
            app = wx.App(None)
            style = wx.FD_OPEN | wx.FD_FILE_MUST_EXIST
            dialog = wx.FileDialog(None, 'Open', wildcard=wildcard, style=style)
            if dialog.ShowModal() == wx.ID_OK:
                path = dialog.GetPath()
            else:
                path = None
            dialog.Destroy()
            self.path_DNA_2 = path
            with open(path, "r") as f:
                text = f.read()
                text = re.sub('\d', '', text)
                text = text.lower()
                if 'U' in text or 'u' in text:
                    text = text.replace('u', 't')
                text = re.sub('[^cagt]', '', text)
                text = re.sub('[^A-Za-z0-9]', '', text)
                text = text.replace(' ', '')
                text = text.replace('\n', '')
                return path, len(text)
        if type == 3:
            app = wx.App(None)
            dialog = wx.DirDialog(None, "Choose input directory", "",
                    wx.DD_DEFAULT_STYLE | wx.DD_DIR_MUST_EXIST)
            if dialog.ShowModal() == wx.ID_OK:
                path = dialog.GetPath()
            else:
                path = None
            dialog.Destroy()
            self.path_res = path
            return path, 0

    @eel.expose
    def run(self, threshold=None):
        if self.path_DNA_1 and self.path_DNA_2 and not threshold:
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
        if self.path_DNA_1 and self.path_DNA_2 and self.path_res and threshold:
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
            DNA_1_list = []
            DNA_2_list = []
            similarity_list = []
            length_dna = min(len(self.DNA_1), len(self.DNA_2))
            max_length = max(len(self.DNA_1), len(self.DNA_2))
            temp_DNA_1 = np.array(list(self.DNA_1))
            temp_DNA_2 = np.array(list(self.DNA_2))
            for i in range(max_length - length_dna):
                if len(self.DNA_1) > len(self.DNA_2):
                    similarity = (np.mean(temp_DNA_2 == temp_DNA_1[i:length_dna + i]) * 100)
                else:
                    similarity = (np.mean(temp_DNA_1 == temp_DNA_2[i:length_dna + i]) * 100)
                if similarity > threshold:
                    if len(self.DNA_1) > len(self.DNA_2):
                        DNA_2_list.append(self.DNA_2)
                        DNA_1_list.append(self.DNA_1[i:length_dna + i + math.ceil(length_dna * 0.2)])
                    else:
                        DNA_1_list.append(self.DNA_1)
                        DNA_2_list.append(self.DNA_2[i:length_dna + i + math.ceil(length_dna * 0.2)])
                    similarity_list.append(similarity)
            res = []
            str_res = ''
            for i in range(len(DNA_1_list)):
                a = StripedSmithWaterman(DNA_1_list[i].upper())
                res.append(a(DNA_2_list[0].upper()).aligned_query_sequence)
            for i in range(len(res)):
                str_res = str_res + f'Процент схожести: {similarity_list[i]}\nDNA_1: {DNA_1_list[i].upper()}\nDNA_2: {DNA_2_list[i].upper()}\n---------\nRES: {res[i].replace("-","")}\n\n'
            if len(str_res) == 0:
                str_res = 'Нет подходящих данных.'
            if self.path_res:
                with open(self.path_res + '\\res_alignment.txt', "w") as f:
                    f.write(str_res)
                os.startfile(self.path_res + '\\res_alignment.txt')
            return True
        return 'Выберете файлы или папку для сохранения результата'

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

    def show_graphic(self):
        if self.path_DNA_1 and self.path_DNA_2:
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
            if len(self.DNA_1) < len(self.DNA_2):
                res = []
                for i in range(len(self.DNA_1)):
                    res.append(self.probability_of_success(1/4, len(self.DNA_1), i))
                if self.path_res:
                    str_res = 'Количество положительный исходов k из n испытаний | Вероятность P, %\n'
                    for i in range(len(self.DNA_1)):
                        str_res = str_res + f'{i} | {res[i]}\n'
                    with open(self.path_res + '\\res_bern.txt', "w") as f:
                        f.write(str_res)
                    os.startfile(self.path_res + '\\res_bern.txt')
                plt.figure(figsize=(8, 6), dpi=80)
                plt.plot([i for i in range(len(self.DNA_1))], res)
                plt.xlabel('Количество положительный исходов k из n испытаний')
                plt.ylabel('Вероятность P, %')
                if self.path_res:
                    plt.savefig(self.path_res + '\\res_bern.png')
                plt.show()
            else:
                res = []
                for i in range(len(self.DNA_2)):
                    res.append(self.probability_of_success(1/4, len(self.DNA_2), i))
                if self.path_res:
                    str_res = 'Количество положительный исходов k из n испытаний | Вероятность P, %\n'
                    for i in range(len(self.DNA_2)):
                        str_res = str_res + f'{i} | {res[i]}\n'
                    with open(self.path_res + '\\res_bern.txt', "w") as f:
                        f.write(str_res)
                    os.startfile(self.path_res + '\\res_bern.txt')
                plt.figure(figsize=(8, 6), dpi=80)
                plt.plot([i for i in range(len(self.DNA_2))], res)
                plt.xlabel('Количество положительный исходов k из n испытаний')
                plt.ylabel('Вероятность P, %')
                if self.path_res:
                    plt.savefig(self.path_res + '\\res_bern.png')
                plt.show()
        elif self.path_DNA_1 or self.path_DNA_2:
            if self.path_DNA_1:
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
                    res = []
                    for i in range(len(self.DNA_1)):
                        res.append(self.probability_of_success(1 / 4, len(self.DNA_1), i))
                    if self.path_res:
                        str_res = 'Количество положительный исходов k из n испытаний | Вероятность P, %\n'
                        for i in range(len(self.DNA_1)):
                            str_res = str_res + f'{i} | {res[i]}\n'
                        with open(self.path_res + '\\res_bern.txt', "w") as f:
                            f.write(str_res)
                        os.startfile(self.path_res + '\\res_bern.txt')
                    plt.figure(figsize=(8, 6), dpi=80)
                    plt.plot([i for i in range(len(self.DNA_1))], res)
                    plt.xlabel('Количество положительный исходов k из n испытаний')
                    plt.ylabel('Вероятность P, %')
                    if self.path_res:
                        plt.savefig(self.path_res + '\\res_bern.png')
                    plt.show()
            else:
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
                    res = []
                    for i in range(len(self.DNA_2)):
                        res.append(self.probability_of_success(1 / 4, len(self.DNA_2), i) * 100)
                    if self.path_res:
                        str_res = 'Количество положительный исходов k из n испытаний | Вероятность P, %\n'
                        for i in range(len(self.DNA_2)):
                            str_res = str_res + f'{i} | {res[i]}\n'
                        with open(self.path_res + '\\res_bern.txt', "w") as f:
                            f.write(str_res)
                        os.startfile(self.path_res + '\\res_bern.txt')
                    plt.figure(figsize=(8, 6), dpi=80)
                    plt.plot([i for i in range(len(self.DNA_2))], res)
                    plt.xlabel('Количество положительный исходов k из n испытаний')
                    plt.ylabel('Вероятность P, %')
                    if self.path_res:
                        plt.savefig(self.path_res + '\\res_bern.png')
                    plt.show()

    def sliding(self):
        if self.path_DNA_1 and self.path_DNA_2:
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
            res = []
            length_dna = min(len(self.DNA_1), len(self.DNA_2))
            max_length = max(len(self.DNA_1), len(self.DNA_2))
            temp_DNA_1 = np.array(list(self.DNA_1))
            temp_DNA_2 = np.array(list(self.DNA_2))
            for i in range(max_length - length_dna):
                if len(self.DNA_1) > len(self.DNA_2):
                    res.append(np.mean(temp_DNA_2 == temp_DNA_1[i:length_dna+i]) * 100)
                else:
                    res.append(np.mean(temp_DNA_1 == temp_DNA_2[i:length_dna+i]) * 100)
            plt.figure(figsize=(8, 6), dpi=80)
            plt.plot([i for i in range(len(res))], res)
            plt.xlabel('Номер позиции')
            plt.ylabel('Процент совпадения, %')
            if self.path_res:
                str_res = 'Номер позиции | Процент совпадения, %\n'
                for i in range(len(res)):
                    str_res = str_res + f'{i} | {res[i]}\n'
                with open(self.path_res + '\\res_sliding.txt', "w") as f:
                    f.write(str_res)
                os.startfile(self.path_res + '\\res_sliding.txt')
            if self.path_res:
                plt.savefig(self.path_res + '\\res_sliding.png')
            plt.show()


    @staticmethod
    def num_of_successes(n, k):
        return Decimal(factorial(n)) / Decimal((factorial(k) * factorial(n - k)))

    def probability_of_success(self, p, n, k):
        C_kn = self.num_of_successes(n, k)
        return C_kn * Decimal((p ** k)) * Decimal((1 - p)) ** Decimal((n - k))
