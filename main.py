import re
import math
import os
import numpy as np
from decimal import Decimal
import matplotlib.pyplot as plt
from math import factorial
from skbio.alignment import StripedSmithWaterman, local_pairwise_align_ssw
from sklearn.cross_decomposition import CCA
from skbio.sequence import DNA
from scipy.stats import pearsonr
import skbio.io.format
import skbio.io
import skbio.io.format.blast6
import skbio.io.format.blast7
import skbio.io.format.clustal
import skbio.io.format.embl
import skbio.io.format.emptyfile
import skbio.io.format.fasta
import skbio.io.format.fastq
import skbio.io.format.genbank
import skbio.io.format.gff3
import skbio.io.format.lsmat
import skbio.io.format.newick
import skbio.io.format.ordination
import skbio.io.format.phylip
import skbio.io.format.qseq
import skbio.io.format.stockholm
import skbio.io.format.tests
import skbio.io.format._base
import skbio.io.format._blast
import skbio.io.format._sequence_feature_vocabulary
import swalign
from scipy.spatial.distance import cosine, cityblock, minkowski

from pandas import DataFrame, ExcelWriter

""" CONST """
PATH = os.path.dirname(os.path.abspath(__file__))


class Program:
    DNA_1 = None
    path_DNA_1 = os.path.join(PATH, 'input', 'DNA_1.txt')
    path_DNA_2 = os.path.join(PATH, 'input', 'DNA_2.txt')
    path_res = PATH
    DNA_2 = None
    cca = CCA(n_components=1)
    sw = swalign.LocalAlignment(swalign.NucleotideScoringMatrix(2, -1))

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
                cb = round((cb_real + cb_imag) / 2, 5)
                mk_real = minkowski(list(map(lambda x: x.real, temp_DNA_1)), list(map(lambda x: x.real, temp_DNA_2)))
                mk_imag = minkowski(list(map(lambda x: x.imag, temp_DNA_1)), list(map(lambda x: x.imag, temp_DNA_2)))
                mk = round((mk_real + mk_imag) / 2, 5)
                return 'Одна из ДНК неизменяема. Очевидно корреляция - 0, коэффициент Отиаи - 1. Выведем другие данные', cb, mk
            p_real = pearsonr(list(map(lambda x: x.real, temp_DNA_1)), list(map(lambda x: x.real, temp_DNA_2)))
            p_imag = pearsonr(list(map(lambda x: x.imag, temp_DNA_1)), list(map(lambda x: x.imag, temp_DNA_2)))
            p = round((p_real[0] + p_imag[0]) / 2, 5)
            if p < 0.0:
                p = round(abs(p), 5)
            p_real = cosine(list(map(lambda x: x.real, temp_DNA_1)), list(map(lambda x: x.real, temp_DNA_2)))
            p_imag = cosine(list(map(lambda x: x.imag, temp_DNA_1)), list(map(lambda x: x.imag, temp_DNA_2)))
            c = round((p_real + p_imag) / 2, 5)
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
            similarity_align = []
            divergence = []
            for i in range(len(DNA_1_list)):
                try:
                    res.append(local_pairwise_align_ssw(DNA(DNA_1_list[i].upper()), DNA(DNA_2_list[i].upper()),
                                                        gap_open_penalty=1, gap_extend_penalty=1)[0].iloc[1])
                except:
                    a = StripedSmithWaterman(DNA_1_list[i].upper())
                    res.append(a(DNA_2_list[i].upper()).aligned_query_sequence)
                if len(self.DNA_1) < len(self.DNA_2):
                    temp_DNA_1 = np.array(list(DNA_1_list[i].upper().replace('-', '')))
                    temp_DNA_2 = np.array(list(str(res[i]).replace('-', '')))
                    temp = str(res[i]).count('-') + DNA_1_list[i].count('-')
                    if len(temp_DNA_1) < len(temp_DNA_2):
                        for j in range(len(temp_DNA_1)):
                            if not temp_DNA_1[j] == temp_DNA_2[j]:
                                temp = temp + 1
                        similarity_align.append(100 - (temp / len(temp_DNA_1)) * 100)
                    else:
                        for j in range(len(temp_DNA_2)):
                            if not temp_DNA_1[j] == temp_DNA_2[j]:
                                temp = temp + 1
                        similarity_align.append(100 - (temp / len(temp_DNA_2)) * 100)
                else:
                    temp_DNA_1 = np.array(list(DNA_2_list[i].upper().replace('-', '')))
                    temp_DNA_2 = np.array(list(str(res[i]).replace('-', '')))
                    temp = str(res[i]).count('-') + DNA_2_list[i].count('-')
                    if len(temp_DNA_1) < len(temp_DNA_2):
                        for j in range(len(temp_DNA_1)):
                            if not temp_DNA_1[j] == temp_DNA_2[j]:
                                temp = temp + 1
                        similarity_align.append(100 - (temp / len(temp_DNA_1)) * 100)
                    else:
                        for j in range(len(temp_DNA_2)):
                            if not temp_DNA_1[j] == temp_DNA_2[j]:
                                temp = temp + 1
                        similarity_align.append(100 - (temp / len(temp_DNA_2)) * 100)
                if len(self.DNA_1) < len(self.DNA_2):
                    divergence.append(round(100 - (len(DNA_1_list[i]) / len(res[i])) * 100, 3))
                else:
                    divergence.append(round(100 - (len(DNA_2_list[i]) / len(res[i])) * 100, 3))
            str_res = ''
            for i in range(len(res)):
                if self.path_res:
                    str_res = str_res + f'Процент схожести: {similarity_align[i]} \nПроцент растяжения: {divergence[i]} \nDNA_1: {DNA_1_list[i].upper()} \nDNA_2: {DNA_2_list[i].upper()} \nRES       : {res[i]}\n\n'
            if self.path_res:
                self.make_xlsx(str_res, 'res_align')
                with open(os.path.join(PATH, 'output', 'stats_alignment.txt'), "w") as f:
                    sum_DNA_len = 0
                    sum_res_len = 0
                    sum_hole_len = 0
                    if len(self.DNA_1) > len(self.DNA_2):
                        for i in range(len(DNA_2_list)):
                            sum_DNA_len = sum_DNA_len + len(DNA_2_list[i])
                    else:
                        for i in range(len(DNA_1_list)):
                            sum_DNA_len = sum_DNA_len + len(DNA_1_list[i])
                    for i in range(len(res)):
                        sum_res_len = sum_res_len + len(res[i])
                    if len(self.DNA_1) > len(self.DNA_2):
                        for i in range(len(DNA_2_list)):
                            sum_hole_len = abs(sum_hole_len + (len(DNA_2_list[i]) - len(res[i])))
                    else:
                        for i in range(len(DNA_1_list)):
                            sum_hole_len = abs(sum_hole_len + (len(DNA_1_list[i]) - len(res[i])))
                    if len(res) == 0:
                        f.write('Нет подходящих данных.')
                    else:
                        f.write(
                            f'При пороге = {threshold}%\n\nОбщая длина ДНК: {sum_DNA_len}\nОбщая длина сходимостей при выравнивании: {sum_res_len}\nОбщая длина пропусков при выравнивании: '
                            f'{sum_hole_len}\n\nПроцент сходимости: {(sum_res_len / sum_DNA_len) * 100}\nПроцент пропусков: {(sum_hole_len / sum_DNA_len) * 100}\n\n'
                            f'Средняя длина ДНК: {sum_DNA_len / len(res)}\nСреднее число сходимостей: {sum_res_len / len(res)}\nСреднее число пропусков: {sum_hole_len / len(res)}')
                os.startfile(os.path.join(PATH, 'output', 'stats_alignment.txt'))
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

    # xlsx
    def make_xlsx(self, data, name):
        if not self.path_res:
            return
        data = data.split('\n')
        result = []
        for i in data:
            if not i: continue
            line = i.split('|')
            result.append(line)

        df = DataFrame(result)
        path = os.path.join(self.path_res, 'output', f'{name}.xlsx')
        writer = ExcelWriter(path, engine='xlsxwriter')
        df.to_excel(writer, sheet_name='result', index=False)
        writer.save()
        return

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
                    res.append(self.probability_of_success(1 / 4, len(self.DNA_1), i))
                if self.path_res:
                    str_res = 'Количество положительный исходов k из n испытаний | Вероятность P, %\n'
                    for i in range(len(self.DNA_1)):
                        str_res = str_res + f'{i} | {res[i]}\n'
                    try:
                        self.make_xlsx(str_res, 'res_bern')
                    except:
                        with open(os.path.join(self.path_res, 'output', 'res_bern.txt'), "w") as f:
                            f.write(str_res)
                        os.startfile(os.path.join(self.path_res, 'output', 'res_bern.txt'))
                plt.figure(figsize=(8, 6), dpi=80)
                plt.plot([i for i in range(len(self.DNA_1))], res)
                plt.xlabel('Количество положительный исходов k из n испытаний')
                plt.ylabel('Вероятность P, %')
                if self.path_res:
                    plt.savefig(os.path.join(self.path_res, 'output', 'res_bern.png'))
                plt.show()
            else:
                res = []
                for i in range(len(self.DNA_2)):
                    res.append(self.probability_of_success(1 / 4, len(self.DNA_2), i))
                if self.path_res:
                    str_res = 'Количество положительный исходов k из n испытаний | Вероятность P, %\n'
                    for i in range(len(self.DNA_2)):
                        str_res = str_res + f'{i} | {res[i]}\n'
                    try:
                        self.make_xlsx(str_res, 'res_bern')
                    except:
                        with open(os.path.join(self.path_res, 'output', 'res_bern.txt'), "w") as f:
                            f.write(str_res)
                        os.startfile(os.path.join(self.path_res, 'output', 'res_bern.txt'))
                plt.figure(figsize=(8, 6), dpi=80)
                plt.plot([i for i in range(len(self.DNA_2))], res)
                plt.xlabel('Количество положительный исходов k из n испытаний')
                plt.ylabel('Вероятность P, %')
                if self.path_res:
                    plt.savefig(os.path.join(self.path_res, 'output', 'res_bern.png'))
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
                            try:
                                self.make_xlsx(str_res, 'res_bern')
                            except:
                                with open(os.path.join(self.path_res, 'output', 'res_bern.txt'), "w") as f:
                                    f.write(str_res)
                                os.startfile(os.path.join(self.path_res, 'output', 'res_bern.txt'))
                    plt.figure(figsize=(8, 6), dpi=80)
                    plt.plot([i for i in range(len(self.DNA_1))], res)
                    plt.xlabel('Количество положительный исходов k из n испытаний')
                    plt.ylabel('Вероятность P, %')
                    if self.path_res:
                        plt.savefig(self.path_res + os.path.join(self.path_res, 'output', 'res_bern.png'))
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
                            try:
                                self.make_xlsx(str_res, 'res_bern')
                            except:
                                with open(os.path.join(self.path_res, 'output', 'res_bern.txt'), "w") as f:
                                    f.write(str_res)
                                os.startfile(os.path.join(self.path_res, 'output', 'res_bern.txt'))
                    plt.figure(figsize=(8, 6), dpi=80)
                    plt.plot([i for i in range(len(self.DNA_2))], res)
                    plt.xlabel('Количество положительный исходов k из n испытаний')
                    plt.ylabel('Вероятность P, %')
                    if self.path_res:
                        plt.savefig(os.path.join(self.path_res, 'output', 'res_bern.png'))
                    plt.show()

    def sliding(self, threshold=None):
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
            if threshold:
                threshold_res = []
                threshold_index = []
                threshold_coordiantes = []
            coordinates = []
            length_dna = min(len(self.DNA_1), len(self.DNA_2))
            max_length = max(len(self.DNA_1), len(self.DNA_2))
            temp_DNA_1 = np.array(list(self.DNA_1))
            temp_DNA_2 = np.array(list(self.DNA_2))
            for i in range(max_length - length_dna):
                if len(self.DNA_1) > len(self.DNA_2):
                    coordinates_temp = []
                    for j in range(len(self.DNA_2)):
                        if temp_DNA_2[j] == temp_DNA_1[i:length_dna + i][j]:
                            coordinates_temp.append(f'({j}, {i + j}, {temp_DNA_2[j].upper()})')
                    coordinates.append(coordinates_temp)
                    mean = np.mean(temp_DNA_2 == temp_DNA_1[i:length_dna + i]) * 100
                    res.append(mean)
                    if threshold:
                        if mean > threshold:
                            threshold_res.append(mean)
                            threshold_index.append(i)
                            threshold_coordiantes.append(coordinates_temp)
                else:
                    coordinates_temp = []
                    for j in range(len(self.DNA_1)):
                        if temp_DNA_1[j] == temp_DNA_2[i:length_dna + i][j]:
                            coordinates_temp.append(f'({j}, {j - i}, {temp_DNA_2[j].upper()})')
                    coordinates.append(coordinates_temp)
                    mean = np.mean(temp_DNA_1 == temp_DNA_2[i:length_dna + i]) * 100
                    res.append(mean)
                    if threshold:
                        if mean > threshold:
                            threshold_res.append(mean)
                            threshold_index.append(i)
                            threshold_coordiantes.append(coordinates_temp)
            plt.figure(figsize=(8, 6), dpi=80)
            plt.plot([i for i in range(len(res))], res)
            plt.xlabel('Номер позиции')
            plt.ylabel('Процент совпадения, %')
            if self.path_res:
                str_res = 'Номер позиции | Процент совпадения, % | Координаты совпадений\n'
                for i in range(len(res)):
                    str_res = str_res + f'{i} | {res[i]} | {coordinates[i]}\n'
                try:
                    self.make_xlsx(str_res, 'res_sliding')
                except:
                    with open(os.path.join(self.path_res, 'output', 'res_sliding.txt'), "w") as f:
                        f.write(str_res)
                    os.startfile(os.path.join(self.path_res, 'output', 'res_sliding.txt'))
            if self.path_res and threshold:
                str_res = f'Номер позиции | Процент совпадения, % | Координаты совпадений\nПорог = {threshold}%\n'
                for i in range(len(threshold_res)):
                    str_res = str_res + f'{threshold_index[i]} | {threshold_res[i]} | {threshold_coordiantes[i]}\n\n'
                try:
                    self.make_xlsx(str_res,'res_sliding_with_threshold')
                except:
                    with open(os.path.join(self.path_res, 'output', 'res_sliding_with_threshold.txt'), "w") as f:
                        f.write(str_res)
                    os.startfile(os.path.join(self.path_res, 'output', 'res_sliding_with_threshold.txt'))
            if self.path_res:
                plt.savefig(os.path.join(self.path_res, 'output', 'res_sliding_with_threshold.png'))
            plt.show()

    @staticmethod
    def num_of_successes(n, k):
        return Decimal(factorial(n)) / Decimal((factorial(k) * factorial(n - k)))

    def probability_of_success(self, p, n, k):
        C_kn = self.num_of_successes(n, k)
        return C_kn * Decimal((p ** k)) * Decimal((1 - p)) ** Decimal((n - k))


if __name__ == '__main__':
    A = Program()
    print('Введите вид обработки данных')
    print('1 - Найти коэффициенты последовательностей (Коэффициент Отиаи, Корреляционный коэффициент')
    print('2 - Открыть и сохранить результаты формулы Бернулли')
    print('3 - Открыть и сохранить результаты скользящего сравнения')
    print('4 - Открыть и сохранить результаты выравнивания')
    print('\n')
    while True:
        type = input('>>')
        if type == '1':
            res = A.run()
            print(f'Корреляционный коэффициент = {res[0]} | Коэффициент Отиаи = {res[1]}')
        if type == '2':
            A.show_graphic()
        if type == '3':
            print('Введите пороговое значение (от 0 до 100). Если оно отсутствует, введите 0')
            tr = input('')
            if tr == '0':
                A.sliding()
            else:
                tr = int(tr)
                A.sliding(tr)
        if type == '4':
            print('Введите пороговое значение (от 0 до 100). Если оно отсутствуе, введите 0')
            tr = input('')
            tr = int(tr)
            res = A.run(tr)
