# -*- coding: utf-8 -*-
#!/usr/bin/python
import numpy as np
from collections import defaultdict
import sys
import argparse

parser = argparse.ArgumentParser(
    description=f"Judgment of belonging to the same species was made based on the value of ANI. \
    The criterion is Similarity < 95% and Similarity < min_similarity (mean-SD), we think that are no the same Species.\
    For example: python {sys.argv[0]} -i input.fasta -o output")
parser.add_argument("-i",type=str,
    help="the result from FastANI")
parser.add_argument("-o",type=str,
    help="the filter result which may represent different species")
args = parser.parse_args()


def load_data(file_path):
    # 从文件加载数据
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            values = line.strip().split()
            data.append((values[0], values[1], float(values[2]), int(values[3]), int(values[4])))
    return data
      
def extract_species_name(file_name):
    # 提取物种名
    species_name = file_name.split("-")[0]
    return species_name

def filter_data(data):
    # 按物种分组
    species_data = defaultdict(list)
    for row in data:
        file1, file2, similarity, _, _ = row
        species1 = extract_species_name(file1)
        species2 = extract_species_name(file2)
        if species1 == species2:  # 判断物种名A是否等于物种名B
            species_data[species1].append(similarity)

    # 计算每个物种的相似度平均值和标准误差
    species_mean_similarity = {}
    species_sd_similarity = {}
    for species, similarities in species_data.items():
        mean_similarity = np.mean(similarities)
        sd_similarity = np.std(similarities)
        species_mean_similarity[species] = mean_similarity
        species_sd_similarity[species] = sd_similarity

    # 根据筛选条件进行数据筛选
    filtered_data = []
    for row in data:
        file1, file2, similarity, _, _ = row
        species1 = extract_species_name(file1)
        species2 = extract_species_name(file2)
        contig_name = file2.strip().split("-")[1]
        if species1 == species2:  # 判断物种名A是否等于物种名B
            min_similarity = species_mean_similarity[species1] - species_sd_similarity[species1]
            if similarity < 95 and similarity < min_similarity:
                filtered_data.append((species2,contig_name,file1,similarity,min_similarity,file2))

    # 查找种子的最高相似性值并添加到结果中
    seed_max_similarity = {}
    # processed_seeds = set()
    for row in data:
        seed_file1, seed, similarity, _, _= row
        if seed_file1 == seed:
            # 跳过相同种子的比较
            continue
        if seed not in seed_max_similarity or similarity > seed_max_similarity[seed][0]:
            seed_max_similarity[seed] = (similarity,seed_file1)

    # 将种子的最高相似性值添加到结果中
    result_data = []
    for row in filtered_data:
        species2, contig_name, file1, similarity, min_similarity, file2 = row
        seed_max_ani,seed_max_file1 = seed_max_similarity[file2]
        if species2 == seed_max_file1.strip().split("-")[0]:
            result = "match"
        else:
            result = "no match"
        result_data.append((species2, contig_name, file1, similarity, min_similarity, seed_max_file1,seed_max_ani,result))

    return result_data
    # return filtered_data

if __name__ == "__main__":
    file_path = args.i #"Dietzia.file.list.ani.txt" 
    data = load_data(file_path)
    result_data = filter_data(data)

    with open(args.o, 'w') as file:
        file.write("Bact_name\tContig\tfile1\tANI\tMin_ANI\tBest_match\tBest_match_ANI\tresult\n")
        for row in result_data:
            file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(row[0], row[1], row[2],row[3],row[4],row[5],row[6],row[7]))