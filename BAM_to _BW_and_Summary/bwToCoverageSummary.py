import pandas as pd
import numpy as np
import pyBigWig
import argparse
import re
from datetime import datetime


class BigWigTOCoverageSum:
    bw_file_list = []
    chrom_list = []
    start_pos_list = []
    end_pos_list = []
    total_positions = 0

    def natural_key(self, string_):
        return [int(s) if s.isdigit() else s for s in re.split(r"(\d+)", string_)]

    def parse_bw_file_list(self, bw_file_name):
        with open(bw_file_name) as fileHandler:
            line = fileHandler.readline()
            while line:
                self.bw_file_list.append(line)
                line = fileHandler.readline()

        print("Number of bigwig files: " + str(len(self.bw_file_list)))

    def parse_bed_file(self, bed_file_name):
        with open(bed_file_name) as fileHandler:
            line = fileHandler.readline()
            while line:
                self.chrom_list.append(line.split("\t")[0])
                self.start_pos_list.append(int(line.split("\t")[1]))
                self.end_pos_list.append(int(line.split("\t")[2]))
                line = fileHandler.readline()
        print("Number of lines in bed file: " + str(len(self.chrom_list)))

        for i in range(len(self.start_pos_list)):
            self.total_positions += self.end_pos_list[i] - self.start_pos_list[i]

        print("Total positions: " + str(self.total_positions))

    def create_coverage_summary(self, output_file):
        chrom_set = set(self.chrom_list)
        chrom_set = sorted(chrom_set, key=self.natural_key)
        print(chrom_set)
        for chrom in chrom_set:
            self.create_coverage_summary_per_chr(chrom, output_file)

        print("Finished generating coverage summary statistics.")

    def create_coverage_summary_per_chr(self, chrom, output_file):
        mean_dict = {}
        over_1_dict = {}
        over_5_dict = {}
        over_10_dict = {}
        over_15_dict = {}
        over_20_dict = {}
        over_25_dict = {}
        over_30_dict = {}
        over_50_dict = {}
        over_100_dict = {}

        for bw_file in self.bw_file_list:
            bw_file = bw_file.strip()
            bw = pyBigWig.open(bw_file)
            for i in range(len(self.chrom_list)):
                if self.chrom_list[i] == chrom:
                    #print("Start reading coverage values from bw ")
                    #print(datetime.now())
                    coverage_values = bw.values(
                        self.chrom_list[i], self.start_pos_list[i], self.end_pos_list[i], numpy=True
                    )
                    coverage_values[np.isnan(coverage_values)] = 0
                    #print("Finish reading coverage values from bw ")
                    #print(datetime.now())

                    over_1_list = np.where(coverage_values >= 1, 1, 0)
                    over_5_list = np.where(coverage_values >= 5, 1, 0)
                    over_10_list = np.where(coverage_values >= 10, 1, 0)
                    over_15_list = np.where(coverage_values >= 15, 1, 0)
                    over_20_list = np.where(coverage_values >= 20, 1, 0)
                    over_25_list = np.where(coverage_values >= 25, 1, 0)
                    over_30_list = np.where(coverage_values >= 30, 1, 0)
                    over_50_list = np.where(coverage_values >= 50, 1, 0)
                    over_100_list = np.where(coverage_values >= 100, 1, 0)
                    #print("Finish creating numpy arr for over_1, over_5 ... values ")
                    #print(datetime.now())

                    if self.start_pos_list[i] not in mean_dict:
                        mean_dict[self.start_pos_list[i]] = coverage_values
                    else:
                        mean_dict[self.start_pos_list[i]] = mean_dict[self.start_pos_list[i]] + coverage_values

                    if self.start_pos_list[i] not in over_1_dict:
                        over_1_dict[self.start_pos_list[i]] = over_1_list
                    else:
                        over_1_dict[self.start_pos_list[i]] = over_1_dict[self.start_pos_list[i]] + over_1_list

                    if self.start_pos_list[i] not in over_5_dict:
                        over_5_dict[self.start_pos_list[i]] = over_5_list
                    else:
                        over_5_dict[self.start_pos_list[i]] = over_5_dict[self.start_pos_list[i]] + over_5_list

                    if self.start_pos_list[i] not in over_10_dict:
                        over_10_dict[self.start_pos_list[i]] = over_10_list
                    else:
                        over_10_dict[self.start_pos_list[i]] = over_10_dict[self.start_pos_list[i]] + over_10_list

                    if self.start_pos_list[i] not in over_15_dict:
                        over_15_dict[self.start_pos_list[i]] = over_15_list
                    else:
                        over_15_dict[self.start_pos_list[i]] = over_15_dict[self.start_pos_list[i]] + over_15_list

                    if self.start_pos_list[i] not in over_20_dict:
                        over_20_dict[self.start_pos_list[i]] = over_20_list
                    else:
                        over_20_dict[self.start_pos_list[i]] = over_20_dict[self.start_pos_list[i]] + over_20_list

                    if self.start_pos_list[i] not in over_25_dict:
                        over_25_dict[self.start_pos_list[i]] = over_25_list
                    else:
                        over_25_dict[self.start_pos_list[i]] = over_25_dict[self.start_pos_list[i]] + over_25_list

                    if self.start_pos_list[i] not in over_30_dict:
                        over_30_dict[self.start_pos_list[i]] = over_30_list
                    else:
                        over_30_dict[self.start_pos_list[i]] = over_30_dict[self.start_pos_list[i]] + over_30_list

                    if self.start_pos_list[i] not in over_50_dict:
                        over_50_dict[self.start_pos_list[i]] = over_50_list
                    else:
                        over_50_dict[self.start_pos_list[i]] = over_50_dict[self.start_pos_list[i]] + over_50_list

                    if self.start_pos_list[i] not in over_100_dict:
                        over_100_dict[self.start_pos_list[i]] = over_100_list
                    else:
                        over_100_dict[self.start_pos_list[i]] = over_100_dict[self.start_pos_list[i]] + over_100_list
                    #print("Finish updating all dict ")
                    #print(datetime.now())

        all_pos_list = []
        all_mean_list = []
        all_over_1_list = []
        all_over_5_list = []
        all_over_10_list = []
        all_over_15_list = []
        all_over_20_list = []
        all_over_25_list = []
        all_over_30_list = []
        all_over_50_list = []
        all_over_100_list = []

        for i in range(len(self.chrom_list)):
            if self.chrom_list[i] == chrom:
                all_pos_list.extend(list(range(self.start_pos_list[i] + 1, self.end_pos_list[i] + 1)))
                all_mean_list.extend(mean_dict[self.start_pos_list[i]].tolist())
                all_over_1_list.extend(over_1_dict[self.start_pos_list[i]].tolist())
                all_over_5_list.extend(over_5_dict[self.start_pos_list[i]].tolist())
                all_over_10_list.extend(over_10_dict[self.start_pos_list[i]].tolist())
                all_over_15_list.extend(over_15_dict[self.start_pos_list[i]].tolist())
                all_over_20_list.extend(over_20_dict[self.start_pos_list[i]].tolist())
                all_over_25_list.extend(over_25_dict[self.start_pos_list[i]].tolist())
                all_over_30_list.extend(over_30_dict[self.start_pos_list[i]].tolist())
                all_over_50_list.extend(over_50_dict[self.start_pos_list[i]].tolist())
                all_over_100_list.extend(over_100_dict[self.start_pos_list[i]].tolist())
        
        #print("Finish merging all values from same chrom ")
        #print(datetime.now())
        cols = [
            "chrom",
            "pos",
            "mean",
            "median",
            "over_1",
            "over_5",
            "over_10",
            "over_15",
            "over_20",
            "over_25",
            "over_30",
            "over_50",
            "over_100",
        ]

        df = pd.DataFrame(columns=cols)
        df["pos"] = all_pos_list
        df["chrom"] = chrom
        df["mean"] = all_mean_list
        df["median"] = -9
        df["over_1"] = all_over_1_list
        df["over_5"] = all_over_5_list
        df["over_10"] = all_over_10_list
        df["over_15"] = all_over_15_list
        df["over_20"] = all_over_20_list
        df["over_25"] = all_over_25_list
        df["over_30"] = all_over_30_list
        df["over_50"] = all_over_50_list
        df["over_100"] = all_over_100_list

        df["mean"] = df["mean"] / len(self.bw_file_list)
        df["over_1"] = df["over_1"] / len(self.bw_file_list)
        df["over_5"] = df["over_5"] / len(self.bw_file_list)
        df["over_10"] = df["over_10"] / len(self.bw_file_list)
        df["over_15"] = df["over_15"] / len(self.bw_file_list)
        df["over_20"] = df["over_20"] / len(self.bw_file_list)
        df["over_25"] = df["over_25"] / len(self.bw_file_list)
        df["over_30"] = df["over_30"] / len(self.bw_file_list)
        df["over_50"] = df["over_50"] / len(self.bw_file_list)
        df["over_100"] = df["over_100"] / len(self.bw_file_list)

        #print("Finish creating dataframe for chromosome " + chrom + " ")
        #print(datetime.now())

        #print("Start writing gzip tsv file ")
        #print(datetime.now())
        df.to_csv(output_file + "_" + chrom + ".tsv.gz", index=False, sep="\t", compression="gzip")
        #print("Finish writing gzip tsv file ")
        #print(datetime.now())
        print("Output is written in " + output_file + "_" + chrom + ".tsv.gz")


def main():
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument("-bw", action="store", type=str, required=True, help="bigWig file list")
    my_parser.add_argument("-bed", action="store", type=str, required=True, help="bed file")
    my_parser.add_argument("-out", action="store", type=str, required=True, help="output file")

    args = my_parser.parse_args()

    my_bigwig_to_cov_sum = BigWigTOCoverageSum()
    my_bigwig_to_cov_sum.parse_bw_file_list(args.bw)
    my_bigwig_to_cov_sum.parse_bed_file(args.bed)
    my_bigwig_to_cov_sum.create_coverage_summary(args.out)


if __name__ == "__main__":
    main()
