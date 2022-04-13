#Creates a manifest file for qiime based on the names of samples in a csv

#Get the list from the user: In the format of Author_Diseased/Controls.csv
sample_list = input("What is the sample list for manifest file? ")

#Takes the name and creates the path based on the name
naming_list = sample_list.split("_")
author = naming_list[0]
sample_type = naming_list[1].strip(".csv")
sample_names = []

#opens the file and creates a list of sample names
with open(sample_list) as samples:
    lines = samples.readlines()
    for line in lines:
        line = line.strip("\n")
        sample_names.append(line)

#Creates a manifest file with the paths on maple
with open(f"{author}_{sample_type}_manifest.txt", "w+") as manifest_file:
    manifest_file.write("sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")
    for sample in sample_names:
        manifest_file.write(f"{sample}\t"
                            f"/vol/mg-binning/Parkinsons/{sample_type}/{author}/{sample}/{sample}_1.fastq\t"
                            f"/vol/mg-binning/Parkinsons/{sample_type}/{author}/{sample}/{sample}_2.fastq")
        #If the sample is NOT the last one on the list, put a new line at the end.
        if sample != sample_names[:-1]:
            manifest_file.write("\n")