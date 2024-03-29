#Creates a manifest file for qiime based on the names of samples in a csv

#Get the list from the user: In the format of Author_Diseased/Controls.csv
full_path = raw_input("What is the sample list for manifest file? ")

#Takes the name and creates the path based on the name
path_list = full_path.split("/")
sample_list = path_list[-1:]
naming_list = sample_list[0].split("_")
author = naming_list[0]
sample_type = naming_list[1].strip(".csv")
sample_names = []

#opens the file and creates a list of sample names
with open(full_path) as samples:
    lines = samples.readlines()
    for line in lines:
        line = line.strip("\n")
        sample_names.append(line)

manifest_path_on_maple = "/vol/mg-binning/Parkinsons/ManifestFiles/"

#Creates a manifest file with the paths on maple
with open(manifest_path_on_maple + author + "_" + sample_type + "_manifest.txt", "w+") as manifest_file:
    manifest_file.write("sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")
    for sample in sample_names:
        manifest_file.write(sample + "\t" + "/vol/mg-binning/Parkinsons/" + sample_type + "/" + author + "/" + sample + "/" + sample + "_1.fastq\t/vol/mg-binning/Parkinsons/" + sample_type + "/" + author + "/" + sample + "/" + sample + "_2.fastq")
        #If the sample is NOT the last one on the list, put a new line at the end.
        if sample != sample_names[:-1]:
            manifest_file.write("\n")