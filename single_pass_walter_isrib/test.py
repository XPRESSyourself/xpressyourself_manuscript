import os



def get_files(
    directory,
    suffix,
    omit=[]):

    # Initialize blank file list to fill
    file_list = []

    # Walk through raw data files within given directory
    for file in next(os.walk(directory))[2]:
        if file.endswith(tuple(suffix)):
            file_list.append(file) # Do not append directory, files in list will be modified and output to different locations

    # Sort files in alphabetical order (helps in formatting the count tables correctly)
    file_list = sorted(file_list)

    # Get rid of bad grabs
    omit_drop = []
    if len(omit) > 0:
        for x in file_list:
            for o in omit:
                if str(o) in x:
                    omit_drop.append(x)

    for x in omit_drop:
        file_list.remove(x)

    return file_list


aligned = get_files('/Users/jordan/Desktop/7155666_export/alignments', ['_1_Aligned.sort.bam'])
deduped = get_files('/Users/jordan/Desktop/7155666_export/alignments', ['dedupRemoved.bam'])


for x in range(0, int(len(aligned)), 2):
    os.system('samtools merge /Users/jordan/Desktop/7155666_export/alignments/aligned_merged/' + str(aligned[x])[:-19] + '_merged.bam /Users/jordan/Desktop/7155666_export/alignments/' + str(aligned[x]) + ' /Users/jordan/Desktop/7155666_export/alignments/' + str(aligned[x+1]))
    os.system('samtools index /Users/jordan/Desktop/7155666_export/alignments/aligned_merged/' + str(aligned[x])[:-19] + '_merged.bam')



for x in range(0, int(len(deduped)), 2):
    os.system('samtools merge /Users/jordan/Desktop/7155666_export/alignments/deduped_merged/' + str(deduped[x])[:-19] + '_merged.bam /Users/jordan/Desktop/7155666_export/alignments/' + str(deduped[x]) + ' /Users/jordan/Desktop/7155666_export/alignments/' + str(deduped[x+1]))
    os.system('samtools index /Users/jordan/Desktop/7155666_export/alignments/deduped_merged/' + str(deduped[x])[:-19] + '_merged.bam')
    
