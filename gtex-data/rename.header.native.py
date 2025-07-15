# Native Python implementation to relabel GTEx columns using tissue names

# Step 1: Read the sample attributes file into a dictionary
sample_to_tissue = {}
with open("GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt", 'r') as attr_file:
    header = attr_file.readline().strip().split('\t')
    sampid_index = header.index('SAMPID')
    smtsd_index = header.index('SMTSD')
    for line in attr_file:
        fields = line.strip().split('\t')
        sampid = fields[sampid_index]
        smtsd = fields[smtsd_index].replace(' ', '')
        sample_to_tissue[sampid] = smtsd

# Step 2: Read the GCT file, skip the first two metadata lines
with open("GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct", 'r') as gct_file:
    _ = gct_file.readline()  # Skip first metadata line
    _ = gct_file.readline()  # Skip second metadata line
    header = gct_file.readline().strip().split('\t')

    # Rename the sample columns
    new_header = []
    for col in header:
        if col in sample_to_tissue:
            new_header.append(f"{sample_to_tissue[col]}-{col}")
        else:
            new_header.append(col)

    # Step 3: Write the modified content to output file
    with open("GTEx_gene_tpm_by_tissue.tsv", 'w') as out_file:
        out_file.write('\t'.join(new_header) + '\n')
        for line in gct_file:
            out_file.write(line)
