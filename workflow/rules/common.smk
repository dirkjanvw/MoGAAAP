def get_all_accessions():
    """
    Return all accession IDs.
    If haplotypes > 1, return "{accessionId}.hap{hap}" for each haplotype.
    """
    return [
        accession if haplotypes == 1 else f"{accession}.hap{hap}"
        for accession, haplotypes in zip(SAMPLES["accessionId"], SAMPLES["haplotypes"])
        for hap in range(1,haplotypes + 1)
    ]

def get_clean_accession_id(asmname):
    """
    Return the clean accession ID for a given asmname.
    """
    return asmname.split(".hap")[0]

def get_all_accessions_from_asmset(asmset, minimum=2):
    """
    Return all accession IDs for a given set, while taking care of haplotypes if needed
    """
    all_accessions = zip(SAMPLES["accessionId"], SAMPLES["haplotypes"])
    accessions = [
        accession if haplotypes == 1 else f"{accession}.hap{hap}"
        for accession, haplotypes in all_accessions
        for hap in range(1,haplotypes + 1)
        if accession in config["set"][asmset]
    ]
    if len(accessions) < minimum:
        print(f'WARNING: set "{asmset}" appears to contain only {len(accessions)} entry! Treating as if empty.')
        accessions = []
    return accessions

def get_hifi(wildcards):
    return SAMPLES[SAMPLES["accessionId"] == get_clean_accession_id(wildcards.asmname)]["hifi"].values.item().split(";")

def get_ont(wildcards):
    return SAMPLES[SAMPLES["accessionId"] == get_clean_accession_id(wildcards.asmname)]["ont"].values.item().split(";")

def get_illumina_1(wildcards):
    return SAMPLES[SAMPLES["accessionId"] == get_clean_accession_id(wildcards.asmname)]["illumina_1"].values.item()

def get_illumina_2(wildcards):
    return SAMPLES[SAMPLES["accessionId"] == get_clean_accession_id(wildcards.asmname)]["illumina_2"].values.item()

def has_hic(wildcards):
    return not SAMPLES[SAMPLES["accessionId"] == get_clean_accession_id(wildcards.asmname)]["hic_1"].isnull().values.item()

def get_hic_1(wildcards):
    return SAMPLES[SAMPLES["accessionId"] == get_clean_accession_id(wildcards.asmname)]["hic_1"].values.item()

def get_hic_2(wildcards):
    return SAMPLES[SAMPLES["accessionId"] == get_clean_accession_id(wildcards.asmname)]["hic_2"].values.item()

def get_haplotypes(wildcards):
    return int(SAMPLES[SAMPLES["accessionId"] == get_clean_accession_id(wildcards.asmname)]["haplotypes"].values.item())

def get_species_name(wildcards):
    return SAMPLES[SAMPLES["accessionId"] == get_clean_accession_id(wildcards.asmname)]["speciesName"].values.item()

def get_taxid(wildcards):
    return SAMPLES[SAMPLES["accessionId"] == get_clean_accession_id(wildcards.asmname)]["taxId"].values.item()

def get_reference_id(asmname):
    return SAMPLES[SAMPLES["accessionId"] == get_clean_accession_id(asmname)]["referenceId"].values.item()
