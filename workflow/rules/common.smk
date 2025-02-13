def singularity_enabled():
    """
    Returns a boolean for whether custom singularity containers may be used.
    """
    return config["custom_singularity"]

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

def get_haplotype_accession_id(asmname):
    """
    Return the haplotype for a given asmname.
    """
    return "hap" + asmname.split(".hap")[1]

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

def get_best_wgstype(asmname):
    """
    Return the best wgstype for a given asmname: if Illumina is available, return "Illumina", otherwise return "HiFi".
    """
    if has_illumina(asmname):
        return "Illumina"
    return "HiFi"

def get_hifi(wildcards):
    return SAMPLES[SAMPLES["accessionId"] == get_clean_accession_id(wildcards.asmname)]["hifi"].values.item().split(";")

def has_ont(asmname):
    return not SAMPLES[SAMPLES["accessionId"] == get_clean_accession_id(asmname)]["ont"].isnull().values.item()

def get_ont(wildcards):
    return SAMPLES[SAMPLES["accessionId"] == get_clean_accession_id(wildcards.asmname)]["ont"].values.item().split(";")

def has_illumina(asmname):
    return not SAMPLES[SAMPLES["accessionId"] == get_clean_accession_id(asmname)]["illumina_1"].isnull().values.item()

def get_illumina_1(wildcards):
    return SAMPLES[SAMPLES["accessionId"] == get_clean_accession_id(wildcards.asmname)]["illumina_1"].values.item()

def get_illumina_2(wildcards):
    return SAMPLES[SAMPLES["accessionId"] == get_clean_accession_id(wildcards.asmname)]["illumina_2"].values.item()

def has_hic(asmname):
    return not SAMPLES[SAMPLES["accessionId"] == get_clean_accession_id(asmname)]["hic_1"].isnull().values.item()

def get_hic_1(wildcards):
    return SAMPLES[SAMPLES["accessionId"] == get_clean_accession_id(wildcards.asmname)]["hic_1"].values.item()

def get_hic_2(wildcards):
    return SAMPLES[SAMPLES["accessionId"] == get_clean_accession_id(wildcards.asmname)]["hic_2"].values.item()

def get_haplotypes(wildcards):
    return int(SAMPLES[SAMPLES["accessionId"] == get_clean_accession_id(wildcards.asmname)]["haplotypes"].values.item())

def get_haplotype_information(asmname):
    return int(SAMPLES[SAMPLES["accessionId"] == get_clean_accession_id(asmname)]["haplotypes"].values.item())

def get_species_name(wildcards):
    return SAMPLES[SAMPLES["accessionId"] == get_clean_accession_id(wildcards.asmname)]["speciesName"].values.item()

def get_taxid(wildcards):
    return SAMPLES[SAMPLES["accessionId"] == get_clean_accession_id(wildcards.asmname)]["taxId"].values.item()

def get_reference_id(asmname):
    return SAMPLES[SAMPLES["accessionId"] == get_clean_accession_id(asmname)]["referenceId"].values.item()

def has_assembly_location(asmname):
    return not SAMPLES[SAMPLES["accessionId"] == get_clean_accession_id(asmname)]["assemblyLocation"].isnull().values.item()

def get_assembly_location(asmname):
    return SAMPLES[SAMPLES["accessionId"] == get_clean_accession_id(asmname)]["assemblyLocation"].values.item()

def has_annotation_location(asmname):
    return not SAMPLES[SAMPLES["accessionId"] == get_clean_accession_id(asmname)]["annotationLocation"].isnull().values.item()

def get_annotation_location(asmname):
    """
    Obtain the annotation location for a given asmname within the annotation module.
    NB: DO NOT USE THIS FUNCTION TO OBTAIN THE ACTUAL ANNOTATION FILE, "final_output/{asmname}.full.gff" SHOULD BE USED INSTEAD.
    """
    if has_annotation_location(asmname):
        return SAMPLES[SAMPLES["accessionId"] == get_clean_accession_id(asmname)]["annotationLocation"].values.item()
    else:
        return f"results/{asmname}/4.annotation/03.combined/{asmname}.gff"
