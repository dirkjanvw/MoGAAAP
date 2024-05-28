def get_all_accessions():
    return SAMPLES["accessionId"].values

def get_hifi(wildcards):
    return SAMPLES[SAMPLES["accessionId"] == wildcards.asmname]["hifi"].values.item()

def get_ont(wildcards):
    return SAMPLES[SAMPLES["accessionId"] == wildcards.asmname]["ont"].values.item()

def get_illumina_1(wildcards):
    return SAMPLES[SAMPLES["accessionId"] == wildcards.asmname]["illumina_1"].values.item()

def get_illumina_2(wildcards):
    return SAMPLES[SAMPLES["accessionId"] == wildcards.asmname]["illumina_2"].values.item()

def get_species_name(wildcards):
    return SAMPLES[SAMPLES["accessionId"] == wildcards.asmname]["speciesName"].values.item()

def get_taxid(wildcards):
    return SAMPLES[SAMPLES["accessionId"] == wildcards.asmname]["taxId"].values.item()

def get_reference_id(wildcards):
    return SAMPLES[SAMPLES["accessionId"] == wildcards.asmname]["referenceId"].values.item()

def get_ref_genome(wildcards):
    referenceId = SAMPLES[SAMPLES["accessionId"] == wildcards.asmname]["referenceId"].values.item()
    return config["reference_genomes"][referenceId]["genome"]

def get_ref_annotation(wildcards):
    referenceId = SAMPLES[SAMPLES["accessionId"] == wildcards.asmname]["referenceId"].values.item()
    return config["reference_genomes"][referenceId]["annotation"]

def get_ref_chr(wildcards):
    referenceId = SAMPLES[SAMPLES["accessionId"] == wildcards.asmname]["referenceId"].values.item()
    return config["reference_genomes"][referenceId]["chromosomes"]
