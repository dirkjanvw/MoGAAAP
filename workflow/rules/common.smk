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

def get_reference_id(asmname):
    return SAMPLES[SAMPLES["accessionId"] == asmname]["referenceId"].values.item()
