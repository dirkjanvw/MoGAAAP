def get_all_accessions():
    all_accessions = set()
    for accession in config["reads"]["hifi"].keys():
        all_accessions.add(accession)
    if "ont" in config["reads"]:
        if config["reads"]["ont"]:
            for accession in config["reads"]["ont"].keys():
                all_accessions.add(accession)
    return all_accessions

def get_ref_genome(wildcards):
    if len(config["ref_genome"]) == 1:
        return [config["ref_genome"][name] for name in config["ref_genome"]]
    else:
        raise ValueError("Exactly one reference genome has to be specified.")
