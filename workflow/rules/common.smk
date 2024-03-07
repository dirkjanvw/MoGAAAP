def get_all_accessions():
    all_accessions = []
    for accession in config["reads"]["hifi"].keys():
        all_accessions.append(accession)
    for accession in config["reads"]["ont"].keys():
        all_accessions.append(accession)
    return all_accessions
