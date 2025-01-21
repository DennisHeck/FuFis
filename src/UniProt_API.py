import requests


def uniprot_domains(protein_names, species='human'):
    """
    Uses the UniProt API to look up the annotated domains and regions for a list of protein names/gene names.
    The domains will only contain the ones from PROSITE, the regions also other annotations. Only entries that are
    reviewed are considered. Also, the proteins must exist like that in UniProt, aliases are not mapped.

    Returns:
        protein_domains: Dictionary of {protein: set(domains)} with the PROSITE domains. Can be empty.
        protein_regions: Dictionary of {protein: set(regions)}. Can be empty.
        missed_proteins: List of proteins where no matching entry was found in UniProt.
    """

    # Prepare the uniprot ID.
    if species == "mouse":
        species_id = '10090'
    elif species == "human":
        species_id = '9606'
    else:
        print("ERROR: Unknown species", species, "currently allowed is human and mouse")
        return

    uniprot_base_url = 'https://rest.uniprot.org/uniprotkb/stream?format=json&query=%28%28%28gene%3A{}%29%29+AND+%28organism_id%3A'+species_id+'%29+AND+%28reviewed%3Atrue%29%29'

    protein_domains = {g: set() for g in protein_names}
    protein_regions = {g: set() for g in protein_names}
    found_proteins = set()
    batches = [protein_names[i:i + 50] for i in range(0, len(protein_names), 50)]
    for batch in batches:
        # Careful with the brackets for a proper OR condition only across the genes.
        uniprot_url = uniprot_base_url.format("%29+OR+%28gene%3A".join(batch))
        if len(uniprot_url) > 2000:
            print("ERROR, URL too long with {} characters".format(len(uniprot_url)))
            # return
        uniprot_call = requests.get(uniprot_url)
        if uniprot_call.status_code == 200:
            uniprot_call = uniprot_call.json()
        else:
            print("ERROR from UniProt's REST API, status code", uniprot_call.status_code)
            continue

        if uniprot_call['results']:
            for result in uniprot_call['results']:
                protein = result['genes'][0]['geneName']['value']
                # For whatever reason there are sometimes multiple entries, if there are others with the same alias.
                if protein in batch:
                    protein_domains[protein] |= set([x['description'] for x in result['features'] if x['type'] == 'Domain'])
                    protein_regions[protein] |= set([x['description'] for x in result['features'] if x['type'] == 'Region'])
                    found_proteins.add(protein)
    missed_proteins = set(protein_names) - found_proteins

    return protein_domains, protein_regions, missed_proteins



