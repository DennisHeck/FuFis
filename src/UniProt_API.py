import pandas as pd
import requests
from multiprocessing import Pool
from timeit import default_timer as clock
from time import sleep


def api_helper(args):
    """
    Function to enable parallelization of the UniProt API. The error handling is very sub-optimal, but there's a
    certain randomness to the API that's hard to capture.
    """
    base_url, n_batch, protein_list = args
    print('Batch', n_batch)
    domains = {g: set() for g in protein_list}
    regions = {g: set() for g in protein_list}
    found_proteins = set()
    failed_reqs = {}

    # Careful with the brackets for a proper OR condition only across the genes.
    uniprot_url = base_url.format("%29+OR+%28gene%3A".join(protein_list))
    if len(uniprot_url) > 2000:
        print("ERROR, URL too long with {} characters".format(len(uniprot_url)))
        return
    uniprot_call = requests.get(uniprot_url)
    if uniprot_call.status_code == 504:  # Means Gateway time out, will wait, will try three times.
        for wait_i in range(3):
            sleep(2)
            uniprot_call = requests.get(uniprot_url)
            if uniprot_call.status_code == 200:
                break
        if uniprot_call.status_code != 200:  # If it still fails, stop here.
            failed_reqs[n_batch] = {"protein_names": protein_list, 'Exception': None,
                                    'status_code': uniprot_call.status_code}
            return None, None, None, failed_reqs
    if uniprot_call.status_code == 200:
        try:
            uniprot_call = uniprot_call.json()
        # Apparently there's a random JSON error showing up randomly. Try once again. The try-except here is really
        # ugly, but I wasn't able to reproduce the error consistently.
        except Exception as e:
            sleep(2)
            uniprot_call = requests.get(uniprot_url)
            if uniprot_call.status_code == 200:
                try:
                    uniprot_call = uniprot_call.json()
                except Exception as e:
                    print("JSON Format Error for batch {}:".format(n_batch), e)
                    failed_reqs[n_batch] = {"protein_names": protein_list, 'Exception': e,
                                            'status_code': uniprot_call.status_code}
                    return None, None, None, failed_reqs
    else:
        print("ERROR from UniProt's REST API, status code", uniprot_call.status_code, 'Batch ' + str(n_batch))
        return None, None, found_proteins, failed_reqs

    if uniprot_call['results']:
        for result in uniprot_call['results']:
            protein = result['genes'][0]['geneName']['value']
            # For whatever reason there are sometimes multiple entries, if there are others with the same alias.
            if protein in protein_list:
                domains[protein] |= set([x['description'] for x in result['features'] if x['type'] == 'Domain'])
                regions[protein] |= set([x['description'] for x in result['features'] if x['type'] == 'Region'])
                found_proteins.add(protein)
    else:
        failed_reqs[n_batch] = {"protein_names": protein_list, 'Exception': "No exception but empty JSON",
                                'status_code': uniprot_call.status_code}
    return domains, regions, found_proteins, failed_reqs


def uniprot_domains(protein_names, species='human', n_cores=1):
    """
    Uses the UniProt API to look up the annotated domains and regions for a list of protein names/gene names.
    The domains will only contain the ones from PROSITE, the regions also other annotations. Only entries that are
    reviewed are considered. Also, the proteins must exist like that in UniProt, aliases are not mapped.
    CARE: The UniProt API is unreliable, time-outs, random mal-formatted JSONs, empty results etc.

    Returns:
        protein_domains: Dictionary of {protein: set(domains)} with the PROSITE domains. Can be empty.
        protein_regions: Dictionary of {protein: set(regions)}. Can be empty.
        missed_proteins: Set of proteins where no matching entry was found in UniProt.
        failed_requests: Dictionary of {batch number: {Failure information}} in cases the request didn't work or where the JSON from the API was malformatted, not sure why and when this happens.
    """
    start = clock()
    # Prepare the uniprot ID.
    if species == "mouse":
        species_id = '10090'
    elif species == "human":
        species_id = '9606'
    else:
        print("ERROR: Unknown species", species, "currently allowed is human and mouse")
        return

    uniprot_base_url = 'https://rest.uniprot.org/uniprotkb/stream?format=json&query=%28%28%28gene%3A{}%29%29+AND+%28organism_id%3A'+species_id+'%29+AND+%28reviewed%3Atrue%29%29'

    batches = [protein_names[i:i + 60] for i in range(0, len(protein_names), 60)]
    print('Total batches:', len(batches))

    process_pool = Pool(processes=n_cores)
    fetched_batches = process_pool.map(api_helper, [[uniprot_base_url, n_b, batch] for n_b, batch in enumerate(batches)])
    process_pool.close()

    protein_domains = {k: v for d in [x[0] for x in fetched_batches if x[0]] for k, v in d.items()}
    protein_regions = {k: v for d in [x[1] for x in fetched_batches if x[1]] for k, v in d.items()}
    found_proteins = set.union(*[x[2] for x in fetched_batches])
    missed_proteins = set(protein_names) - found_proteins
    failed_requests = {k: v for d in [x[3] for x in fetched_batches if x[3]] for k, v in d.items()}

    print('UniProtAPI queried', clock() - start)
    return protein_domains, protein_regions, missed_proteins, failed_requests



