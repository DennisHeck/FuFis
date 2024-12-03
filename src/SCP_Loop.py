import subprocess

"""Based on filesystem string matching, download files from the servers."""

fetch_samples = ['IHECRE00000009', 'IHECRE00000010.3', 'IHECRE00000027.3', 'IHECRE00000102.3', 'IHECRE00000174.3',
                 'IHECRE00000187.3,', 'IHECRE00000283.2', 'IHECRE00000308.2', 'IHECRE00000335.2', 'IHECRE00001457.1',
                 'IHECRE00003213.2']

fetch_samples = [f.split('.')[0] for f in fetch_samples]

# And sent the respective bigwigs to local.
local_folder = "/Users/dennis/Dev/IHEC_ABC/MultipleMS/SuspiciousInteractions/MS_CD6_filteredInteractions/"
server_folder = "/projects/apog/work/IHEC/K27_data/"
for sample in fetch_samples:
    print(sample)
    subprocess.call('scp dhecker@biocontact1:/'+server_folder+'/*'+sample+'*.pval0.01.500K.narrowPeak.gz ' + local_folder, shell=True)




