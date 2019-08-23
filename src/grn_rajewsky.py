import pandas as pd
from distributed import Client, LocalCluster
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2
from itertools import izip_longest
import sys

def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return izip_longest(*args, fillvalue=fillvalue)

in_file  = 'rawdata/rajewsky.dge.txt'
tf_file  = 'rawdata/go-regulators_rajewsky.txt'

# ex_matrix is a DataFrame with gene names as column names
sys.stderr.write("\nReading count matrix...")
ex_matrix = pd.read_csv(in_file, sep='\t', index_col=0, header=None).T
sys.stderr.write("done.\n")

# tf_names is read using a utility function included in Arboreto
sys.stderr.write("\nLoading putative transcription factors...")
tf_names = load_tf_names(tf_file)
sys.stderr.write("done.\n")

sys.stderr.write("\nStarting Dusk cluster...") 
local_cluster = LocalCluster(n_workers=32,
    threads_per_worker=1,
    memory_limit=8e10)
custom_client = Client(local_cluster)
sys.stderr.write("done.\n")

sys.stderr.write("\nPredicting co-expression network in chunks...\n") 
i = 0
for chunk in grouper(tf_names, 20):
    sys.stderr.write("Working on chunk %s\n" % str(i))
    network = grnboost2(expression_data=ex_matrix, tf_names=chunk, client_or_address=custom_client)
    network.to_csv("network_rajewsky_" + str(i) + ".csv", sep=",", header=False, index=False)
    i += 1
sys.stderr.write("done.\n")

sys.stderr.write("\n\n# All done\n")
