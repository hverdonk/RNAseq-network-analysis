import CleanData as cda
import py4cytoscape as p4c

# # Test that cytoscape is connected
# p4c.cytoscape_ping()
# print(p4c.cytoscape_version_info())

# import differential expression data and identify differentially expressed genes
gene_datafile = ''
data = cda.CleanData(gene_datafile)

print(data.genes())

# list_of_genes = []

# query STRING database to build a network of protein-protein interactions using the differentially expressed genes
string_cmd = 'string protein query query=Rgs4 cutoff=0.4 species="Mus musculus" limit=150' # % list_of_genes
p4c.commands_run(string_cmd)

test_table = p4c.get_table_columns('node','display name')
print(test_table)

p4c.network_views.export_image(resolution=300)

# # exclude low degree nodes, e.g., those with only 0, 1 or 2 connections
# p4c.create_degree_filter('degree filter', [0,2], 'IS_NOT_BETWEEN')

# # create a subnetwork of the selected set of nodes and all relevant edges
# p4c.create_subnetwork(subnetwork_name='highly connected nodes')



# You MUST open cytoscape on your local machine to run network_finder.py


# Cyni Toolbox is the app that contains the arachne algorithm
# There is also a "NetworkAnalyzer" bundle that exists

# stringApp is it's own bundle (aka app)

# STEP 1: build a regulatory network for XX and for XY individuals (and for female and male individuals)
    # FOR PPI: use the MusMusculus_binary_hq.txt file as the interaction network input (or download a different interaction network file 
    # from HINT: High quality INTeractomes)
        # you may have to reformat the file so that cytoscape will accept it
    # FOR GRN: you'll need a list of raw gene expression data (which we have) and a list of transcription factors to use as hubs (which I can ask Nora for)