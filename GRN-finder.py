import utils
import subprocess
import pandas as pd
import sys
import traceback
import py4cytoscape as p4c
from bioinfokit import visuz

# BEFORE RUNNING THIS PROGRAM, fill in appropriate filepaths for run-ARACHNE-AP.sh

# check that cytoscape is connected (make sure the cytoscape app is open before running GRN-finder.py)
p4c.cytoscape_ping()
p4c.cytoscape_version_info()

# import differential gene expression data
male_DEGs = ''
male_data = utils.DEgenes(DE_genes=male_DEGs)
mgenes = len(male_data.get_DEGs())

female_DEGs = ''
female_data = utils.DEgenes(DE_genes=female_DEGs)
fgenes = len(female_data.get_DEGs())

print('There are %s differentially expressed genes in the male dataset and %s differentially expressed genes in the female dataset.' % (mgenes, fgenes))

# make volcano plot
print('Generating volcano plots.')
# Genes with low count numbers and outlier samples are prefiltered by DESeq2 and have a padj of "NA"
female = female_data.full_data.dropna(subset=['padj'])
male = male_data.full_data.dropna(subset=['padj'])
visuz.gene_exp.volcano(df=female, lfc="log2FoldChange", pv="padj", sign_line=True, color=("#FFA500", "#bcbcbc", "#2986CC"), figname='female-DEGs.png')
visuz.gene_exp.volcano(df=male, lfc="log2FoldChange", pv="padj", sign_line=True, color=("#FFA500", "#bcbcbc", "#2986CC"), figname='male-DEGs.png')

# Estimate gene regulatory network using default ARACHNe algorithm settings
print('Using ARACHNe-AP to estimate gene regulatory network.')
subprocess.run(['bash', 'run-ARACHNe-AP.sh'])

# Import inferred gene regulatory network into Cytoscape
# cytoscape column types to assign to imported network columns (ordered by column index)
# 'source,target,edge attribute,edge attribute'
attr='s,t,ea,ea'
# given a large enough network, you may have to go into the cytoscape app itself and force the software to create a network view
# if you get an error at any point, try to import the network manually, create a network view, and run the rest of this script
print('Importing inferred network...')
p4c.import_network_from_tabular_file('network.txt', first_row_as_column_names=True, column_type_list=attr)
grn_suid = p4c.get_network_suid()
p4c.layout_network(layout_name='force-directed', network=grn_suid)

# clone network to do male and female DEG analyses
p4c.clone_network(network=grn_suid)
grn_female = 'network.txt_1'
p4c.clone_network(network=grn_suid)
grn_male = 'network.txt_2'

# set node appearance defaults
p4c.delete_style_mapping(style_name='default', visual_prop='NODE_LABEL')
p4c.set_node_shape_default('ELLIPSE')
p4c.set_node_color_default('#808080')
p4c.set_node_size_default(20)

# create subnetworks
print("Creating subnetworks...")

# select female DEGs (XYF-vs-XXF comparison)
p4c.clear_selection(network=grn_female)
fem_deg = list(female_data.get_DEGs())
fem_nodes = p4c.select_nodes(nodes=fem_deg, by_col='name', network=grn_female)
try:
    fem_nodes = fem_nodes['nodes']
except Exception:
    if len(fem_nodes) == 0:
        print("There are no differentially expressed genes present in the female network. Cannot build a female subnetwork.")
    else:
        traceback.print_exc()
    sys.exit(1)
# identify genes associated with female DEGs via network propagation
# large networks MUST have a view before running the code below
p4c.set_current_network(grn_female)
p4c.diffusion_basic()
p4c.create_subnetwork(nodes="selected", subnetwork_name="female DEG diffusion")
p4c.layout_network('force-directed', 'female DEG diffusion')
female_subnetwork = 'female DEG diffusion'
# make female DEGs more visible
p4c.clear_selection(network='female DEG diffusion')
female_node_names = p4c.node_suid_to_node_name(fem_nodes, network='female DEG diffusion') 
p4c.set_node_label_bypass(node_names=fem_nodes, new_labels=female_node_names, network='female DEG diffusion')
p4c.set_node_size_bypass(node_names=fem_nodes, new_sizes=80, network='female DEG diffusion')

# select male DEGs (XXM-vs-XYM comparison)
p4c.clear_selection(network=grn_male)
male_deg = list(male_data.get_DEGs())
male_nodes = p4c.select_nodes(nodes=male_deg, by_col='name', network=grn_male)
try:
    male_nodes = male_nodes['nodes']
except Exception:
    if len(male_nodes) == 0:
        print("There are no differentially expressed genes present in the male network. Cannot build a male subnetwork.")
    else:
        traceback.print_exc()
    sys.exit(1)
# identify genes associated with male DEGs via network propagation
# large networks MUST have a view before running the code below
p4c.set_current_network(grn_male)
p4c.diffusion_basic()
p4c.create_subnetwork(nodes="selected", subnetwork_name="male DEG diffusion")
p4c.layout_network('force-directed', 'male DEG diffusion')
male_subnetwork = 'male DEG diffusion'
# make male DEGs more visible
p4c.clear_selection(network='male DEG diffusion')
male_node_names = p4c.node_suid_to_node_name(male_nodes, network='male DEG diffusion') 
p4c.set_node_label_bypass(node_names=male_nodes, new_labels=male_node_names, network='male DEG diffusion')
p4c.set_node_size_bypass(node_names=male_nodes, new_sizes=80, network='male DEG diffusion')

# Map male and female log2 fold changes onto the subnetworks
print('Mapping log2 fold changes onto subnetworks...')
print('SettingWithCopyWarnings are generated by cytoscape API and can be ignored.\n')
female_lfc_padj = female_data.reduce_data(data='full', cols=['geneNames', 'log2FoldChange', 'padj'])
p4c.load_table_data(data=female_lfc_padj, data_key_column='geneNames', table='node', 
                    table_key_column='name', network=female_subnetwork)
male_lfc_padj = male_data.reduce_data(data='full', cols=['geneNames', 'log2FoldChange', 'padj'])
p4c.load_table_data(data=male_lfc_padj, data_key_column='geneNames', table='node', 
                    table_key_column='name', network=male_subnetwork)

# Map color gradient onto over- and underexpressed genes.
networks = [female_subnetwork, male_subnetwork]
for network in networks:
    p4c.set_node_color_mapping(**p4c.gen_node_color_map('log2FoldChange', p4c.palette_color_brewer_d_RdBu(reverse=True), 
                                                        mapping_type='c', network=network))

# export subnetwork differential expression images
p4c.export_image(filename='female-deg-subnetwork-mapping.png', type='PNG', network=female_subnetwork)
p4c.export_image(filename='male-deg-subnetwork-mapping.png', type='PNG', network=male_subnetwork)

### Network analysis ###
# Calculate betweenness centrality for both subnetworks
print('\nCalculating betweenness centralities...\n')
p4c.set_current_network(network=female_subnetwork)
p4c.analyze_network(directed=False)
p4c.set_current_network(network=male_subnetwork)
p4c.analyze_network(directed=False)

# clear old network appearances (label & size of DEGs, color mapping)
p4c.delete_style_mapping(style_name='default', visual_prop='NODE_FILL_COLOR')
p4c.clear_node_property_bypass(node_names=female_node_names, visual_property='NODE_LABEL', network=female_subnetwork)
p4c.clear_node_property_bypass(node_names=female_node_names, visual_property='NODE_SIZE', network=female_subnetwork)        
p4c.clear_node_property_bypass(node_names=male_node_names, visual_property='NODE_LABEL', network=male_subnetwork)
p4c.clear_node_property_bypass(node_names=male_node_names, visual_property='NODE_SIZE', network=male_subnetwork)

# map subnetwork betweenness centrality
networks = [female_subnetwork, male_subnetwork]
for network in networks:
    p4c.set_node_size_mapping(**p4c.gen_node_size_map('BetweennessCentrality', p4c.scheme_c_number_continuous(0, 100), network=network))
    p4c.set_node_color_mapping(**p4c.gen_node_color_map('BetweennessCentrality', p4c.palette_color_brewer_d_PuOr(), mapping_type='c', network=network))

# Select top 5 genes with highest betweenness centrality for each subnetwork
# female top 5 genes
fem_bcent = p4c.get_table_columns(columns='BetweennessCentrality', network=female_subnetwork)
top_fem_bcent=fem_bcent.sort_values(by=['BetweennessCentrality'], ascending=False)
top_fem_bcent=top_fem_bcent.index.tolist()[:5]
node_names = p4c.node_suid_to_node_name(top_fem_bcent, network=female_subnetwork) 
p4c.set_node_label_bypass(node_names=top_fem_bcent, new_labels=node_names, network=female_subnetwork)
print('Top five female subnetwork nodes with highest betweenness centrality:')
print(node_names)
print()

# male top 5 genes
male_bcent = p4c.get_table_columns(columns='BetweennessCentrality', network=male_subnetwork)
top_male_bcent=male_bcent.sort_values(by=['BetweennessCentrality'], ascending=False)
top_male_bcent=top_male_bcent.index.tolist()[:5]
node_names = p4c.node_suid_to_node_name(top_male_bcent, network=male_subnetwork) 
p4c.set_node_label_bypass(node_names=top_male_bcent, new_labels=node_names, network=male_subnetwork)
print('Top five male subnetwork nodes with highest betweenness centrality:')
print(node_names)

# export subnetwork betweenness centrality images
p4c.export_image(filename='female-betweenness-subnetwork-mapping.png', type='PNG', network=female_subnetwork)
p4c.export_image(filename='male-betweenness-subnetwork-mapping.png', type='PNG', network=male_subnetwork)

