import CleanData
from py2cytoscape import cyrest
from py2cytoscape.cyrest.base import api

# You MUST open cytoscape on your local machine to run network_finder.py

# Set up cytoscape connection
cytoscape=cyrest.cyclient()
cytoscape.session.new()

# Cyni Toolbox is the app that contains the arachne algorithm
# There is also a "NetworkAnalyzer" bundle that exists

# stringApp is it's own bundle (aka app)

# STEP 1: build a regulatory network for XX and for XY individuals (and for female and male individuals)
    # FOR PPI: use the MusMusculus_binary_hq.txt file as the interaction network input (or download a different interaction network file 
    # from HINT: High quality INTeractomes)
        # you may have to reformat the file so that cytoscape will accept it
    # FOR GRN: you'll need a list of raw gene expression data (which we have) and a list of transcription factors to use as hubs (which I can ask Nora for)