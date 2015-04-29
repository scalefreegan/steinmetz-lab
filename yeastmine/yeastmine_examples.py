#!/usr/bin/env python

"""Examples of querying yeastmine with intermine webservice"""

__author__ = "Aaron Brooks"
__copyright__ = "Copyright 2015, "
__credits__ = ["Aaron Brooks"]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Aaron Brooks"
__email__ = "aaron.brooks@embl.de"
__status__ = "Development"

from intermine.webservice import Service
import pandas as pd

service = Service("http://yeastmine.yeastgenome.org/yeastmine/service")

#-------------------------------------------------------------------#
# Gene Info
#-------------------------------------------------------------------# 
gene = service.model.Gene.where(symbol = 'HFA1').first()

print gene.symbol + "\n" + gene.description
print gene

#-------------------------------------------------------------------#
# Model templates
#-------------------------------------------------------------------#
template = service.get_template("Gene_Pathways")
for row in template.results(A={"symbol":"HFA1"}):
	print row

#-------------------------------------------------------------------#
# Query
#-------------------------------------------------------------------#
query = service.new_query("Gene")
query.add_view("primaryIdentifier","name","symbol","pathways.name")
query.add_constraint("Gene", "LOOKUP", "HFA1")
for row in query.rows():
	print row

# The view specifies the output columns
query.add_view(
    "primaryIdentifier", "secondaryIdentifier", "symbol", "name", "sgdAlias",
    "regulationSummary.summaryParagraph",
    "regulationSummary.publications.pubMedId",
    "regulationSummary.publications.citation"
)

# Uncomment and edit the line below (the default) to select a custom sort order:
# query.add_sort_order("Gene.primaryIdentifier", "ASC")

# You can edit the constraint values below
query.add_constraint("Gene", "IN", "ALL_Verified_Uncharacterized_Dubious_ORFs", code = "A")

# Uncomment and edit the code below to specify your own custom logic:
# query.set_logic("A")

# Outer Joins
# (display properties of these relations if they exist,
# but also show objects without these relationships)
query.outerjoin("regulationSummary.publications")

for row in query.rows():
    print row["primaryIdentifier"], row["secondaryIdentifier"], row["symbol"], row["name"], \
        row["sgdAlias"], row["regulationSummary.summaryParagraph"], \
        row["regulationSummary.publications.pubMedId"], \
        row["regulationSummary.publications.citation"]