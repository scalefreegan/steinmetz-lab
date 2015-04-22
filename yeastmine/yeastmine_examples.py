#!/usr/bin/env python

"""Examples of querying yeastmine with intermine webservice"""

__author__ = "Aaron Brooks"
__copyright__ = "Copyright 2014, "
__credits__ = ["Aaron Brooks"]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Aaron Brooks"
__email__ = "aaron.brooks@embl.de"
__status__ = "Development"

from intermine.webservice import Service

yeastmine = Service('http://yeastmine.yeastgenome.org/yeastmine/')

#-------------------------------------------------------------------#
# Gene Info
#-------------------------------------------------------------------# 
rad54 = yeastmine.model.Gene.where(symbol = 'rad54').first()

print(rad54.name)
print(rad54)

#-------------------------------------------------------------------#
# Model templates
#-------------------------------------------------------------------#
template = yeastmine.get_template("Gene_Pathways")
for row in template.results(A={"symbol":"FAS1"}):
	print row

#-------------------------------------------------------------------#
# Query
#-------------------------------------------------------------------#
query = yeastmine.new_query()
query.add_view("Gene.symbol", "Gene.pathway.name")
query.add_constraint("Gene", "LOOKUP", "rad54")
for row in query.results():
	print row

new_list = service.create_list("some/file/with.ids", "Gene")
list_on_server = service.get_list("On server")
in_both = new_list & list_on_server
in_both.name = "Intersection of these lists"
for row in in_both:
do_something_with(row)