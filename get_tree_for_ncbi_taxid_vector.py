def get_tree_for_ncbi_taxid_vector(x):

   from ete3 import NCBITaxa
   from os import path
   import shutil
   
   # update database as required 
   # ncbi.update_taxonomy_database() # creates  /Users/paul/.etetoolkit/*
   # shutil.move("/Users/paul/taxa.tab", "/Volumes/HGST1TB/Users/paul/Sequences/References/ete3/taxa.tab")
   # shutil.move("/Users/paul/taxdump.tar.gz", "/Volumes/HGST1TB/Users/paul/Sequences/References/ete3/taxdump.tar.gz")
  
   # set path to database location
   # ncbi = NCBITaxa(taxdump_file=path.abspath("/Volumes/HGST1TB/Users/paul/Sequences/References/ete3/taxdump.tar.gz")
  
   ncbi = NCBITaxa()
   # ncbi = NCBITaxa(dbfile="/path/to/taxa.sqlite") https://github.com/etetoolkit/ete/issues/295
   
   # for testing only
   x = [9606, 9598, 10090, 7707, 8782]
   
   # create tree as per
   #  http://etetoolkit.org/docs/latest/tutorial/tutorial_ncbitaxonomy.html#id1
   tree = ncbi.get_topology(x)
   print tree.get_ascii(attributes=["sci_name", "rank"])
   
   # return Newick tree
   return tree.write(format=1)
   
