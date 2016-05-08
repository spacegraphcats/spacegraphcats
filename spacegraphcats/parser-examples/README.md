# GXT format definition

The `gxt` format (short for graph txt) is two csv files for the vertices and the edges concatenated. Between these two is a single blank line. The only required fields (in this order) are `id` and `size` for nodes and `src` and `dest` for edges. Node ids (and thus the src and dest in edges) have to be integers.
