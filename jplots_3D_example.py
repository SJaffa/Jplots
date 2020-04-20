"""
Example script for a simple Jplots analysis

"""
import j3d as jp

#Read in data to a 2D or 3D numpy array
data=jp.read_data(source='example_3d')

#Extract the structures of interest to  list of masks
structs=jp.get_structures(data, method='dendrogram_ppp')

# Find J values of the structures
results=jp.j3d(structs)

#Output results
jp.out_to_file(results)
jp.out_to_plot(results)
