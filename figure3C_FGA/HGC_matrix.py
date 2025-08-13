import time
import numpy 

start_time = time.time()
	
########Loading full HGC to memory########
def load_connectome(max_add, dist_array, num_edges): #Load connectome ranking and distance data
	print("\nLoading ranking and distance connectomes data to memory\n")
	# file1 = open("/Users/yuval/Desktop/Connectome/9.05/combined_all_condensed_9.05.txt", 'r')
	file1 = open("combined_all_dist_rank_p_10.00.txt", 'r')
	counter = 0
	for line in file1:
		line1 = line.strip()
		line2 = line1.split("\t")
		#print line2[0] + "\t" + line2[1]
		if (line2[0] in genes_dic and line2[1] in genes_dic):
			print(line2[0] + "\t" + line2[1])
			t = (line2[0], line2[1])
			dist = float(line2[2])
			dist_array[counter] = dist
			if (dist > max_add):
				max_add = dist
			s = [dist]
			dic[t] = s
		counter = counter + 1
	temp_time = time.time()
	print("\nLoaded in " + str(time.time() - start_time) + "\tseconds\n")
	file1.close()
	return max_add, dist_array

def count_lines():
	print("\nCounting number of edges")
	# file1 = open("/Users/yuval/Desktop/Connectome/9.05/combined_all_condensed_9.05.txt", 'r')
	file1 = open("combined_all_dist_rank_p_10.00.txt", 'r')
	counter = 0
	for line in file1:
		counter = counter + 1
	file1.close()
	return counter

def verify_gene_name(gene_name): #check if gene is in full list
	is_there = 0 #Is gene in full list
	file2 = open("Nodes_binding_10.txt", "r")
	for temp_line in file2:
		temp_line2 = temp_line.rstrip("\n")
		temp_line3 = temp_line2.rstrip("\r")
		#print temp_line2 + "\t" + gene_name
		if (temp_line3 == gene_name):
			#print "yay"
			file2.close()
			return 1 #gene found
	file2.close()
	return 0 #gene not found

########Distance between two genes########
def find_distance(temp_gene1, temp_gene2):
	dist = -1
	t1 = (temp_gene1, temp_gene2)
	t2 = (temp_gene2, temp_gene1)
	if (temp_gene2 == temp_gene1):
		return 0
	if (t1 in dic):
		dist = dic[temp_gene1, temp_gene2][0]
	if (t2 in dic):
		dist = dic[temp_gene2, temp_gene1][0]
	return dist
	
###########Main program#########
print("\nCreating distance matix and edges list")
dic = {} #full connectome data
dic_genes = {}
genes_dic = {}

file1 = open("gene_list.txt")
file2 = open("distance_matrix.txt", 'w')

for line1 in file1:
	src = line1.rstrip("\n")
	src1 = src.rstrip("\r")
	is_there = verify_gene_name(src1)
	if (is_there == 1):
		#header = header + "\t" + src1
		genes_dic[src1] = 1
	else:
		exit("Gene " + src1 + " not found")

header = ''
for gene in genes_dic:
	header = header + "\t" + gene

header = header + "\n" 
print(header)
file2.write(header)
list_len = len(genes_dic)

counter = 0
max_len = 300 #maximum length for non-connected genes
num_edges = count_lines()
dist_array = numpy.zeros(num_edges, float)	
max_len, temp_dist_array = load_connectome(max_len, dist_array, num_edges)
final_dist_array = numpy.sort(temp_dist_array)
print("\nConnectome loaded, starting run\n")

counter1 = 0
for i in genes_dic:
	counter2 = 0
	row_str_dist = i
	row_str_p = i
	print(i + "\t" + str(counter1+1) + "\tout of\t" + str(list_len) + "\ttime in seconds:\t" + str(time.time()-start_time))
	for j in genes_dic:
		gene1 = i
		gene2 = j
		if (gene1 != gene2):
			try:
				dist = find_distance(gene1, gene2)
				row_str_dist = row_str_dist + "\t"+ str(dist) 
			except:
				dist = max_len
				row_str_dist = row_str_dist + "\t"+ str(dist) 
		else:
			row_str_dist = row_str_dist + "\t0"
		counter2 = counter2 + 1
	row_str_dist = row_str_dist + "\n"
	file2.write(row_str_dist)
	counter1 = counter1 + 1
	
file1.close()
file2.close()

