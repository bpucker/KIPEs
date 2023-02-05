### Boas Pucker ###
### b.pucker@tu-braunschweig.de ###
__version__ = "v3.2.0"	#converted to Python3

__reference__ = "Pucker et al., 2020: https://doi.org/10.3390/plants9091103 and Rempel&Pucker, 2023: https://doi.org/10.1101/2022.06.30.498365"

__usage__ = """
					KIPEs """ + __version__ + """("""+ __reference__ +""")
					
					Usage:
					python3 KIPEs3.py
					--baits <FOLDER_WITH_BAIT_SEQ_FILES>
					--positions <FOLDER_WITH_POSITION_FILES>|--residues
					--out <OUTPUT_DIR>
					--subject <SUBJECT_FILE (peptide,transcript,genomic sequences)> | --subjectdir <SUBJECT_FOLDER_WITH_SEQ_FILES>
					
					optional:
					--seqtype <TYPE_OF_SUBJECT_SEQUENCES(pep|rna|dna)>[pep]
					--checks <VALIDATE_INPUT_BEFORE_RUNNING_PIPELINE (on|off)>[on]
					--cpus <INT, NUMBER_OF_BLAST_THREADS>[10]
					--scoreratio <FLOAT, BLAST_SCORE_RATIO_CUTOFF>[0.3]
					--simcut <FLOAT, MINIMAL_BLAST_HIT_SIMILARITY_IN_PERCENT>[40.0]
					--genesize <INT, SIZE_OF_GENE_FOR_GROUPING_OF_EXON_HITS>[5000]
					--minsim <FLOAT, MINIMAL_SIMILARITY_IN_GLOBAL_ALIGNMENT>[0.4]
					--minres <FLOAT, MINIMAL_PROPORTION_OF_CONSERVED_RESIDUES>[off]
					--minreg <FLOAT, MINIMAL_PROPORTION_OF_CONSERVED_REGIONS>[off]
					--pathway <TEXT_FILE_SPECIFYING_GENE_ORDER>
					--possibilities <INT, NUMBER_OF_CONSIDERED_ENZYME_FUNCTIONS_PER_SEQ>[3]
					
					--mafft <PATH_TO_MAFFT>[mafft]
					--blastp <PATH_TO_AND_INCLUDING_BINARY>[blastp]
					--tblastn <PATH_TO_AND_INCLUDING_BINARY>[tblastn]
					--makeblastdb <PATH_TO_AND_INCLUDING_BINARY>[makeblastdb]
					
					--fasttree <PATH_TO_FASTTREE>(recommended for tree building)
					--iqtree <PATH_TO_IQ-TREE>
					
					--forester <ACTIVATES_GENE_TREE_CONSTRUCTION>[off]
					
					--exp <GENE_EXPRESSION_FILE_ACTIVATES_COEXPRESSION_ANALYSIS>[off]
					--rcut <CORRELATION_CUTOFF>[0.3]
					--pcut <P_VALUE_CUTOFF>[0.05]
					--minexp <MIN_EXP_PER_GENE>[30]
					
					bug reports and feature requests: b.pucker@tu-bs.de
					Complete documentation: https://github.com/bpucker/KIPEs
					"""

import os, glob, sys, time, re, math, subprocess, dendropy
from operator import itemgetter
try:
	from scipy import stats
except:
	sys.stdout.write(  "WARNING: scipy import failed. Analysis of co-expression will not be possible. Please ensure that scipy is installed to enable this option.\n" )
	sys.stdout.flush()
try:
	import hashlib
except ImportError:
	sys.stdout.write(  "WARNING: hashlib import failed. Calculation of md5sums for inputs files not possible. Please install hashlib to improve the documentation quality.\n" )
	sys.stdout.flush()

# --- end of imports --- #

def generate_query( fasta, blast_query_dir, nmt ):
	"""! @brief generate query file """
	
	ID = fasta.split('/')[-1].split('.')[0]
	output_file = blast_query_dir + ID + ".fasta"
	with open( output_file, "w" ) as out:
		with open( fasta, "r" ) as f:
			counter = 1
			nmt.write( f.readline().strip()[1:].replace( "\t", "   " ) + '\t' )	#remove header
			
			line = f.readline()
			seq = []
			while line:
				if line[0] == ">":
					out.write( '>' + ID + "_%_" + str( counter ).zfill( 3 ) + '\n' + "".join( seq ) + '\n' )
					nmt.write( ID + "_%_" + str( counter ).zfill( 3 ) + '\n' + line.strip()[1:].replace( "\t", "   " ) + "\t" )
					counter += 1
					seq = []
				else:
					seq.append( line.strip() )
				line = f.readline()
			out.write( '>' + ID + "_%_" + str( counter ).zfill( 3 ) + '\n' + "".join( seq ) + '\n' )
			nmt.write( ID + "_%_" + str( counter ).zfill( 3 ) + '\n' )
	return output_file


def load_self_BLAST_hit_scores( self_blast_result_file ):
	"""! @brief load self BLAST hit scores """
	
	self_BLAST_hit_scores = {}
	with open( self_blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				self_BLAST_hit_scores[ parts[0] ]
			except KeyError:
				self_BLAST_hit_scores.update( { parts[0]: float( parts[-1] ) } )
			line = f.readline()
	return self_BLAST_hit_scores


def load_BLAST_results( blast_result_file, self_scores, score_ratio_cutoff, similarity_cutoff, possibility_cutoff ):
	"""! @brief load BLAST results """
	
	valid_blast_hits = {}
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if float( parts[2] ) > similarity_cutoff:	#similarity is sufficient
				if float( parts[-1] ) > score_ratio_cutoff * self_scores[ parts[0] ]:	#substantial part of query is matched
					try:
						valid_blast_hits[ parts[1] ].append( { 'gene': parts[0].split('_%_')[0], 'score': float( parts[-1] ) } )
					except KeyError:
						valid_blast_hits.update( { parts[1]: [ { 'gene': parts[0].split('_%_')[0], 'score': float( parts[-1] ) } ] } )
			line = f.readline()
	
	# --- reduce BLAST hit number to given number of candidate possibilities ---- #
	final_valid_blast_hits = {}
	for key in list(valid_blast_hits.keys()):
		hits = sorted( valid_blast_hits[ key ], key=itemgetter( 'score' ) )[::-1]
		genes = []
		for hit in hits:
			if hit['gene'] not in genes:
				if len( genes ) < possibility_cutoff:
					genes.append( hit['gene'] )
		final_valid_blast_hits.update( { key: genes } )
	
	return final_valid_blast_hits


def  find_sisters_in_tree( tree_file, all_seq_IDs_in_tree ):
	"""! @brief find sister clade in phylogenetic tree 
			@note This function is based on the MYB_annotator: https://doi.org/10.1101/2021.10.16.464636
	"""
	
	if len( all_seq_IDs_in_tree ) > 0:
		# --- find node objects of reference genes --- #
		tree = dendropy.Tree.get_from_path( tree_file, "newick" )
		pdm = dendropy.PhylogeneticDistanceMatrix.from_tree( tree )
		my_mean_nearest_taxon_distance = pdm.mean_nearest_taxon_distance()
		
		ref_node_objects = []
		new_node_objects = []
		for node in tree.taxon_namespace:
			if node.label == "CANDIDATE":
				new_node_objects.append( node )
			else:
				ref_node_objects.append( node )
		data = []
		t1 = new_node_objects[0]
		for t2 in ref_node_objects:	#calculate distance to all other sequences in tree
			edge_distance = pdm.path_edge_count( t1, t2)
			patr_distance = pdm.patristic_distance( t1, t2 )
			data.append( { 'ref_gene': (t2.label).split('-%-')[0], 'edge': edge_distance, 'patr': patr_distance } )
		sisters = []
		for entry in sorted( data, key=itemgetter('edge', 'patr') ):
			if entry['ref_gene'] not in sisters:
				sisters.append( entry['ref_gene'] )
		return sisters	#returns a sorted list of reference gene function names (closest homolog first)
	else:
		sys.stdout.write( "ERROR: " + tree_file + "\n" )
		sys.stdout.flush()
		return []


def alignment_trimming( aln_file, cln_aln_file, occupancy ):
	"""! @brief remove all alignment columns with insufficient occupancy """
	
	alignment = load_alignment( aln_file, {} )
	# --- if there is an alignment (expected case) 
	if len( list(alignment.keys()) ) > 0:
		# --- identify valid residues in aligned sequences (columns with sufficient occupancy) --- #
		valid_index = []
		for idx, aa in enumerate( list(alignment.values())[0] ):
			counter = 0
			for key in list(alignment.keys()):
				if alignment[ key ][ idx ] != "-":
					counter += 1
			if counter / float( len( list(alignment.keys()) ) ) > occupancy:
				valid_index.append( idx )
		
		# --- generate new sequences --- #
		with open( cln_aln_file, "w" ) as out:
			for key in list(alignment.keys()):
				seq = alignment[ key ]
				new_seq = []
				for idx in valid_index:
					new_seq.append( seq[ idx ] )
				out.write( ">" + key + '\n' + "".join( new_seq ) + '\n' )
	# --- just in case the alignment file is empyt (is this possible?) ---#
	else:
		with open( cln_aln_file, "w" ) as out:
			out.write( "" )


def tree_based_classification( blast_result_file, self_scores, score_ratio_cutoff, similarity_cutoff, tree_tmp, treemethod, fasttree, iqtree, mafft, peps, ref_seqs, possibility_cutoff ):
	"""! @brief load BLAST results """
	
	final_results = {}
	# --- get all valid seqs --- #
	valid_blast_hits = {}
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if float( parts[2] ) > similarity_cutoff:	#similarity is sufficient
				if float( parts[-1] ) > score_ratio_cutoff * self_scores[ parts[0] ]:	#substantial part of query is matched
					try:
						valid_blast_hits[ parts[1] ].append( { 'id': parts[0], 'score': float( parts[-1] ), 'gene': parts[0].split('_%_')[0] } )
					except KeyError:
						valid_blast_hits.update( { parts[1]: [ { 'id': parts[0], 'score': float( parts[-1] ), 'gene': parts[0].split('_%_')[0] } ] } )
			line = f.readline()
	
	# --- prepare refseq data --- #
	ref_seq_seqs = {}
	ref_seq_ID_to_gene_mapping = {}
	for value in list(ref_seqs.values()):
		for entry in value:
			ref_seq_seqs.update( { entry['id']: entry['seq'] } )	#id or name?
			ref_seq_ID_to_gene_mapping.update( { entry['id']: entry['gene'] } )
	
	# --- process all valid BLAST results --- #
	for key in list(valid_blast_hits.keys()):
		if len( valid_blast_hits[ key ] ) == 1:
			final_results.update( { key: [ valid_blast_hits[ key ][0]['gene'] ] } )
		else:
			# --- build tree --- #
			seq_file = tree_tmp + key + ".fasta"
			all_seq_IDs_in_tree = []
			with open( seq_file, "w" ) as out:
				out.write( '>CANDIDATE\n' + peps[ key ] + '\n'  )
				for each in valid_blast_hits[ key ]:
					out.write( '>' + each['id'].replace('_',"-") + '\n' + ref_seq_seqs[ each['id'] ] + '\n'  )
					all_seq_IDs_in_tree.append( each['id'].replace('_',"-") )
			aln_file = seq_file + ".aln"
			p = subprocess.Popen( args= " ".join( [ mafft, seq_file, ">", aln_file, "2>", aln_file+".err" ] ), shell=True )
			p.communicate()
			
			cln_aln_file = aln_file + ".cln"
			occupancy = 0.3	#could be implemented as option later
			alignment_trimming( aln_file, cln_aln_file, occupancy )
			
			
			tree_file = cln_aln_file + ".tree"
			if treemethod == "fasttree":	#construct tree with FastTree2
				p = subprocess.Popen( args= " ".join( [ fasttree, "-wag -nosupport <", cln_aln_file, ">", tree_file, "2>", tree_file+".err" ] ), shell=True )
				p.communicate()
			elif treemethod == "iqtree":	#construct tree with IQ-TREE
				tmp_tree_file = cln_aln_file + ".treefile"	#".treefile" is appended to provided alignment file name
				p = subprocess.Popen( args= iqtree + " -nt 1 -s " + cln_aln_file, shell=True )	#bootstrapping: -alrt 1000 -bb 1000
				p.communicate()
				p = subprocess.Popen( args= "mv " + tmp_tree_file + " " + tree_file, shell=True )
				p.communicate()
			
			# --- find reference seqs sister --- #
			sisters = find_sisters_in_tree( tree_file, all_seq_IDs_in_tree )
			if len( sisters ) > 0:
				genes = []
				for s, sister in enumerate( sisters ):
					if s < possibility_cutoff:
						genes.append( sister )	#sister.split('_%_')[0]
				final_results.update( { key: genes } )
	return final_results


def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip()
		# if " " in header:
			# header = header.split(' ')[0]
			# if "\t" in header:
				# header = header.split('\t')[0]
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:]
					# if " " in header:
						# header = header.split(' ')[0]
						# if "\t" in header:
							# header = header.split('\t')[0]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def load_ref_seqs( query_file, name_mapping_table ):
	"""! @brief load all reference sequences """
	
	mapping_table = {}
	with open( name_mapping_table, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			mapping_table.update( { parts[1]: parts[0] } )
			line = f.readline()
	
	# --- load sequences --- #
	seqs = load_sequences( query_file )
	ref_seqs = {}
	for key in list(seqs.keys()):
		gene = key.split('_%_')[0]
		try:
			ref_seqs[ gene ].append( { 'id': key, 'gene': gene, 'seq': seqs[ key ], 'name': mapping_table[ key ] } )
		except KeyError:
			ref_seqs.update( { gene: [ { 'id': key, 'gene': gene, 'seq': seqs[ key ], 'name': mapping_table[ key ] } ] } )
	return ref_seqs


def load_alignment( aln_file, tmp_mapping ):
	"""! @brief load alignment and replace query IDs by real sequence names """
	
	sequences = {}
	
	with open( aln_file ) as f:
		header = f.readline()[1:].strip()
		try:
			header = tmp_mapping[ header ]
		except KeyError:
			pass
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:]
					try:
						header = tmp_mapping[ header ]
					except KeyError:
						pass
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def calculate_similarity( seq1, seq2 ):
	"""! @brief calculate similarity of two sequences with seq2 as reference """
	
	counter = 0
	for idx, aa in enumerate( seq1 ):
		if aa == seq2[ idx ] and aa != "-":
			counter += 1
	ref_len = float( len( seq2 ) - seq2.count( "-" ) )
	return counter / ref_len


def calculate_sim_matrix_per_gene( query_names_by_gene, candidates_by_gene, alignment_per_candidate ):
	"""! @brief calculate alignment-based similarity matrix per gene """
	
	sim_matrix_per_gene = {}
	#sim_matrix_per_gene: { 'CHS': { 'candidate1': { 'q1': xx, 'q2': xx, 'q3': xxx, ... }, 'candidate2': { ... }, .... }, 'CHI': { ... }, ... }
	#query_names_by_gene: { 'CHS': [ q1, q2, ... ], 'CHI': [...], ... }
	#candidates_by_gene: { 'CHS': [ candidate1, candidate2, ... ], 'CHI': [...], ... }
	#alignment_per_candidate: { 'candidate1': { 'candidate1': xxxxx, 'q1': xxxxxx, 'q2': xxxxx, ... }, 'candidate2': {...}, .... }
	for gene in list(candidates_by_gene.keys()):
		gene_sim_matrix = {}
		for candidate in candidates_by_gene[ gene ]:
			candiate_sim_matrix = {}
			for query in query_names_by_gene[ gene ]:
				try:
					sim = calculate_similarity( alignment_per_candidate[ candidate ][ gene ][ candidate ], alignment_per_candidate[ candidate ][ gene ][ query ] )
				except KeyError:
					sim = 0
				candiate_sim_matrix.update( { query: sim } )
			gene_sim_matrix.update( { candidate: candiate_sim_matrix } )
		sim_matrix_per_gene.update( { gene: gene_sim_matrix } )
	return sim_matrix_per_gene


def generate_global_alignments( mafft, peps, blast_hits, tmp_dir, ref_seqs ):
	"""! @brief generate and load global alignments, calculate similarity matrix """
	
	alignment_per_candidate = {}
	candidates_by_gene = {}
	for candidate in list(blast_hits.keys()):
		genes = blast_hits[ candidate ]
		for gene in genes:
			try:
				candidates_by_gene[ gene ].append( candidate )
			except KeyError:
				candidates_by_gene.update( { gene: [ candidate ] } )
			tmp_mapping = {}
			tmp_seq_file = tmp_dir + candidate + ".fasta"
			aln_file = tmp_dir + candidate + ".fasta.aln"
			# --- prepare multiple FASTA file --- #
			with open( tmp_seq_file, "w" ) as out:
				out.write( '>' + candidate + '\n' + peps[ candidate ] + '\n' )
				for ref in ref_seqs[ gene ]:
					out.write( '>' + ref['id'] + '\n' + ref['seq'] + '\n' )
					tmp_mapping.update( { ref['id']: ref['name'] } )
			# --- run alignment --- #
			p = subprocess.Popen( args= mafft + " " + tmp_seq_file + " > " + aln_file + " 2> " + aln_file+".err", shell=True )
			p.communicate()
			try:
				alignment_per_candidate[ candidate ].update( { gene: load_alignment( aln_file, tmp_mapping ) } )
			except KeyError:
				alignment_per_candidate.update( { candidate: { gene: load_alignment( aln_file, tmp_mapping ) } } )
	
	# --- get all query sequence names per gene --- #
	query_names_by_gene = {}
	for ref in [ x for sublist in list(ref_seqs.values()) for x in sublist]:	#walk through all query sequences
		try:
			query_names_by_gene[ ref['gene'] ].append( ref['name'] )
		except KeyError:
			query_names_by_gene.update( { ref['gene']: [ ref['name'] ] } )
	
	# --- calculate similarity matrix per gene --- #
	sim_matrix_per_gene = calculate_sim_matrix_per_gene( query_names_by_gene, candidates_by_gene, alignment_per_candidate )
	
	return alignment_per_candidate, sim_matrix_per_gene, candidates_by_gene


def generate_sim_matrix_output_files( sim_matrix_folder, sim_matrix_per_gene, subject_name_mapping_table ):
	"""! @brief generate one output file per gene """
	
	sim_per_pep = {}
	for gene in list(sim_matrix_per_gene.keys()):
		output_file = sim_matrix_folder + gene + "_sim_matrix.txt"
		with open( output_file, "w" ) as out:
			data = sim_matrix_per_gene[ gene ]	#similarity matrices for all candidates
			queries = sorted( list(data.values())[0].keys() )
			out.write( "\t".join( [ "candidate" ] + queries ) + '\n' )
			for candidate in sorted( data.keys() ):
				new_line = [ subject_name_mapping_table[ candidate ] ]
				for query in queries:
					new_line.append( 100.0*data[ candidate ][ query ] )
				try:
					sim_per_pep[ candidate ].update( { gene: sum( new_line[1:] ) / len( new_line[1:] ) } )
				except KeyError:
					sim_per_pep.update( { candidate: { gene: sum( new_line[1:] ) / len( new_line[1:] ) } } )
				out.write( "\t".join( map( str, new_line ) ) + '\n' )
	return sim_per_pep


def load_pos_data_per_gene( pos_data_files ):
	"""! @brief load conserved residue positions from given data file """
	
	pos_per_gene = {}
	regions_per_gene = {}
	for filename in pos_data_files:
		gene = filename.split('/')[-1].split('.')[0]
		with open( filename, "r" ) as f:
			#SEQ_NAME
			#R	X,Y,Z	1
			#R	X2
			#R	X3
			#	... 
			#R	X10
			#D	DOMAIN_SEQ	START	END
			#D	DOMAIN_SEQ	START	END	
			#...
			regions = []
			residues = []
			line = f.readline()
			while line:
				if line[0] != "#":
					if line[0] == "!":
						ref_seq = line.strip()[1:]
					else:
						parts = line.strip().split('\t')
						if len( parts ) > 0:
							if parts[0] == "D":
								if len( parts ) > 3:
									if len( parts ) > 4:
										comment = "".join( parts[4:] )
									else:
										comment = ""
									regions.append( { 'name': parts[1], 'start': int( parts[2] ), 'end': int( parts[3] ), 'comment': comment } )
							elif parts[0] == "R":
								if "," in parts[1]:
									res = parts[1].split(',')
								else:
									res = [ parts[1] ]
								if len( parts ) > 3:
									comment = "".join( parts[3:] )
								else:
									comment = ""
								residues.append( { 'aa': res, 'pos': int( parts[2] ), 'comment': comment } )		
				line = f.readline()
			regions_per_gene.update( { gene: { 'seq': ref_seq, 'regions': regions } } )
			pos_per_gene.update( { gene: { 'seq': ref_seq, 'residues': residues } } )
	return pos_per_gene, regions_per_gene


def get_alignment_pos( ref_aln, target_index ):
	"""! @brief get index in alignment based on target """
	
	index = -1
	counter = 0
	while counter < len( ref_aln ):
		if ref_aln[ counter ] != "-":
			index += 1
			if index == target_index:
				return counter
		counter += 1
	return "ERROR: target index exceeds sequence length!"


def check_alignment_for_cons_res( can_aln, ref_aln, residues ):
	"""! @brief inspect alignment of candidate at conserved residue positions """
	
	results = []
	extra_results = []
	for res in residues:
		alignment_pos = get_alignment_pos( ref_aln, res['pos']-1 )
		if can_aln[ alignment_pos ] == "-":
			results.append( "-" )
			extra_results.append( "-" )
		else:	
			results.append( can_aln[ alignment_pos ] in res['aa'] )	#multiple different amino acids might be permitted at one position
			extra_results.append( can_aln[ alignment_pos ] )
	return results, extra_results


def check_cons_res( cons_res_matrix_folder, pos_data_per_gene, alignment_per_candidate, candidates_by_gene, subject_name_mapping_table ):
	"""! @brief check all candidate sequences for conserved residues and generate result tables """
	
	cons_pos_per_pep = {}
	cons_pos_per_pep_extra = {}
	for gene in list(candidates_by_gene.keys()):
		candidates = candidates_by_gene[ gene ]
		try:
			info = pos_data_per_gene[ gene ]
			output_file = cons_res_matrix_folder + gene + "_conserved_residues.txt"
			residues = sorted( info['residues'], key=itemgetter('pos') )
			with open( output_file, "w" ) as out:
				header = [ "candidate" ]
				for each in residues:
					header.append( "/".join( each['aa'] ) + str( each['pos'] ) )
				out.write( "\t".join( header ) + '\n' )
				for candidate in candidates:
					can_aln = alignment_per_candidate[ candidate ][ gene ][ candidate ]
					ref_aln = alignment_per_candidate[ candidate ][ gene ][ info['seq'] ]
					results, extra_results = check_alignment_for_cons_res( can_aln, ref_aln, residues )
					out.write( "\t".join( map( str, [ subject_name_mapping_table[ candidate ] ] + results  ) ) + '\n' )
					try:
						if len( results ) > 0:
							cons_pos_per_pep[ candidate ].update( { gene: 100.0* results.count( True ) / len( results ) } )
							cons_pos_per_pep_extra[ candidate ].update( { gene: extra_results } )
						else:
							cons_pos_per_pep[ candidate ].update( { gene: 0.0 } )
							cons_pos_per_pep_extra[ candidate ].update( { gene: extra_results } )
					except KeyError:
						if len( results ) > 0:
							cons_pos_per_pep.update( { candidate: { gene: 100.0* results.count( True ) / len( results ) } } )
							cons_pos_per_pep_extra.update( { candidate: { gene: extra_results } } )
						else:
							cons_pos_per_pep.update( { candidate: { gene: 0.0 } } )
							cons_pos_per_pep_extra.update( { candidate: { gene: extra_results } } )
		except KeyError:
			sys.stdout.write( "ERROR: no information (conserved residues) available about gene: " + gene + "\n" )
			sys.stdout.flush()
	return cons_pos_per_pep, cons_pos_per_pep_extra


def check_alignment_for_cons_reg( can_aln, ref_aln, regions ):
	"""! @brief check conserved regions """
	
	results = []
	for reg in regions:
		alignment_start_pos = get_alignment_pos( ref_aln, reg['start']-1 )
		alignment_end_pos = get_alignment_pos( ref_aln, reg['end']-1 ) 
		
		#calculate similarity
		ref = ref_aln[ alignment_start_pos:alignment_end_pos+1 ]
		can = can_aln[ alignment_start_pos:alignment_end_pos+1 ]
		
		results.append( 100.0*calculate_similarity( can, ref ) )
	return results


def check_cons_reg( cons_reg_matrix_folder, regions_per_gene, alignment_per_candidate, candidates_by_gene, subject_name_mapping_table ):
	"""! @brief check all candidate sequences for conserved residues and generate result tables """
	
	cons_reg_per_pep = {}
	for gene in list(candidates_by_gene.keys()):
		candidates = candidates_by_gene[ gene ]
		try:
			info = regions_per_gene[ gene ]
			output_file = cons_reg_matrix_folder + gene + "_conserved_residues.txt"
			regions = sorted( info['regions'], key=itemgetter('start') )
			with open( output_file, "w" ) as out:
				header = [ "candidate" ]
				for each in regions:
					header.append( each['name'] + "," +  str( each['start'] ) + "-" + str( each['end'] ) )
				out.write( "\t".join( header ) + '\n' )
				for candidate in candidates:
					can_aln = alignment_per_candidate[ candidate ][ gene ][ candidate ]
					ref_aln = alignment_per_candidate[ candidate ][ gene ][ info['seq'] ]
					results = check_alignment_for_cons_reg( can_aln, ref_aln, regions )
					out.write( "\t".join( map( str, [ subject_name_mapping_table[ candidate ] ] + results  ) ) + '\n' )
					try:
						try:
							cons_reg_per_pep[ candidate ].update( { gene: sum( results ) / len( results ) } )
						except ZeroDivisionError:
							cons_reg_per_pep[ candidate ].update( { gene: 0.0 } )
					except KeyError:
						try:
							cons_reg_per_pep.update( { candidate: { gene: sum( results ) / len( results ) } } )
						except ZeroDivisionError:
							cons_reg_per_pep.update( { candidate: { gene: 0.0 } } )
		except KeyError:
			sys.stdout.write( "ERROR: no information (conserved regions) available about gene: " + gene + "\n" )
			sys.stdout.flush()
	return cons_reg_per_pep


def generate_subject_file( peptide_file, subject_name_file, subject_file ):
	"""! @brief rename all subject sequences at input """
	
	with open( peptide_file, "w" ) as out:
		with open( subject_name_file, "w" ) as out2:
			with open( subject_file, "r" ) as f:
				counter = 1
				out2.write( f.readline().strip()[1:].replace( "\t", "   " ) + '\t' )	#remove header
				line = f.readline()
				seq = []
				while line:
					if line[0] == ">":
						out.write( '>X' + str( counter ) + '\n' + "".join( seq ) + '\n' )
						out2.write( "X" + str( counter ) + '\n' + line.strip()[1:].replace( "\t", "   " ) + "\t" )
						counter += 1
						seq = []
					else:
						seq.append( line.strip() )
					line = f.readline()
				out.write( '>X' + str( counter ) + '\n' + "".join( seq ) + '\n' )
				out2.write(  "X" + str( counter ) + '\n' )


def load_subject_name_mapping_table( subject_name_file ):
	"""! @brief load subject name mapping table """
	
	mapping_table = {}
	with open( subject_name_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			mapping_table.update( { parts[1]: parts[0] } )
			line = f.readline()
	return mapping_table


# --- functions to handle RNA (transcriptome sequence) input --- #

def translate( seq, genetic_code ):
	"""! @brief translates the given nucleotide sequence into peptide and splits at each star (stop codon) """
	
	seq = seq.upper()
	
	peptide = []
	
	for i in range( int( len( seq ) / 3.0 ) ):
		codon = seq[i*3:i*3+3]
		try:
			peptide.append( genetic_code[ codon ] )
		except:
			peptide.append( "*" )
	return "".join( peptide )


def revcomp( seq ):
	"""! @brief construct reverse complement of sequence """
	
	new_seq = []
	
	bases = { 'a':'t', 't':'a', 'c':'g', 'g':'c' }
	for nt in seq.lower():
		try:
			new_seq.append( bases[nt] )
		except:
			new_seq.append( 'n' )
	return ''.join( new_seq[::-1] ).upper()


def translator( input_file, min_len_cutoff ):
	"""! @brief translate all sequences in given file in all frames """
	
	genetic_code = {	'CTT': 'L', 'ATG': 'M', 'AAG': 'K', 'AAA': 'K', 'ATC': 'I',
								'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T',
								'AGC': 'S', 'ACA': 'T', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N',
								'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H',
								'ACG': 'T', 'CCG': 'P', 'AGT': 'S', 'CAG': 'Q', 'CAA': 'Q',
								'CCC': 'P', 'TAG': '*', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C',
								'CGA': 'R', 'CCA': 'P', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R',
								'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G',
								'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'GAG': 'E', 'TCG': 'S',
								'TTA': 'L', 'GAC': 'D', 'TCC': 'S', 'GAA': 'E', 'TCA': 'S',
								'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A',
								'GTG': 'V', 'TTC': 'F', 'GTT': 'V', 'GCT': 'A', 'ACC': 'T',
								'TTG': 'L', 'CGT': 'R', 'TAA': '*', 'CGC': 'R'
							}
	
	sequences = load_sequences( input_file )
	peptides = []
	for seq in list(sequences.values()):
		# --- translate in all frames --- #
		frame1 = translate( seq, genetic_code ).split('*')
		frame2 = translate( seq[1:], genetic_code ).split('*')
		frame3 = translate( seq[2:], genetic_code ).split('*')
		rev_seq = revcomp( seq )
		frame4 = translate( rev_seq, genetic_code ).split('*')
		frame5 = translate( rev_seq[1:], genetic_code ).split('*')
		frame6 = translate( rev_seq[2:], genetic_code ).split('*')
		
		# --- pick all peps over length cutoff --- #
		for p in frame1+frame2+frame3+frame4+frame5+frame6:
			if len( p ) >= min_len_cutoff:
				peptides.append( p )
	return peptides


def translate_to_generate_pep_file( peptide_file, subject_name_file, subject, min_len_cutoff ):
	"""! @brief run in silico translation to generate proper input """
	
	#peptide_file needs to be final output
	#mapping table (subject_name_file) needs two columns with identical names
	#subject contains input sequences
	peptides = translator( subject, min_len_cutoff )
	with open( peptide_file, "w" ) as out:
		with open( subject_name_file, "w" ) as out2:
			for idx, pep in enumerate( peptides ):
				out.write( '>X' + str( idx ) + '\n' + pep + '\n' )
				out2.write( "X" + str( idx ) + '\t' + "X" + str( idx ) + '\n' )


# --- functions to handle DNA (genome sequence) input --- #

def get_subject_overlap( existing_parts, parts ):
	"""! @brief overlap calculator """
	
	start, end = int( parts[6] ), int( parts[7] )
	for each in existing_parts:
		if each['chr'] == parts[1]:
			if each['sstart'] < end:
				if each['send'] > start:
					return True	#overlap with existing HSPs detected
	return False	#no overlap with existing HSPs


def adjust_input_file( input_file, output_file, name_mapping_table_file ):
	"""! @brief generate clean input file """
	
	seqs = load_sequences( input_file )
	with open( output_file, "w" ) as out:
		with open( name_mapping_table_file, "w" ) as out2:
			for kidx, key in enumerate( seqs.keys() ):
				out.write( '>Z' + str( kidx+1 ) + '\n' + seqs[ key ].upper() + '\n' )
				out.write( '>Z' + str( kidx+1 ) + '_revcomp\n' + revcomp( seqs[ key ] ).upper() + '\n' )
				out2.write( key + '\t' + 'Z' + str( kidx+1 ) + '\tZ' + str( kidx+1 ) + '_revcomp\n'  )


def query_overlap( existing_parts, new_part, cutoff=3 ):
	"""! @brief check for query overlap with existing parts """
	
	overlap = 0
	for xeach in existing_parts:
		if xeach['qstart'] < new_part['qend']:
			if xeach['qend'] > new_part['qstart']:
				pos = sorted( [ xeach['qstart'], xeach['qend'], new_part['qend'], new_part['qstart'] ] )
				overlap = pos[2]-pos[1]
				if overlap > cutoff:
					return True		#significant overlap detected
	return False	#no significant overlap detected


def get_gene_groups( parts_per_query, max_gene_size ):
	"""! @brief build gene groups based on all BLAST hits per query """
	
	#gene: sstart, send, parts
	genes = {}
	for hit in parts_per_query:	#go through all HSPs
		best_gene = []	#list of potential genes to assign the HSP to
		for key in list(genes.keys()):
			gene = genes[ key ]
			if not query_overlap( gene['parts'], hit ):
				if hit['sstart'] < gene['send'] and hit['send'] > gene['sstart']:	#new part inside gene
					best_gene.append( { 'id': key, 'dist': 0 } )
				else:
					pos = sorted( [ hit['sstart'], gene['send'], hit['send'], gene['sstart'] ] )
					dist = pos[2]-pos[1]
					if dist <= max_gene_size:
						best_gene.append( { 'id': key, 'dist': dist } )
		if len( best_gene ) == 0:
			genes.update( { len( list(genes.keys()) )+1: { 	'chr': hit['chr'],
																			'sstart': hit['sstart'],
																			'send': hit['send'],
																			'parts': [ hit ]																			
																		} } )
		else:
			best_gene = sorted( best_gene, key=itemgetter('dist') )[0]['id']
			genes[ best_gene ]['parts'].append( hit )
			genes[ best_gene ]['sstart'] = min( [ genes[ best_gene ]['sstart'], hit['sstart'] ] )
			genes[ best_gene ]['send'] = min( [ genes[ best_gene ]['send'], hit['send'] ] )
	return genes


def load_blast_result( blast_result_file, max_gene_size ):
	"""! @brief load BLAST results """
	
	# --- load data per query --- #
	data = {}
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if int( parts[9] ) > int( parts[8] ):	#only consider fw hits
				try:
					if not get_subject_overlap( data[ parts[0] ], parts ):	#include only hits at new places
						data[ parts[0] ].append( { 'chr': parts[1],
																'qstart': int( parts[6] ),
																'qend': int( parts[7] ),
																'sstart': int( parts[8] ),
																'send': int( parts[9] ),
																'score': float( parts[-1] )
															} )
				except KeyError:
					data.update( { parts[0]: [ { 	'chr': parts[1],
																'qstart': int( parts[6] ),
																'qend': int( parts[7] ),
																'sstart': int( parts[8] ),
																'send': int( parts[9] ),
																'score': float( parts[-1] )
															} ] } )
			line = f.readline()
	
	# --- process data per query --- #
	gene_groups_per_query = {}
	for query in list(data.keys()):
		gene_groups_per_query.update( { query: get_gene_groups( data[ query ], max_gene_size ) } )
	return gene_groups_per_query


def extend_to_start( sstart, qstart, send, seq ):
	"""! @brief extend sequence to putative start M """
	
	core_seq = seq[ sstart-1:send ]
	
	# --- get sequence up to homolog M --- #
	# include complete homologous sequence if no stop codon is present
	minimal_extension = []
	start = sstart-4
	while start >= sstart-1 - 3*(qstart-1):
		if seq[ start:start+3 ] not in [ "TAA", "TAG", "TGA" ]:
			if not 'N' in seq[ start:start+3 ]:
				minimal_extension.append( seq[ start:start+3 ] )
			else:
				return "".join( minimal_extension[::-1] ) + core_seq
		else:
			return "".join( minimal_extension[::-1] ) + core_seq
		start -= 3
	
	# --- try to get longer sequence --- #
	# continue N-terminal part until M or stop codon is reached
	# fall back to M in minimal extension if not M reached before stop codon
	additional_extension = []
	while start >= 0:
		if seq[ start:start+3 ] not in [ "TAA", "TAG", "TGA" ]:
			if not 'N' in seq[ start:start+3 ]:
				additional_extension.append( seq[ start:start+3 ] )
			else:
				if "ATG" in additional_extension:
					return "ATG" + "".join( additional_extension[:additional_extension.index('ATG') ] ) + "".join( minimal_extension[::-1] ) + core_seq
				else:
					if "ATG" in minimal_extension:
						"ATG" + "".join( minimal_extension[:minimal_extension.index('ATG') ] ) + core_seq
					else:
						"".join( minimal_extension[::-1] ) + core_seq
		else:
			if "ATG" in additional_extension:
				return "ATG" + "".join( additional_extension[:additional_extension.index('ATG') ] ) + "".join( minimal_extension[::-1] ) + core_seq
			else:
				if "ATG" in minimal_extension:
					"ATG" + "".join( minimal_extension[:minimal_extension.index('ATG') ] ) + core_seq
				else:
					"".join( minimal_extension[::-1] ) + core_seq
		start -= 3
	return "".join( minimal_extension[::-1] ) + core_seq


def extend_to_stop( sstart, send, qlen, qend, seq ):
	"""! @brief extent to stop cdon """
	
	core_seq = seq[ sstart-1:send ]
	
	# --- minimal extension --- #
	minimal_extension = []
	start = send
	while start+3 < len( seq )-3:
		if seq[ start:start+3 ] not in [ "TAA", "TAG", "TGA" ]:
			if not 'N' in seq[ start:start+3 ]:
				minimal_extension.append( seq[ start:start+3 ] )
			else:
				return core_seq + "".join( minimal_extension ) + seq[start:start+3]
		else:
			return core_seq + "".join( minimal_extension )+ seq[start:start+3]
		start += 3


def get_seq_to_stop( start, seq ):
	"""! @brief get additional sequence to stop codon """
	
	extension = []
	while start < len( seq )-3:
		if seq[ start:start+3 ] not in [ "TAA", "TAG", "TGA" ]:
			if not 'N' in seq[ start:start+3 ]:
				extension.append( seq[ start:start+3 ] )
			else:
				return "".join( extension ) + seq[start:start+3]
		else:
			return "".join( extension ) + seq[start:start+3]
		start += 3


def get_from_start_to_stop( sstart, qstart, send, qlen, qend, seq ):
	"""! @brief get from start to stop """
	
	core_seq = seq[ sstart-1:send ]
	
	# --- get sequence up to homolog M --- #
	# include complete homologous sequence if no stop codon is present
	minimal_extension = []
	start = sstart-4
	while start >= sstart-1 - 3*(qstart-1):
		if seq[ start:start+3 ] not in [ "TAA", "TAG", "TGA" ]:
			if not 'N' in seq[ start:start+3 ]:
				minimal_extension.append( seq[ start:start+3 ] )
			else:
				return "".join( minimal_extension[::-1] ) + core_seq + get_seq_to_stop( send, seq )
		else:
			return "".join( minimal_extension[::-1] ) + core_seq + get_seq_to_stop( send, seq )
		start -= 3
	
	# --- try to get longer sequence --- #
	# continue N-terminal part until M or stop codon is reached
	# fall back to M in minimal extension if not M reached before stop codon
	additional_extension = []
	while start >= 0:
		if seq[ start:start+3 ] not in [ "TAA", "TAG", "TGA" ]:
			if not 'N' in seq[ start:start+3 ]:
				additional_extension.append( seq[ start:start+3 ] )
			else:
				if "ATG" in additional_extension:
					return "ATG" + "".join( additional_extension[:additional_extension.index('ATG') ] ) + "".join( minimal_extension[::-1] ) + core_seq + get_seq_to_stop( send, seq )
				else:
					if "ATG" in minimal_extension:
						"ATG" + "".join( minimal_extension[:minimal_extension.index('ATG') ] ) + core_seq
					else:
						"".join( minimal_extension[::-1] ) + core_seq
		else:
			if "ATG" in additional_extension:
				return "ATG" + "".join( additional_extension[:additional_extension.index('ATG') ] ) + "".join( minimal_extension[::-1] ) + core_seq + get_seq_to_stop( send, seq )
			else:
				if "ATG" in minimal_extension:
					"ATG" + "".join( minimal_extension[:minimal_extension.index('ATG') ] ) + core_seq + get_seq_to_stop( send, seq )
				else:
					"".join( minimal_extension[::-1] ) + core_seq + get_seq_to_stop( send, seq )
		start -= 3
	return "".join( minimal_extension[::-1] ) + core_seq + get_seq_to_stop( send, seq )


def calculate_alignment_sim( candidate, reference ):
	"""! @brief calculate alignment similarity (identical residues per ref) """
	
	counter = 0.0
	for z, aa in enumerate( candidate ):
		if aa != "-" and aa == reference[ z ]:
			counter += 1
	return counter / len( reference )


def find_missing_seq_between_hits( amino_acids, DNA, genetic_code, tmp_folder, mafft ):
	"""! @brief find missing sequence between BLAST hits """
	
	#get all blocks till GT splice site
	all_GT_starts = [ m.start() for m in re.finditer( "GT", DNA ) ]
	prev_GT_blocks = []
	for GTidx in all_GT_starts:
		if len( DNA[:GTidx] ) > 0:
			prev_GT_blocks.append( DNA[:GTidx] )
	
	#get all blocks following AG splice site
	all_AG_starts = [ m.start() for m in re.finditer( "AG", DNA ) ]
	following_AG_blocks = []
	for ATidx in all_AG_starts:
		if len( DNA[ATidx+2:] ) > 0:
			prev_GT_blocks.append( DNA[ATidx+2:] )
	
	# --- generate all potential missing sequences --- #
	all_seq_combis = []
	for x, GT in enumerate( prev_GT_blocks ):
		for y, AG in enumerate( following_AG_blocks ):
			out_fasta_file = tmp_folder + str( x ) + "_" + str( y ) + ".fasta"
			with open( out_fasta_file, "w" ) as out:
				out.write( '>refaa\n' + amino_acids + '\ncandidate\n' + GT+AG + '\n' )
			aln_file = out_fasta_file + ".aln"
			p = subprocess.Popen( args= mafft + " " + out_fasta_file + " > " + aln_file + " 2> " + aln_file+".err", shell=True )
			p.communicate()
			alignment = load_sequences( aln_file )
			sim = calculate_alignment_sim( alignment['candidate'], alignment['refaa'] )
			all_seq_combis.append( { 'seq': GT + AG, 'lendiff': 99999999- abs( len( GT + AG ) / 3.0 - amino_acids ), 'sim': sim } )
			
	try:
		return sorted( all_seq_combis, keys=itemgetter('sim', 'lendiff') )[-1]['seq']
	except:
		return ""


def construct_pep_seq( parts, seq, genetic_code, qlen, qpep, tmp_folder, mafft ):
	"""! @brief construct peptide sequence based on parts """
	
	parts = sorted( parts, key=itemgetter('sstart') )
	seq_parts = []
	for idx, part in enumerate( parts ):
		#first coding exon
		if idx == 0:
			#>1 coding exons
			if len( parts ) > 1:
				#complete
				if part['qstart'] == 1:
					seq_parts.append( seq[ part['sstart']-1:part['send'] ] )
				#start is missing
				else:
					if part['sstart'] > 3*part['qstart']:	#sstart reveales contig length left of hit
						seq_parts.append( extend_to_start( part['sstart'], part['qstart'], part['send'], seq ) )
						#seq_parts.append( seq[ part['sstart']-1 - 3*(part['qstart']-1):part['send'] ] )
						#reduce qstart by one due to 1-offset
						#include additional codons to enable detection of putative start M
					else:	#contig/chromosome end causes truncated pep sequence
						sys.stdout.write( "WARNING: contig end causes truncated peptide sequence (N-terminal end)\n")
						sys.stdout.flush()
						seq_parts.append( seq[ part['sstart']-1:part['send'] ] )
			#single coding exon gene
			else:
				#exon detected completely
				if part['qstart'] == 1 and part['qend'] == qlen:	#complete
					seq_parts.append( seq[ part['sstart']-1:part['send']+3 ] )
				#start and/or end is missing
				else:
					#sufficient upstream DNA sequence
					if part['sstart'] > 3*part['qstart']:	#sstart reveales contig length left of hit
						#sufficient downstream DNA sequence
						if len( seq )-part['send'] > 3*( qlen-part['qend'] ):
							seq_parts.append( get_from_start_to_stop( part['sstart'], part['qstart'], part['send'], qlen, part['qend'], seq ) )
							#seq_parts.append( seq[ part['sstart']-1 - 3*(part['qstart']-1):part['send']+3*( qlen-part['qend'] ) ] )	#NEW FUNCTION
						#insufficient downstream DNA sequence
						else:
							sys.stdout.write( "WARNING: contig end causes truncated peptide sequence (C-terminal end)\n" )
							sys.stdout.flush()
							seq_parts.append( extend_to_start( part['sstart'], part['qstart'], part['send'], seq ) )
					#insufficient upstream sequence
					else:
						#sufficient downstream DNA sequence
						if len( seq )-part['send'] > 3*( qlen-part['qend'] ):
							seq_parts.append( extend_to_stop( part['sstart'], part['send'], qlen, part['qend'], seq ) )
							#seq_parts.append( seq[ part['sstart']-1:part['send']+3*( qlen-part['qend'] ) ] )
							sys.stdout.write( "WARNING: contig end causes truncated peptide sequence (N-terminal end)\n" )
							sys.stdout.flush()
						#insufficient downstream DNA sequence
						else:
							sys.stdout.write( "WARNING: contig end causes truncated peptide sequence (N- and C-terminal end)\n" )
							sys.stdout.flush()
							seq_parts.append( seq[ part['sstart']-1:part['send'] ] )
		#last coding exon (more than 1 coding exons)
		elif idx == len( parts )-1:
			#direct hit next to previous one
			if parts[idx-1]['qend']+1 == part['qstart']:
				#seq_parts.append( "NNN" )
				#alignment until end
				if part['qend'] == qlen:
					seq_parts.append( seq[ part['sstart']-1:part['send']+3 ] )
				#alignment not complete at end
				else:
					#sufficient downstream DNA sequence
					if len( seq )-part['send'] > 3*( qlen-part['qend'] ):
						seq_parts.append( extend_to_stop( part['sstart'], part['send'], qlen, part['qend'], seq ) )
					#insufficient downstream DNA sequence
					else:
						sys.stdout.write( "WARNING: contig end causes truncated peptide sequence (C-terminal end)\n" )
						sys.stdout.flush()
						seq_parts.append( seq[ part['sstart']-1:part['send'] ] )
			#gap to previous coding exon hit
			elif parts[idx-1]['qend']+1 < part['qstart']:
				#seq_parts.append( "NNNNNN" )
				missing_seq = find_missing_seq_between_hits( 	qpep[ parts[idx-1]['qend']+1:part['qstart']+1 ],	#missing amino acids
																									seq[ parts[idx-1]['send']: part['sstart'] ],	#corresponding DNA sequence																									
																									genetic_code, tmp_folder, mafft
																								)
				#alignment until end
				if part['qend'] == qlen:
					seq_parts.append( missing_seq + seq[ part['sstart']-1:part['send']+3 ] )
				#alignment not complete at end
				else:
					#sufficient downstream DNA sequence
					if len( seq )-part['send'] > 3*( qlen-part['qend'] ):
						seq_parts.append( missing_seq + extend_to_stop( part['sstart'], part['send'], qlen, part['qend'], seq ) )
					#insufficient downstream DNA sequence
					else:
						sys.stdout.write( "WARNING: contig end causes truncated peptide sequence (C-terminal end)\n" )
						sys.stdout.flush()
						seq_parts.append( missing_seq + seq[ part['sstart']-1:part['send'] ] )
			#overlap with previous coding exon
			else:
				#seq_parts.append( "NNNNNNNNN" )
				#alignment until end
				if part['qend'] == qlen:
					seq_parts.append( seq[ part['sstart']-1+3*(1+parts[idx-1]['qend']-part['qstart']):part['send']+3 ] )
				#alignment not complete at end
				else:
					#sufficient downstream DNA sequence
					if len( seq )-part['send'] > 3*( qlen-part['qend'] ):
						seq_parts.append( extend_to_stop( part['sstart']+3*(1+parts[idx-1]['qend']-part['qstart']), part['send'], qlen, part['qend'], seq ) )
						#seq_parts.append( seq[ part['sstart']-1+3*(1+parts[idx-1]['qend']-part['qstart']):part['send']+3*( qlen-part['qend'] ) ] )
					#insufficient downstream DNA sequence
					else:
						sys.stdout.write( "WARNING: contig end causes truncated peptide sequence (C-terminal end)\n" )
						sys.stdout.flush()
						seq_parts.append( seq[ part['sstart']-1+3*(1+parts[idx-1]['qend']-part['qstart']):part['send'] ] )
		#all enclosed coding exons
		else:
			#direct hit next to previous one
			if parts[idx-1]['qend']+1 == part['qstart']:
				#seq_parts.append( "NNN" )
				### ADD SPLICE SITE CHECK
				seq_parts.append( seq[ part['sstart']-1:part['send'] ] )
			
			#gap
			elif parts[idx-1]['qend']+1 < part['qstart']:
				#seq_parts.append( "NNNNNN" )
				# --- search for missing sequence --- #
				missing_seq = find_missing_seq_between_hits( 	qpep[ parts[idx-1]['qend']+1:part['qstart']+1 ],	#missing amino acids
																									seq[ parts[idx-1]['send']: part['sstart'] ],	#corresponding DNA sequence																									
																									genetic_code, tmp_folder, mafft
																								)		
				seq_parts.append( missing_seq + seq[ part['sstart']-1:part['send'] ] )
			
			#overlap
			else:
				#seq_parts.append( "NNNNNNNNN" )
				seq_parts.append( seq[ part['sstart']-1+3*(1+parts[idx-1]['qend']-part['qstart']):part['send'] ] )
			
	#find start M
	#find stop codon
	lengths = []
	for each in seq_parts:
		lengths.append( len( each ) / 3.0 )
	return seq_parts


def reconstruct_sequences( gene_groups_per_query, genome_seq, querys, tmp_folder, mafft ):
	"""! @brief reconstruct peptide sequences """
	
	genetic_code = {	'CTT': 'L', 'ATG': 'M', 'AAG': 'K', 'AAA': 'K', 'ATC': 'I',
								'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T',
								'AGC': 'S', 'ACA': 'T', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N',
								'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H',
								'ACG': 'T', 'CCG': 'P', 'AGT': 'S', 'CAG': 'Q', 'CAA': 'Q',
								'CCC': 'P', 'TAG': '*', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C',
								'CGA': 'R', 'CCA': 'P', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R',
								'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G',
								'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'GAG': 'E', 'TCG': 'S',
								'TTA': 'L', 'GAC': 'D', 'TCC': 'S', 'GAA': 'E', 'TCA': 'S',
								'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A',
								'GTG': 'V', 'TTC': 'F', 'GTT': 'V', 'GCT': 'A', 'ACC': 'T',
								'TTG': 'L', 'CGT': 'R', 'TAA': '*', 'CGC': 'R'
							}
	
	all_pep_seqs = {}
	for query in list(gene_groups_per_query.keys()):
		genes = gene_groups_per_query[ query ]
		pep_seqs = {}
		for gene in list(genes.keys()):
			chromosome = genes[ gene ]['chr']
			parts = genes[ gene ]['parts']
			pep_seq = construct_pep_seq( parts, genome_seq[ chromosome ], genetic_code, len( querys[ query ] ), querys[ query ], tmp_folder, mafft )
			pep_seqs.update( { gene: translate( "".join( pep_seq ), genetic_code ) } )	#TRANSLATION
		all_pep_seqs.update( { query: pep_seqs } )
	return all_pep_seqs


def dna_screener( subject, peptide_file, subject_name_file, dna_folder, query_file, mafft, cpus, max_gene_size, makeblastdb, tblastn ):
	"""! @brief screen genome sequence for encoded peptides matching the query """
	
	if not os.path.exists( dna_folder ):
		os.makedirs( dna_folder )
	
	### all input sequences are also included as revcomps >> only forward hits are processed ###
	dna_file = dna_folder + "subject_dna.fasta"
	name_mapping_table_file = dna_folder + "subject_name_mapping.txt"
	if not os.path.isfile( dna_file ):
		adjust_input_file( subject, dna_file, name_mapping_table_file )
	
	### run BLAST search ###
	blast_db = dna_folder + "db"
	blast_result_file = dna_folder + "results.txt"
	
	if not os.path.isfile( blast_result_file ):
		p = subprocess.Popen( args= makeblastdb + " -in " + dna_file + " -out " + blast_db + " -dbtype nucl", shell=True )
		p.communicate()
		p = subprocess.Popen( args= tblastn + " -query " + query_file + " -db " + blast_db + " -out " + blast_result_file + " -outfmt 6 -evalue 0.01 -word_size 4 -num_threads " + str( cpus ), shell=True )
		p.communicate()
		
	gene_groups_per_query = load_blast_result( blast_result_file, max_gene_size )
	
	genome_seq = load_sequences( dna_file )
	query_sequences = load_sequences( query_file )
	
	tmp_folder = dna_folder + "tmp/"
	if not os.path.exists( tmp_folder ):
		os.makedirs( tmp_folder )
	pep_seqs = reconstruct_sequences( gene_groups_per_query, genome_seq, query_sequences, tmp_folder, mafft )
	
	no_dup_status = True	#this option could be integrated later
	black_list = {}
	
	with open( subject_name_file, "w" ) as out2:
		with open( peptide_file, "w" ) as out:	#this file will serve as input for the remaining pipeline
			for query in list(pep_seqs.keys()):
				for gene in list(pep_seqs[ query ].keys()):
					if no_dup_status:
						try:
							black_list[ pep_seqs[ query ][ gene ] ]
						except KeyError:
							out.write( '>' + str( query ) + "_%_" + str( gene ) + '\n' + pep_seqs[ query ][ gene ] + '\n' )
							out2.write( str( query ) + "_%_" + str( gene ) + "\t" + str( query ) + "_%_" + str( gene ) + '\n' )
							black_list.update( { pep_seqs[ query ][ gene ]: None } )
					else:
						out.write( '>' + str( query ) + "_%_" + str( gene ) + '\n' + pep_seqs[ query ][ gene ] + '\n' )
						out2.write( str( query ) + "_%_" + str( gene ) + "\t" + str( query ) + "_%_" + str( gene ) + '\n' )


### --- final summary section --- ###

def generate_final_pep_files( 	peps, final_pep_folder, candidates_by_gene,
													xsimcut, xconsrescut, xconsregcut,
													sim_per_pep, cons_pos_per_pep, cons_reg_per_pep,
													summary_file, subject_name_mapping_table
												):
	"""! @brief generate multiple FASTA files per gene for all candidates passing the filter """
	
	complete_summary = []
	final_file_per_gene = {}
	with open( summary_file, "w") as summary:
		summary.write( "ID\tGene\tSimilarity\tConservedResidues\tConservedRegions\n" )
		for gene in list(candidates_by_gene.keys()):
			output_file = final_pep_folder + gene + ".fasta"
			final_file_per_gene.update( { gene: output_file } )
			with open( output_file, "w" ) as out:
				candidates = candidates_by_gene[ gene ]
				values_for_sorting = []
				for candidate in candidates:
					try:
						psim = sim_per_pep[ candidate ][ gene ]
					except KeyError:
						psim = 0
					try:
						pres = cons_pos_per_pep[ candidate ][ gene ]
					except KeyError:
						pres = 0
					try:
						preg = cons_reg_per_pep[ candidate ][ gene ]
					except KeyError:
						preg = 0
					
					values_for_sorting.append( { 	'id': candidate,
																		'sim': psim,
																		'res': pres,
																		'reg': preg,
																		'gene': gene,
																		'label': subject_name_mapping_table[ candidate ]
																	} )
				values_for_sorting = sorted( values_for_sorting, key=itemgetter('res', 'reg', 'sim') )[::-1]
				complete_summary += values_for_sorting
				for y, candidate in enumerate( values_for_sorting ):
					if candidate['sim'] >= xsimcut and candidate['res'] >= xconsrescut and candidate['reg'] >= xconsregcut:
						out.write( '>' + subject_name_mapping_table[ candidate['id'] ] + '\n' + peps[ candidate['id'] ] + '\n'  )
						summary.write( "\t".join( map( str, [ 	subject_name_mapping_table[ candidate['id'] ],
																						gene + "_" + str( y+1 ),
																						candidate['sim'],
																						candidate['res'],
																						candidate['reg']
																					] ) ) + '\n' )
	return complete_summary, final_file_per_gene


def which( cmd, mode=os.F_OK, path=None ):
	"""! @brief checks if software is available
	code is based on whichcraft (Daniel Roy Greenfield): https://github.com/cookiecutter/whichcraft/blob/master/whichcraft.py#L20
	"""
	
	def _access_check( fn, mode ):
		return os.path.exists(fn) and os.access(fn, mode) and not os.path.isdir(fn)
	
	if os.path.dirname( cmd ):
		if _access_check( cmd, mode ):
			return cmd
		return None
	if path is None:
		path = os.environ.get( "PATH", os.defpath )
	if not path:
		return None
	
	path = path.split( os.pathsep )
	files = [ cmd ]
	seen = set()
	for directory in path:
		normdir = os.path.normcase( directory )
		if normdir not in seen:
			seen.add( normdir )
			for thefile in files:
				name = os.path.join( directory, thefile )
				if _access_check( name, mode ):
					return name
	return None


def validate_input( pos_data_dir, bait_seq_data_dir, makeblastdb, blastp, tblastn, mafft, treemethod, fasttree, iqtree ):
	"""! @brief validate input """
	
	### tool dependency checks ###
	tool_labels = [ "makeblastdb", "blastp", "tblastn", "mafft", "fasttree" ]
	tool_state = True
	for t, each_tool in enumerate( [ makeblastdb, blastp, tblastn, mafft, fasttree, iqtree ] ):
		each_tool_state = which( each_tool )
		if each_tool_state == None:
			if len( each_tool ) > 0:
				sys.stdout.write("ERROR: specified binary of " + tool_labels[ t ] + " not detected - " + each_tool + "\n")
				sys.stdout.flush()
				tool_state = False
	if not tool_state:
		sys.exit( "ERROR: specified binaries not detected. Please see error messages above for details." )
	
	
	### data consistency checks ###
	pos_data_files = glob.glob( pos_data_dir + "*.txt" ) + glob.glob( pos_data_dir + "*.res" )
	fasta_files = glob.glob( bait_seq_data_dir + "*.fa" ) + glob.glob( bait_seq_data_dir + "*.fasta" )
	baits = {}
	for fasta in fasta_files:
		baits.update( { fasta.split('/')[-1].split('.')[0]: load_sequences( fasta ) } )
	errors = []
	for filename in pos_data_files:
		ID = filename.split('/')[-1].split('.')[0]
		seq_id = False
		len_error_state = False
		try:
			with open( filename, "r" ) as f:
				seq_id = ""
				positions = []
				line = f.readline()
				while line:
					if line[0] == "!":
						seq_id = line.strip()[1:]
					else:
						parts = line.strip().split('\t')
						if parts[0] == "R":
							positions.append( int( parts[2] ) )
						elif parts[0] == "D":
							positions.append( int( parts[3] ) )
					line = f.readline()
			length = len( baits[ ID ][ seq_id ] )
			if len( positions ) > 0:
				if max( positions )-1 > length:
					len_error_state = True
					errors.append( { 'seq': seq_id, 'gene': ID, 'len': len_error_state } )
		except:
			errors.append( { 'seq': seq_id, 'gene': ID, 'len': len_error_state } )
	return errors


def get_pathway( pathway_file, candidates_by_gene ):
	"""! @brief load order of genes from given pathway file """
	
	# --- load file content if any --- #
	if os.path.isfile( pathway_file ):
		with open( pathway_file, "r" ) as f:
			pathway = f.read().strip().split('\n')
	else:
		pathway = []
	
	# --- check for existence of genes and completeness of gene list --- #
	final_pathway = []
	for gene in pathway:
		try:
			candidates_by_gene[ gene ]
			final_pathway.append( gene )
		except KeyError:
			sys.stdout.write( "ERROR: gene specified in pathway file was not investigated (baits missing) - " + gene + "\n" )
			sys.stdout.flush()
	
	for gene in list(candidates_by_gene.keys()):
		if gene not in final_pathway:
			final_pathway.append( gene )
	return final_pathway


def generate_summary_html( html_file, summary, peps, cons_pos_per_pep_extra, pathway, pos_data_per_gene ):
	"""! @brief generate final summary HTML """
	
	#cons_pos_per_pep_extra = { 'pep1': [ A, T, A, S, -, T, S, -, -, - ], 'pep2': [ A, S, A, S, T, -, - ] }	contains residues at key positions
	#summary = [ { 'id': candidate, 'sim': psim, 'res': pres, 'reg': preg, 'gene': gene, 'label': subject_name_mapping_table[ candidate ] }, {...}, ... ]
	#supply pathway file for order of genes in output => additional input file required
	
	# --- preparing data and final filtering --- #
	entry_order = []
	status_per_gene = {}
	for gene in pathway:
		for entry in summary:
			if entry['gene'] == gene:
				try:
					status_per_gene[ gene ]
					if entry['res'] == 100.0:
						entry_order.append( entry )
				except KeyError:
					entry_order.append( entry )
					status_per_gene.update( { gene: True } )	#only take perfect hits (or best hit if no perfect hit available)
	
	# --- generate HTML file --- #
	#Gene		SeqID	Similarity	Residues		LinkToData
	with open( html_file, "w" ) as out:
		out.write( "<h1>SUMMARY</h1>\n")
		out.write( "These results were produced by KIPEs " + __version__ + ". " )
		out.write('Please see the instructions in the <a href="https://github.com/bpucker/KIPEs">KIPEs</a> github repository for details.\n' )
		out.write( "<table>\n<tr> <th>Gene</th><th>SeqID</th><th>Similarity</th><th>Residues</th><th>PeptideSequence</th> </tr>" )
		for entry in entry_order:
			new_line = []
			new_line.append( "<td>" + entry['gene'] + "</td>" )	#gene name (function in pathway)
			new_line.append( "<td>" + entry['label'] + "</td>" )	#sequence ID
			new_line.append( "<td>" + str( round( entry['sim'], 3 ) ) + "%</td>" )	#sequence ID
			
			try:
				residues_per_gene = pos_data_per_gene[ entry['gene'] ]['residues']	#[ { 'aa': X, 'pos': int( res[1:] ) }, { 'aa': X, 'pos': int( res[1:] ) }, ...]
				present_residues = cons_pos_per_pep_extra[ entry['id']  ][ entry['gene']]
				tmp_blocks = []
				for zz, residue in enumerate( present_residues ):
					if residue in residues_per_gene[ zz ]['aa']:	#match (conserved residue)
						tmp_blocks.append( "/".join( residues_per_gene[ zz ]['aa'] )  + str( residues_per_gene[ zz ]['pos']  ) + residue )
					else:	#missmatch
						tmp_blocks.append( '<strong style="color: red;">' + "/".join( residues_per_gene[ zz ]['aa'] )  + str( residues_per_gene[ zz ]['pos']  ) + residue + "</strong>" )
				new_line.append( "<td>" + ", ".join( tmp_blocks ) + "</td>" )
			except KeyError:
				new_line.append( "<td>.</td>" )
			
			chunks = [ peps[ entry['id'] ][i:i+100] for i in range( 0, len( peps[ entry['id'] ] ), 100 ) ]
			new_line.append( "<td>" + "<br>".join( chunks ) + "<br></td>" )	#add link to final peptide files
			
			
			out.write( "<tr>" + "".join( new_line ) + "</tr>\n" )
		
		out.write( "</table>\n" )


### ---- FORESTER ###


def modify_FASTA_names( input_file, output_file, mapping_table_file ):
	"""! @brief modify names in FASTA file """
	
	seqs = load_sequences( input_file )
	mapping_table = {}
	with open( output_file, "w" ) as out:
		for idx, key in enumerate( seqs.keys() ):
			out.write( '>seq' + str( idx ).zfill(5) + '\n' + seqs[ key ] + '\n' )
			mapping_table.update( { key: 'seq' + str( idx ).zfill(5) } )
	
	with open( mapping_table_file, "w" ) as out:
		for key in list(mapping_table.keys()):
			if len( key ) == 0:
				sys.exit( "ERROR: " + seqs[ key ] )
			out.write( key + '\t' + mapping_table[ key ] + '\n' )


def load_mapping_table( mapping_table_file ):
	"""! @brief load mapping table """
	
	mapping_table = {}
	with open( mapping_table_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			mapping_table.update( { parts[1]: parts[0] } )
			line = f.readline()
	return mapping_table


def modify_names_in_tree( input_tree_file, output_tree_file, mapping_table_file ):
	"""! @brief replace short names in tree by original names """
	
	mapping_table = load_mapping_table( mapping_table_file )
	
	with open( input_tree_file, "r" ) as f:
		tree = f.read()
	
	for key in list(mapping_table.keys()):
		tree = tree.replace( key, mapping_table[ key ] )
	
	with open( output_tree_file, "w" ) as out:
		out.write( tree )


def construct_tree( input_file,  output_dir, mafft, treemethod, fasttree, iqtree, occupancy, name ):
	"""! @brief handle tree construction """
		
	if output_dir[-1] != '/':
		output_dir += "/"
	
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )
	
	mod_FASTA = output_dir + "names_modified.fasta"
	mapping_table_file = output_dir + "seq_names_mapping_table.txt"
	modify_FASTA_names( input_file, mod_FASTA, mapping_table_file )
	
	alignment_file = mod_FASTA + ".aln"
	p = subprocess.Popen( args= " ".join( [ mafft, mod_FASTA, ">", alignment_file, "2>", alignment_file+".log" ] ), shell=True )
	p.communicate()
	
	clean_alignment_file = alignment_file + ".cln"
	alignment_trimming( alignment_file, clean_alignment_file, occupancy )
	
	tree_file = clean_alignment_file + ".tre"
	if treemethod == "fasttree":		#construct tree with FastTree
		p = subprocess.Popen( args= " ".join( [ fasttree, "-wag -nosupport <", clean_alignment_file, ">", tree_file, "2>", tree_file+".log" ] ), shell=True )
		p.communicate()
	elif treemethod == "iqtree":	#construct tree with IQ-TREE
		tmp_tree_file = clean_alignment_file + ".treefile"	#".treefile" is appended to provided alignment file name
		p = subprocess.Popen( args= iqtree + " -nt 1 -alrt 1000 -bb 1000 -s " + clean_alignment_file, shell=True )
		p.communicate()
		if not os.path.isfile( tmp_tree_file ):	#run analysis without bootstrapping if there are not enough sequences
			p = subprocess.Popen( args= iqtree + " -nt 1 -s " + clean_alignment_file, shell=True )	#bootstrapping: -alrt 1000 -bb 1000 
			p.communicate()
		p = subprocess.Popen( args= "mv " + tmp_tree_file + " " + tree_file, shell=True )
		p.communicate()
	
	output_tree_file = output_dir + name + "FINAL_TREE.tre"
	modify_names_in_tree( tree_file, output_tree_file, mapping_table_file )
	
	return output_tree_file


def forester( input_folder, output_folder, refseq_folder, mafft, treemethod, fasttree, iqtree, occupancy ):
	"""! @brief run everything """
	
	clean = False	#temporary files are stored
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	input_files = glob.glob( input_folder + "*.fasta" )
	
	for m, filename in enumerate( input_files ):
		ID = filename.split('/')[-1].replace( ".fasta", "" )
		sys.stdout.write( "processing: " +  ID + " - (" + str( m+1 ) + "/" + str( len( input_files ) ) + ")\n" )
		sys.stdout.flush()
		ref_seq_file = refseq_folder + ID + ".fasta"
		tree_input_file = output_folder + ID + ".fasta"
		if os.path.isfile( ref_seq_file ):
			p = subprocess.Popen( args= "cat " + filename + " " + ref_seq_file + " > " + tree_input_file, shell=True )
			p.communicate()
			tree_tmp_folder = output_folder + ID + "_tmp/"

			tree = construct_tree( tree_input_file, tree_tmp_folder, mafft, treemethod, fasttree, iqtree, occupancy, ID )
			
			p = subprocess.Popen( args= "cp " + tree + " " + output_folder + ID + ".tre", shell=True )
			p.communicate()
			if clean:
				p = subprocess.Popen( args= "rm -r " + tree_tmp_folder, shell=True )	#remove temp folders
				p.communicate()



### --- COEXPRESSION --- ###

def load_expression_values( filename ):
	"""! @brief load all expression values """
	
	expression_data = {}
	with open( filename, "r" ) as f:
		tissues = f.readline().strip().split('\t')[1:]
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			expression = {}
			for idx, each in enumerate( parts[1:] ):
				expression.update( { tissues[  idx ] : float( parts[ idx+1 ] ) } )
			line = f.readline()
			expression_data.update( { parts[0]: expression } )
	return expression_data


def one_vs_all_coexp( candidate, candidate_genes, gene_expression, single_coexp_output_file, rcutoff, pcutoff, minexpression ):
	"""! @brief compare candidate gene expression against all genes to find co-expressed genes """
	
	tissues = sorted( gene_expression[ list(gene_expression.keys())[0] ].keys() )
	coexpressed_genes = {}
	output_data = []
	for i, gene2 in enumerate( candidate_genes ):
		if candidate != gene2:
			values = []
			total_expression = 0
			for tissue in tissues:
				try:
					x = gene_expression[ candidate ][ tissue ]
					y = gene_expression[ gene2 ][ tissue ]
					total_expression += y
					if not math.isnan( x ) and not math.isnan( y ) :
						values.append( [ x, y ] )
				except KeyError:
					sys.stdout.write( "ERROR (coexp): tissue not found - " + tissue+ "\n" )
					sys.stdout.flush()
			r, p = stats.spearmanr( values )
			if not math.isnan( r ) and total_expression > minexpression:
				if r > rcutoff and p < pcutoff:
					coexpressed_genes.update( { gene2: r } )
					output_data.append( { 'id': gene2, 'r': r, 'p': p } )
	
	# --- write data into output file --- #
	with open( single_coexp_output_file, "w" ) as out:
		for each in sorted( output_data, key=itemgetter('r', 'p') ):
			new_line = [ candidate, each['id'], str( each['r'] ), str( each['p'] ) ]
			out.write( "\t".join( new_line ) + "\n" )	
	
	return coexpressed_genes


def load_genes_of_interest( KIPEs_summary_file ):
	"""! @brief load genes per step in pathway """
	
	genes_per_function = {}
	with open( KIPEs_summary_file, "r" ) as f:
		f.readline()	#remove header
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			step = parts[1].split('_')[0]
			try:
				genes_per_function[ step ].append( parts[0] )
			except KeyError:
				genes_per_function.update( { step: [ parts[0] ] } )
			line = f.readline()
	return genes_per_function


def construct_html_output( html_file, genes_per_function, function_order, coexp_per_gene ):
	"""! @brief generate HTML heatmap """
	
	html_string = '<html><style>#rotator {transform: rotate(-90deg); height: 250px; width: 10px;} </style><table>'
	
	html_string += '<tr style="width:10px" ><th id="rotator">Function_GeneID</th>'
	for function in function_order:
		genes = sorted( genes_per_function[ function ] )
		for gene in genes:
			html_string += '<th id="rotator">'+function + "_" + gene+"</th>\n"
	html_string += "</tr>\n"
	
	for function1 in function_order:
		genes = sorted( genes_per_function[ function1 ] )
		for gene in genes:
			html_string += "<tr><td>" + function1 + "_" + gene + "</td>\n"
			for function2 in function_order:
				genes2 = sorted( genes_per_function[ function2 ] )
				for gene2 in genes2:
					try:
						html_string += '<td style="width:10px" bgcolor="' + calculate_color( coexp_per_gene[ gene ][ gene2 ] ) + '"><span title="' + function1 +"_"+ gene + " & " + function2 +"_"+ gene2 + '">' + str( round( coexp_per_gene[ gene ][ gene2 ], 3 ) ) + "</span></td>\n"
					except KeyError:
						html_string += '<td style="width:10px" bgcolor="#FFFFFF"><span title="' + function1 +"_"+ gene + " & " + function2 +"_"+ gene2 + '">0</span></td>\n'
			html_string += "</tr>\n"
			
	html_string += "</table></html>\n"

	with open( html_file, "w" ) as out:
		out.write( html_string )


def calculate_color( value ):
	"""! @brief calculates color for given value """
	
	r = 255
	g = 255-int( 255*value )
	b = 255-int( 255*value )
	
	color = '#%02x%02x%02x' % (r, g, b)
	
	return color


def automatic_coexp( expression_file, output_file, KIPEs_summary_file, coexp_output_folder, rcutoff, pcutoff, minexpression ):
	"""! @brief run automatic co-expression analysis for all candidates identified by KIPEs """

	gene_expression = load_expression_values( expression_file )
	genes_per_function = load_genes_of_interest( KIPEs_summary_file )
	
	candidate_genes = [ x for function in list(genes_per_function.values()) for x in function ]
	
	single_coexp_folder = coexp_output_folder + "individual_genes/"
	if not os.path.exists( single_coexp_folder ):
		os.makedirs( single_coexp_folder )
	
	coexp_per_gene = {}
	for gene in candidate_genes:
		single_coexp_output_file = single_coexp_folder + gene + ".txt"
		coexp_per_gene.update( { gene: one_vs_all_coexp( gene, candidate_genes, gene_expression, single_coexp_output_file, rcutoff, pcutoff, minexpression ) } )
	
	
	function_order = sorted( genes_per_function.keys() )	#steps in the pathway
	with open( output_file, "w" ) as out:
		out.write( "\t".join( [ "Function", "GeneID" ] + function_order ) + "\n" )
		for function1 in function_order:
			candidates1 = genes_per_function[ function1 ]
			data_per_candidate = []
			for candidate1 in candidates1:
				values_per_function = []
				for function2 in function_order:
					candidates2 = genes_per_function[ function2 ]
					values = []
					for candidate2 in candidates2:
						try:
							values.append( coexp_per_gene[ candidate1 ][ candidate2 ] )
						except KeyError:
							pass
					if len( values ) > 0:
						values_per_function.append( round( float( sum( values ) ) / len( values ), 3 ) )
					else:
						values_per_function.append( 0 )
				data_per_candidate.append( { 'ID': candidate1, 'val': float( sum( values_per_function ) ) / len( values_per_function ), 'rawval': values_per_function } )
			best_candidate = sorted( data_per_candidate, key=itemgetter('val') )[::-1]	[0]
			best_values = data_per_candidate[ candidates1.index( best_candidate['ID'] ) ]['rawval']
			out.write( "\t".join( [ function1, best_candidate['ID'] ] + list(map( str, best_values )) ) + "\n" )
	
	html_file = output_file + ".html"
	construct_html_output( html_file, genes_per_function, function_order, coexp_per_gene )


def md5_calculator( input_file ):
	"""! @brief calculate md5sum of given file """
	
	if len( input_file ) > 0:
		if os.path.isfile( input_file ):
			with open( input_file, "rb" ) as f:
				content = f.read()
			try:
				return hashlib.md5( content ).hexdigest()
			except NameError:
				return "n/a"
		else:
			return "n/a"
	else:
		return "n/a"


def write_general_input_to_doc_file( fulldoc, bait_seq_data_dir, output_dir, subject, pos_data_dir, seqtype, mafft, blastp, tblastn, makeblastdb, treemethod, fasttree, iqtree, pathway_file, cpus, score_ratio_cutoff, similarity_cutoff, max_gene_size, xsimcut, xconsrescut, xconsregcut, checks, possibility_cutoff, exp_file, forester_state ):
	"""! @brief write all input into the documentation file 
		@note This function is based on the MYB_annotator: https://doi.org/10.1101/2021.10.16.464636
	"""
	
	# --- KIPEs version and general info --- #
	fulldoc.write( "Please cite 'Pucker, B.; Reiher, F.; Schilbert, H.M. Automatic Identification of Players in the Flavonoid Biosynthesis with Application on the Biomedicinal Plant Croton tiglium. Plants 2020, 9, 1103. https://doi.org/10.3390/plants9091103 ' when using KIPEs3.py.\n\n" )
	fulldoc.write( "KIPEs3.py version: " + __version__ + "\n\n" )
	
	fulldoc.write( "This documentation file contains the input file names followed by their md5sums. Modification of the file content will result in a different md5sum. All settings are documented as well to enable reproduction.\n\n" )
	
	# --- bait sequences --- #
	fasta_files = sorted( glob.glob( bait_seq_data_dir + "*.fa" ) + glob.glob( bait_seq_data_dir + "*.fasta" ) )
	fulldoc.write( "bait sequence files:\n" )
	for xfasta in fasta_files:
		xfasta_md5 = md5_calculator( xfasta )
		fulldoc.write( xfasta + "\t" + xfasta_md5 + "\n" )
	fulldoc.write( "\n" )
	
	# --- output folder and input file --- #
	fulldoc.write( "output folder: "+ output_dir +"\n" )
	subject_md5 = md5_calculator( subject )
	fulldoc.write( "subject file: " + subject + "\t" + subject_md5 + "\n" )
	fulldoc.write( "input sequence type: " + seqtype + "\n" )
	pathway_file_md5 = md5_calculator( pathway_file )
	fulldoc.write( "pathway file (optional): " + pathway_file + "\t" + pathway_file_md5 + "\n" )
	exp_file_md5 = md5_calculator( exp_file )
	fulldoc.write( "expression file (optional): " + exp_file + "\t" + exp_file_md5 + "\n\n" )
	
	# --- residue information --- #
	pos_data_files = sorted( glob.glob( pos_data_dir + "*.txt" ) + glob.glob( pos_data_dir + "*.res" ) )
	fulldoc.write( "residue info files:\n" )
	for xpos in pos_data_files:
		xpos_md5 = md5_calculator( xpos )
		fulldoc.write( xpos + "\t" + xpos_md5 + "\n" )
	fulldoc.write( "\n" )
	
	# --- provided paths to tools --- #
	fulldoc.write( "MAFFT path: " + mafft + " " )
	try:
		mafft_version_raw = subprocess.Popen( args=mafft + " --version", stderr=subprocess.PIPE, shell=True )
		mafft_version = mafft_version_raw.stderr.read()
		fulldoc.write ( "(" + str( mafft_version )[2:-3] + ")\n" )	#remove characters introduced through binary
	except:
		fulldoc.write ( "(version detection failed)\n" )	#if no MAFFT installation was detected
	fulldoc.write( "BLASTp path: " + blastp + "\n" )
	fulldoc.write( "tBLASTn path: " + tblastn + "\n" )
	fulldoc.write( "makeblastdb path: " + makeblastdb + "\n" )
	fulldoc.write( "Tree method: " + treemethod + "\n" )
	try:
		fasttree_version_raw = subprocess.Popen( args=fasttree + " -help", stderr=subprocess.PIPE, shell=True )
		fasttree_version = fasttree_version_raw.stderr.read()
		fulldoc.write( "FastTree2 path: " + fasttree + " (" + str( fasttree_version )[10:18] + ")\n" )
	except:
		fulldoc.write( "FastTree2 path: " + fasttree + " (Please add version manually)\n" )
	try:
		iqtree_version_raw = subprocess.Popen( args=iqtree, stderr=subprocess.PIPE, shell=True )
		iqtree_version = iqtree_version_raw.stderr.read()
		fulldoc.write( "IQ-TREE path: " + iqtree + " (" + str( iqtree_version )[26:32] + ")\n" )
	except:
		fulldoc.write( "IQ-TREE path: " + iqtree + " (Please add version manually)\n\n" )
	
	# --- cutoffs and other settings --- #
	fulldoc.write( "--cpus: " + str( cpus ) + "\n" )
	fulldoc.write( "--scoreratio: " + str( score_ratio_cutoff ) + "\n" )
	fulldoc.write( "--simcut: " + str( similarity_cutoff ) + "\n" )
	fulldoc.write( "--genesize: " + str( max_gene_size ) + "\n" )
	fulldoc.write( "--minsim: " + str( xsimcut ) + "\n" )
	fulldoc.write( "--minres: " + str( xconsrescut ) + "\n" )
	fulldoc.write( "--minreg: " + str( xconsregcut ) + "\n" )
	fulldoc.write( "--possibilities: " + str( possibility_cutoff ) + "\n\n" )
	
	# --- options --- #
	fulldoc.write( "check input completeness: " + str( checks ) + "\n" )
	fulldoc.write( "construction of final gene trees: " + str( forester_state ) + "\n" )


def KIPEs( bait_seq_data_dir, output_dir, subject, pos_data_dir, seqtype, mafft, blastp, tblastn, makeblastdb, treemethod, fasttree, iqtree, pathway_file, cpus, score_ratio_cutoff, similarity_cutoff, max_gene_size, xsimcut, xconsrescut, xconsregcut, checks, possibility_cutoff, exp_file, forester_state, rcutoff, pcutoff, minexpression ):
	"""! @brief run whole KIPEs analysis for one subject sequence file """
	
	errors = validate_input( pos_data_dir, bait_seq_data_dir, makeblastdb, blastp, tblastn, mafft, treemethod, fasttree, iqtree, )
	if len( errors ) > 0:
		for error in errors:
			if error['seq']:
				if error['len']:
					sys.stdout.write( "ERROR: conserved residue position in sequence " + error['seq'] + " of " + error['gene'] + " exceeds sequence length.\n" )
					sys.stdout.flush()
				else:
					sys.stdout.write( "ERROR: conserved residue sequence " + error['seq'] + " of " + error['gene'] + " is not matching bait sequences.\n" )
					sys.stdout.flush()
			else:
				sys.stdout.write( "ERROR: no conserved residue information detected for " + error['gene'] + ".\n" )
				sys.stdout.flush()
		if checks == "on":
			sys.exit( "ERROR: Execution of script terminated due to errors. Fix errors or set '--checks' to 'off' in order to proceed." )
	
	if output_dir[-1] != "/":
		output_dir += "/"
	
	if len( pos_data_dir ) > 0:
		if pos_data_dir[-1] != "/":
			pos_data_dir += "/"
	
	if bait_seq_data_dir[-1] != "/":
		bait_seq_data_dir += "/"
		
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )
	
	if subject[-3:] in [ ".gz", ".GZ" ]:	#KIPEs3 can handle gzip compressed input files (subject only!)
		try:
			p = subprocess.Popen( args="cp " + subject + " " + output_dir, shell=True  )
			p.communicate()
			if "/" in subject:
				subject = output_dir + subject.split('/')[-1]
			else:
				subject = output_dir + subject
			p = subprocess.Popen( args="gunzip " + subject, shell=True  )
			p.communicate()
		except:
			sys.exit( "ERROR: subject file appears to be compressed, but gzip is not available." )
		subject = subject[:-3]
	
	# --- start documentation --- #
	documentation_file = output_dir + "documentation_of_parameters_and_inputs.txt"
	with open( documentation_file, "w" ) as fulldoc:
		# --- generate basic documentation of input --- #
		write_general_input_to_doc_file( fulldoc, bait_seq_data_dir, output_dir, subject, pos_data_dir, seqtype, mafft, blastp, tblastn, makeblastdb, treemethod, fasttree, iqtree, pathway_file, cpus, score_ratio_cutoff, similarity_cutoff, max_gene_size, xsimcut, xconsrescut, xconsregcut, checks, possibility_cutoff, exp_file, forester_state )
	
		# --- prepare subject --- #
		peptide_file = output_dir + "subject.fasta"
		subject_name_file = output_dir + "subject_names.txt"
		if seqtype == "pep":
			sys.stdout.write( "INFO: Input sequences are peptides >> candidate detection can start directly.\n" )
			sys.stdout.flush()
			if not os.path.isfile( peptide_file ):
				generate_subject_file( peptide_file, subject_name_file, subject )
		elif seqtype == "rna":
			sys.stdout.write( "INFO: Input sequences are DNA (transcripts expected) >> in silico translation will be performed in all 6 frames ...\n" )
			sys.stdout.flush()
			min_len_cutoff = 50	#minimal number of amino acids to keep a sequence as peptide
			translate_to_generate_pep_file( peptide_file, subject_name_file, subject, min_len_cutoff )
		
		
		# --- prepare query files --- #
		fasta_files = glob.glob( bait_seq_data_dir + "*.fa" ) + glob.glob( bait_seq_data_dir + "*.fasta" )
		blast_query_dir = output_dir + "blast_query/"
		if not os.path.exists( blast_query_dir ):
			os.makedirs( blast_query_dir )
		
		query_file = output_dir + "QUERY.fasta"
		name_mapping_table = output_dir + "NAME_MAPPING.txt"
		if not os.path.isfile( query_file ):
			with open( name_mapping_table, "w" ) as nmt:
				query_files = []
				for fasta in fasta_files:
					qf = generate_query( fasta, blast_query_dir, nmt )
					query_files.append( qf )
			p = subprocess.Popen( args="cat " + " ".join( query_files ) + ' > ' + query_file, shell=True  )
			p.communicate()
		
		# --- get subject peptide sequences if DNA is provided --- #
		if seqtype == "dna":
			dna_folder = output_dir + "DNA_screen/"
			dna_screener( subject, peptide_file, subject_name_file, dna_folder, query_file, mafft, cpus, max_gene_size, makeblastdb, tblastn )
		
		
		subject_name_mapping_table = load_subject_name_mapping_table( subject_name_file )
		
		# --- run BLAST searches --- #
		blast_result_file = output_dir + "blast_results.txt"
		blast_db = output_dir + "blastdb"
		if not os.path.isfile( blast_result_file ):
			p = subprocess.Popen( args= makeblastdb + " -in " + peptide_file + " -out " + blast_db + " -dbtype prot", shell=True )
			p.communicate()
			p = subprocess.Popen( args= blastp + " -query " + query_file + " -db " + blast_db + " -out " + blast_result_file + " -outfmt 6 -evalue 0.00001 -num_threads " + str( cpus ), shell=True )
			p.communicate()
		self_blast_result_file = output_dir + "self_blast_results.txt"
		self_blast_db = output_dir + "self_blastdb"
		if not os.path.isfile( self_blast_result_file ):
			p = subprocess.Popen( args= makeblastdb + " -in " + query_file + " -out " + self_blast_db + " -dbtype prot", shell=True )
			p.communicate()
			p = subprocess.Popen( args= blastp + " -query " + query_file + " -db " + self_blast_db + " -out " + self_blast_result_file + " -outfmt 6 -evalue 0.00001 -num_threads " + str( cpus ), shell=True )
			p.communicate()
		
		
		# --- load sequence data --- #
		peps = load_sequences( peptide_file )	#target peptide sequences
		ref_seqs = load_ref_seqs( query_file, name_mapping_table )
		
		
		# --- load BLAST results or classify based on phylogenetic tree --- #
		self_scores = load_self_BLAST_hit_scores( self_blast_result_file )
		
		if treemethod != "none":
			tree_tmp = output_dir + "tree_tmp/"
			if not os.path.exists( tree_tmp ):
				os.makedirs( tree_tmp )
			blast_hits = tree_based_classification( blast_result_file, self_scores, score_ratio_cutoff, similarity_cutoff,
																			tree_tmp, treemethod, fasttree, iqtree, mafft, peps, ref_seqs, possibility_cutoff
																			)
		else:
			blast_hits = load_BLAST_results( blast_result_file, self_scores, score_ratio_cutoff, similarity_cutoff, possibility_cutoff )
		
		# --- analyse candidates in global alignments --- #
		tmp_dir = output_dir + "tmp/"
		if not os.path.exists( tmp_dir ):
			os.makedirs( tmp_dir )
		
		alignment_per_candidate, sim_matrix_per_gene, candidates_by_gene = generate_global_alignments( mafft, peps, blast_hits, tmp_dir, ref_seqs )
		sim_matrix_folder = output_dir + "similarity_matrix/"
		if not os.path.exists( sim_matrix_folder ):
			os.makedirs( sim_matrix_folder )
		sim_per_pep = generate_sim_matrix_output_files( sim_matrix_folder, sim_matrix_per_gene, subject_name_mapping_table )
		
		
		# --- check conserved residues in alignment --- #
		if len( pos_data_dir ) > 1:
			pos_data_files = glob.glob( pos_data_dir + "*.txt" ) + glob.glob( pos_data_dir + "*.res" )
		else:
			pos_data_files = []
		pos_data_per_gene, regions_per_gene = load_pos_data_per_gene( pos_data_files )
		cons_res_matrix_folder = output_dir + "conserved_residues/"
		if not os.path.exists( cons_res_matrix_folder ):
			os.makedirs( cons_res_matrix_folder )
		cons_pos_per_pep, cons_pos_per_pep_extra = check_cons_res( cons_res_matrix_folder, pos_data_per_gene, alignment_per_candidate, candidates_by_gene, subject_name_mapping_table )
		
		
		# --- check conserved regions in alignment --- #
		cons_reg_matrix_folder = output_dir + "conserved_regions/"
		if not os.path.exists( cons_reg_matrix_folder ):
			os.makedirs( cons_reg_matrix_folder )
		cons_reg_per_pep = check_cons_reg( cons_reg_matrix_folder, regions_per_gene, alignment_per_candidate, candidates_by_gene, subject_name_mapping_table )
		
		# --- generate final peptide files and summary file --- #
		final_pep_folder = output_dir + "final_pep_files/"
		if not os.path.exists( final_pep_folder ):
			os.makedirs( final_pep_folder )
		
		summary_file = output_dir + "summary.txt"
		summary, final_file_per_gene = generate_final_pep_files( peps, final_pep_folder, candidates_by_gene, xsimcut, xconsrescut, xconsregcut, sim_per_pep, cons_pos_per_pep, cons_reg_per_pep, summary_file, subject_name_mapping_table )
		
		# --- generate summary HTML --- #
		html_file = output_dir + "SUMMARY.html"
		
		pathway = get_pathway( pathway_file, candidates_by_gene )
		generate_summary_html( html_file, summary, peps, cons_pos_per_pep_extra, pathway, pos_data_per_gene )

		#build phylogenetic tree with landmark sequences
		if forester_state:
			forester_output_folder = output_dir + "gene_trees/"
			if treemethod == "none":
				if len( fasttree ) == 0:
					fasttree = "FastTree"
			forester( final_pep_folder, forester_output_folder, bait_seq_data_dir, mafft, treemethod, fasttree, iqtree, occupancy=0.1 )
			
		if len( exp_file ) > 0:
			coexp_output_folder = output_dir + "coexpression/"
			if not os.path.exists( coexp_output_folder ):
				os.makedirs( coexp_output_folder )
			output_file = coexp_output_folder + "coexp_summary.txt"
			automatic_coexp( exp_file, output_file, summary_file, coexp_output_folder, rcutoff, pcutoff, minexpression )


def seqtype_check( sequences ):
	"""! @brief checking the type of provided input sequences """
	
	sequences = list( sequences.values() )
	one_seq_string = "".join( sequences[:min( [ len(sequences), 100 ] )] ).upper()
	acgtn = one_seq_string.count( "A" ) + one_seq_string.count( "C" ) + one_seq_string.count( "G" ) + one_seq_string.count( "T" ) + one_seq_string.count( "N" )
	if len( one_seq_string ) > 0:
		if float( acgtn ) / len( one_seq_string ) > 0.9:
			return "rna"
		else:
			return "pep"
	else:
		return "pep"


def main( arguments ):
	"""! @brief run everything """
	
	bait_seq_data_dir = arguments[ arguments.index('--baits')+1 ]
	
	if '--subjectdir' in arguments:
		subject_dir = arguments[ arguments.index('--subjectdir')+1 ]
		subjects = []
		for ext in [ "*.faa", "*.fa", "*.fasta", "*.faa.gz", "*.fa.gz", "*.fasta.gz", "*.FAA", "*.FA", "*.FASTA", "*.FAA.GZ", "*.FA.GZ", "*.FASTA.GZ" ]:
			subjects += glob.glob( subject_dir + ext )
		output_dir = arguments[ arguments.index('--out')+1 ]
		if output_dir[-1] != '/':
			output_dir += "/"
		output_dirs = []
		for subject in subjects:
			ID = subject.split('/')[-1].split('.')[0]
			output_dirs.append( output_dir + ID + '/' )
	else:
		subjects = [ arguments[ arguments.index('--subject')+1 ] ]
		output_dirs = [ arguments[ arguments.index('--out')+1 ] ]
	
	if '--positions' in arguments:
		pos_data_dir = arguments[ arguments.index('--positions')+1 ]
	else:
		if '--residues' in arguments:
			pos_data_dir = arguments[ arguments.index('--residues')+1 ]
		else:
			pos_data_dir = ""
	
	if '--seqtype' in arguments:
		seqtype = arguments[ arguments.index('--seqtype')+1 ]
		#pep = peptide collection
		#rna = transcriptome sequence
		#dna = genome sequence
	else:
		seqtype = "pep"
	
	if '--mafft' in arguments:
		mafft = arguments[ arguments.index('--mafft')+1 ]
	else:
		mafft = "mafft"
	
	if '--blastp' in arguments:
		blastp = arguments[ arguments.index('--blastp')+1 ]
	else:
		blastp = "blastp"
	
	if '--tblastn' in arguments:
		tblastn = arguments[ arguments.index('--tblastn')+1 ]
	else:
		tblastn = "tblastn"
	
	if '--makeblastdb' in arguments:
		makeblastdb = arguments[ arguments.index('--makeblastdb')+1 ]
	else:
		makeblastdb = "makeblastdb"
	
	if '--fasttree' in arguments:
		fasttree = arguments[ arguments.index('--fasttree')+1 ]
		treemethod = "fasttree"
		sys.stdout.write(  "INFO: classification of candidates will be based on phylogenetic trees.\n" )
		sys.stdout.flush()
	else:
		fasttree = ""
	
	if '--iqtree' in arguments:
		iqtree = arguments[ arguments.index('--iqtree')+1 ]
		treemethod = "iqtree"
		sys.stdout.write(  "INFO: classification of candidates will be based on phylogenetic trees.\n" )
		sys.stdout.flush()
	else:
		iqtree = ""
	
	if fasttree == "" and iqtree == "":
		treemethod = "none"
		sys.stdout.write(  "INFO: classification of candidates will be based on BLAST hit similarity.\n" )
		sys.stdout.flush()
	
	if '--cpus' in arguments:
		cpus = int( arguments[ arguments.index('--cpus')+1 ] )
	else:
		cpus = 10
	
	# --- parameters for filtering of BLAST results --- #
	if '--scoreratio' in arguments:
		score_ratio_cutoff = float( arguments[ arguments.index('--scoreratio')+1 ] )
	else:	
		score_ratio_cutoff = 0.3
	
	if '--simcut' in arguments:
		similarity_cutoff = float( arguments[ arguments.index('--simcut')+1 ] )
	else:
		similarity_cutoff = 40.0	#value in percent
	
	if '--genesize' in arguments:
		max_gene_size = int( arguments[ arguments.index('--genesize')+1 ] )
	else:
		max_gene_size = 5000
	
	if '--minsim' in arguments:
		xsimcut = float( arguments[ arguments.index('--minsim')+1 ] )
	else:
		xsimcut = 0.4	#minimal similarity in global alignment to keep candidate
	
	if '--minres' in arguments:
		xconsrescut = float( arguments[ arguments.index('--minres')+1 ] )
	else:
		xconsrescut = -1	#minimal ratio of conserved residues to keep candidate
	
	if '--minreg' in arguments:
		xconsregcut = float( arguments[ arguments.index('--minreg')+1 ] )
	else:
		xconsregcut = -1	#minimal similarity of conserved regions to keep candidate (deactivated by default)
	
	if '--possibilities' in arguments:
		possibility_cutoff = int( arguments[ arguments.index('--possibilities')+1 ] )
	else:
		possibility_cutoff = 3
	
	if '--checks' in arguments:
		checks = arguments[ arguments.index('--checks')+1 ]
	else:
		checks = "on"
	
	if '--pathway' in arguments:
		pathway_file = arguments[ arguments.index('--pathway')+1 ]
	else:
		pathway_file = ""
	
	for subject in subjects:
		if not os.path.isfile( subject ):
			sys.stdout.write( "ERROR: subject file not detected: " + subject + "\n" )
			sys.stdout.flush()
	
	if not os.path.exists( bait_seq_data_dir ):
		sys.exit( "ERROR: bait sequence folder not detected!" )
	
	time.sleep( 10 )
	
	if '--exp' in arguments:
		exp_file = arguments[ arguments.index('--exp')+1 ]
		sys.stdout.write( "Co-expression analysis activated.\n" )
		sys.stdout.flush()
		if not os.path.isfile( exp_file ):
			sys.stdout.write( "ERROR: gene expression file not detected - " + exp_file + "\n" )
			sys.stdout.flush()
	else:
		exp_file = ""
	
	if '--rcut' in arguments:
		rcutoff = float( arguments[ arguments.index('--rcut')+1 ] )
	else:
		rcutoff = 0.3
	
	if '--pcut' in arguments:
		pcutoff = float( arguments[ arguments.index('--pcut')+1 ] )
	else:
		pcutoff = 0.05
	
	if '--minexp' in arguments:
		minexpression = int( arguments[ arguments.index('--minexp')+1 ] )
	else:
		minexpression = 30
	
	if '--forester' in arguments:
		forester_state = True
		sys.stdout.write( "Gene tree construction activated.\n" )
		sys.stdout.flush()
	else:
		forester_state = False
	
	# --- run KIPEs for each supplied subject file --- #
	sys.stdout.write( "number of subjects to process: " + str( len( subjects ) ) + "\n" )
	sys.stdout.flush()
	for xxx, subject in enumerate( subjects ):
		# --- check seqtype --- #
		if seqtype == "pep":
			if os.path.isfile( subject ):
				sequences =  load_sequences( subject )
				real_seqtype = seqtype_check( sequences )
				if real_seqtype != "pep":
					sys.stdout.write( "WARNING: PROVIDED SEQUENCES APPEAR TO BE NUCLEOTIDES NOT POLYPPEPTIDES! Please reconsider your --seqtype settings. \n" )
					sys.stdout.flush()
		
		# --- start KIPEs --- #
		KIPEs( 	bait_seq_data_dir, output_dirs[ xxx ], subject, pos_data_dir, seqtype,
					mafft, blastp, tblastn, makeblastdb, treemethod, fasttree, iqtree, pathway_file,
					cpus, score_ratio_cutoff, similarity_cutoff, max_gene_size,
					xsimcut, xconsrescut, xconsregcut, checks, possibility_cutoff,
					exp_file, forester_state, rcutoff, pcutoff, minexpression
				)


if '--baits' in sys.argv and '--out' in sys.argv and '--subject' in sys.argv:
	main( sys.argv )
elif '--baits' in sys.argv and '--out' in sys.argv and '--subjectdir' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
