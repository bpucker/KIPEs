### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.18 ###

__usage__ = """
					python cgi.py
					--baits <FOLDER_WITH_BAIT_SEQ_FILES>
					--positions <FOLDER_WITH_POSITION_FILES>
					--out <OUTPUT_DIR>
					--subject <SUBJECT_FILE (peptide,transcript,genomic sequences)>
					
					optional:
					--seqtype <TYPE_OF_SUBJECT_SEQUENCES(pep|rna|dna)>[pep]
					--mafft <PATH_TO_MAFFT>[mafft]
					--cpus <INT, NUMBER_OF_BLAST_THREADS>[10]
					--scoreratio <FLOAT, BLAST_SCORE_RATIO_CUTOFF>[0.3]
					--simcut <FLOAT, MINIMAL_BLAST_HIT_SIMILARITY_IN_PERCENT>[40.0]
					--genesize <INT, SIZE_OF_GENE_FOR_GROUPING_OF_EXON_HITS>[5000]
					--minsim <FLOAT, MINIMAL_SIMILARITY_IN_GLOBAL_ALIGNMENT>[0.4]
					--minres <FLOAT, MINIMAL_PROPORTION_OF_CONSERVED_RESIDUES>[0.8]
					--minreg <FLOAT, MINIMAL_PROPORTION_OF_CONSERVED_REGIONS>[off]
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import os, glob, sys, time, re
from operator import itemgetter

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
					out.write( '>' + ID + "_%_" + str( counter ) + '\n' + "".join( seq ) + '\n' )
					nmt.write( ID + "_%_" + str( counter ) + '\n' + line.strip()[1:].replace( "\t", "   " ) + "\t" )
					counter += 1
					seq = []
				else:
					seq.append( line.strip() )
				line = f.readline()
			out.write( '>' + ID + "_%_" + str( counter ) + '\n' + "".join( seq ) + '\n' )
			nmt.write( ID + "_%_" + str( counter ) + '\n' )
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


def load_BLAST_results( blast_result_file, self_scores, score_ratio_cutoff, similarity_cutoff ):
	"""! @brief load BLAST results """
	
	valid_blast_hits = {}
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if float( parts[2] ) > similarity_cutoff:	#similarity is sufficient
				if float( parts[-1] ) > score_ratio_cutoff * self_scores[ parts[0] ]:	#substantial part of query is matched
					try:
						if float( parts[-1] ) > valid_blast_hits[ parts[1] ]['score']:
							valid_blast_hits[ parts[1] ] = { 'gene': parts[0].split('_%_')[0], 'score': float( parts[-1] ) }
					except KeyError:
						valid_blast_hits.update( { parts[1]: { 'gene': parts[0].split('_%_')[0], 'score': float( parts[-1] ) } } )
			line = f.readline()
	return valid_blast_hits


def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip()
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:]
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
	for key in seqs.keys():
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
	for gene in candidates_by_gene.keys():
		gene_sim_matrix = {}
		for candidate in candidates_by_gene[ gene ]:
			candiate_sim_matrix = {}
			for query in query_names_by_gene[ gene ]:
				try:
					sim = calculate_similarity( alignment_per_candidate[ candidate ][ candidate ], alignment_per_candidate[ candidate ][ query ] )
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
	for candidate in blast_hits.keys():
		gene = blast_hits[ candidate ]['gene']
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
		os.popen( mafft + " " + tmp_seq_file + " > " + aln_file )
		alignment_per_candidate.update( { candidate: load_alignment( aln_file, tmp_mapping ) } )
	
	# --- get all query sequence names per gene --- #
	query_names_by_gene = {}
	for ref in [ x for sublist in ref_seqs.values() for x in sublist]:	#walk through all query sequences
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
	for gene in sim_matrix_per_gene.keys():
		output_file = sim_matrix_folder + gene + "_sim_matrix.txt"
		with open( output_file, "w" ) as out:
			data = sim_matrix_per_gene[ gene ]	#similarity matrices for all candidates
			queries = sorted( data.values()[0].keys() )
			out.write( "\t".join( [ "candidate" ] + queries ) + '\n' )
			for candidate in sorted( data.keys() ):
				new_line = [ subject_name_mapping_table[ candidate ] ]
				for query in queries:
					new_line.append( 100.0*data[ candidate ][ query ] )
				sim_per_pep.update( { candidate: sum( new_line[1:] ) / len( new_line[1:] ) } )
				out.write( "\t".join( map( str, new_line ) ) + '\n' )
	return sim_per_pep


def load_pos_data_per_gene( pos_data_files ):
	"""! @brief load conserved residue positions from given data file """
	
	pos_per_gene = {}
	regions_per_gene = {}
	for filename in pos_data_files:
		gene = filename.split('/')[-1].split('.')[0]
		with open( filename, "r" ) as f:
			#SEQ_NAME		X1	X2	X3	... X10, ... X100
			#SEQ_NAME		REGION1,START,END		REGION2,START,END		...
			line = f.readline()
			while line:
				data_row = line.strip().split('\t')
				if "," in "".join( data_row[1:] ):
					regions = []
					for reg in data_row[1:]:
						parts = reg.split(',')
						regions.append( { 'name': parts[0], 'start': int( parts[1] ), 'end': int( parts[2] ) } )
					regions_per_gene.update( { gene: { 'seq': data_row[0], 'regions': regions } } )
				else:
					residues = []
					for res in data_row[1:]:
						residues.append( { 'aa': res[0], 'pos': int( res[1:] ) } )
					pos_per_gene.update( { gene: { 'seq': data_row[0], 'residues': residues } } )
				line = f.readline()
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
	for res in residues:
		alignment_pos = get_alignment_pos( ref_aln, res['pos']-1 )
		results.append( can_aln[ alignment_pos ] == res['aa'] )
	return results


def check_cons_res( cons_res_matrix_folder, pos_data_per_gene, alignment_per_candidate, candidates_by_gene, subject_name_mapping_table ):
	"""! @brief check all candidate sequences for conserved residues and generate result tables """
	
	cons_pos_per_pep = {}
	for gene in candidates_by_gene.keys():
		candidates = candidates_by_gene[ gene ]
		try:
			info = pos_data_per_gene[ gene ]
			output_file = cons_res_matrix_folder + gene + "_conserved_residues.txt"
			residues = sorted( info['residues'], key=itemgetter('pos') )
			with open( output_file, "w" ) as out:
				header = [ "candidate" ]
				for each in residues:
					header.append( each['aa'] + str( each['pos'] ) )
				out.write( "\t".join( header ) + '\n' )
				for candidate in candidates:
					print "gene: " + gene
					can_aln = alignment_per_candidate[ candidate ][ candidate ]
					ref_aln = alignment_per_candidate[ candidate ][ info['seq'] ]
					results = check_alignment_for_cons_res( can_aln, ref_aln, residues )
					out.write( "\t".join( map( str, [ subject_name_mapping_table[ candidate ] ] + results  ) ) + '\n' )
					cons_pos_per_pep.update( { candidate: sum( results ) / len( results ) } )
		except KeyError:
			print "no information (conserved residues) available about gene: " + gene
	return cons_pos_per_pep


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
	for gene in candidates_by_gene.keys():
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
					can_aln = alignment_per_candidate[ candidate ][ candidate ]
					ref_aln = alignment_per_candidate[ candidate ][ info['seq'] ]
					results = check_alignment_for_cons_reg( can_aln, ref_aln, regions )
					out.write( "\t".join( map( str, [ subject_name_mapping_table[ candidate ] ] + results  ) ) + '\n' )
					cons_reg_per_pep.update( { candidate: sum( results ) / len( results ) } )
		except KeyError:
			print "no information (conserved regions) available about gene: " + gene
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
	for seq in sequences.values():
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
		for key in genes.keys():
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
			genes.update( { len( genes.keys() )+1: { 	'chr': hit['chr'],
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
	for query in data.keys():
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
			os.popen( mafft + " " + out_fasta_file + " > " + aln_file )
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
						print "WARNING: contig end causes truncated peptide sequence (N-terminal end)"
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
							print "WARNING: contig end causes truncated peptide sequence (C-terminal end)"
							seq_parts.append( extend_to_start( part['sstart'], part['qstart'], part['send'], seq ) )
					#insufficient upstream sequence
					else:
						#sufficient downstream DNA sequence
						if len( seq )-part['send'] > 3*( qlen-part['qend'] ):
							seq_parts.append( extend_to_stop( part['sstart'], part['send'], qlen, part['qend'], seq ) )
							#seq_parts.append( seq[ part['sstart']-1:part['send']+3*( qlen-part['qend'] ) ] )
							print "WARNING: contig end causes truncated peptide sequence (N-terminal end)"
						#insufficient downstream DNA sequence
						else:
							print "WARNING: contig end causes truncated peptide sequence (N- and C-terminal end)"
							seq_parts.append( seq[ part['sstart']-1:part['send'] ] )
		#last coding exon (more than 1 coding exons)
		elif idx == len( parts )-1:
			#direct hit next to previous one
			if parts[idx-1]['qend']+1 == part['qstart']:
				seq_parts.append( "NNN" )
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
						print "WARNING: contig end causes truncated peptide sequence (C-terminal end)"
						seq_parts.append( seq[ part['sstart']-1:part['send'] ] )
			#gap to previous coding exon hit
			elif parts[idx-1]['qend']+1 < part['qstart']:
				seq_parts.append( "NNNNNN" )
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
						print "WARNING: contig end causes truncated peptide sequence (C-terminal end)"
						seq_parts.append( missing_seq + seq[ part['sstart']-1:part['send'] ] )
			#overlap with previous coding exon
			else:
				seq_parts.append( "NNNNNNNNN" )
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
						print "WARNING: contig end causes truncated peptide sequence (C-terminal end)"
						seq_parts.append( seq[ part['sstart']-1+3*(1+parts[idx-1]['qend']-part['qstart']):part['send'] ] )
		#all enclosed coding exons
		else:
			#direct hit next to previous one
			if parts[idx-1]['qend']+1 == part['qstart']:
				seq_parts.append( "NNN" )
				### ADD SPLICE SITE CHECK
				seq_parts.append( seq[ part['sstart']-1:part['send'] ] )
			
			#gap
			elif parts[idx-1]['qend']+1 < part['qstart']:
				seq_parts.append( "NNNNNN" )
				# --- search for missing sequence --- #
				missing_seq = find_missing_seq_between_hits( 	qpep[ parts[idx-1]['qend']+1:part['qstart']+1 ],	#missing amino acids
																									seq[ parts[idx-1]['send']: part['sstart'] ],	#corresponding DNA sequence																									
																									genetic_code, tmp_folder, mafft
																								)		
				seq_parts.append( missing_seq + seq[ part['sstart']-1:part['send'] ] )
			
			#overlap
			else:
				seq_parts.append( "NNNNNNNNN" )
				seq_parts.append( seq[ part['sstart']-1+3*(1+parts[idx-1]['qend']-part['qstart']):part['send'] ] )
			
	#find start M
	#find stop codon
	lengths = []
	for each in seq_parts:
		lengths.append( len( each ) / 3.0 )
	#print lengths
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
	for query in gene_groups_per_query.keys():
		genes = gene_groups_per_query[ query ]
		pep_seqs = {}
		for gene in genes.keys():
			chromosome = genes[ gene ]['chr']
			parts = genes[ gene ]['parts']
			pep_seq = construct_pep_seq( parts, genome_seq[ chromosome ], genetic_code, len( querys[ query ] ), querys[ query ], tmp_folder, mafft )
			pep_seqs.update( { gene: translate( "".join( pep_seq ), genetic_code ) } )	#TRANSLATION
		all_pep_seqs.update( { query: pep_seqs } )
	return all_pep_seqs


def dna_screener( subject, peptide_file, subject_name_file, dna_folder, query_file, mafft, cpus, max_gene_size ):
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
		os.popen( "makeblastdb -in " + dna_file + " -out " + blast_db + " -dbtype nucl" )
		os.popen( "tblastn -query " + query_file + " -db " + blast_db + " -out " + blast_result_file + " -outfmt 6 -evalue 0.01 -word_size 4 -num_threads " + str( cpus ) )
	
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
			for query in pep_seqs.keys():
				for gene in pep_seqs[ query ].keys():
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
	with open( summary_file, "w") as summary:
		summary.write( "ID\tGene\tSimilarity\tConservedResidues\tConservedRegions\n" )
		for gene in candidates_by_gene.keys():
			output_file = final_pep_folder + gene + ".fasta"
			with open( output_file, "w" ) as out:
				candidates = candidates_by_gene[ gene ]
				values_for_sorting = []
				for candidate in candidates:
					try:
						psim = sim_per_pep[ candidate ]
					except KeyError:
						psim = 0
					try:
						pres = cons_pos_per_pep[ candidate ]
					except KeyError:
						pres = 0
					try:
						preg = cons_reg_per_pep[ candidate ]
					except KeyError:
						preg = 0
					
					values_for_sorting.append( { 	'id': candidate,
																			'sim': psim,
																			'res': pres,
																			'reg': preg,
																			'gene': gene
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
	return complete_summary


def main( arguments ):
	"""! @brief run everything """
	
	bait_seq_data_dir = arguments[ arguments.index('--baits')+1 ]
	pos_data_dir = arguments[ arguments.index('--positions')+1 ]
	output_dir = arguments[ arguments.index('--out')+1 ]
	subject = arguments[ arguments.index('--subject')+1 ]
	
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
		xconsrescut = 0.8	#minimal ratio of conserved residues to keep candidate
	
	if '--minreg' in arguments:
		xconsregcut = float( arguments[ arguments.index('--minreg')+1 ] )
	else:
		xconsregcut = -1	#minimal similarity of conserved regions to keep candidate (deactivated by default)
	
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )
	
	# --- prepare subject --- #
	peptide_file = output_dir + "subject.fasta"
	subject_name_file = output_dir + "subject_names.txt"
	if seqtype == "pep":
		print "Input sequences are peptides >> candidate detection can start directly."
		if not os.path.isfile( peptide_file ):
			generate_subject_file( peptide_file, subject_name_file, subject )
	elif seqtype == "rna":
		print "Input sequences are DNA (transcripts expected) >> in silico translation will be performed in all 6 frames ..."
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
		os.popen( "cat " + " ".join( query_files ) + ' > ' + query_file  )
	
	# --- get subject peptide sequences if DNA is provided --- #
	if seqtype == "dna":
		dna_folder = output_dir + "DNA_screen/"
		dna_screener( subject, peptide_file, subject_name_file, dna_folder, query_file, mafft, cpus, max_gene_size )
	
	
	subject_name_mapping_table = load_subject_name_mapping_table( subject_name_file )
	
	# --- run BLAST searches --- #
	blast_result_file = output_dir + "blast_results.txt"
	blast_db = output_dir + "blastdb"
	if not os.path.isfile( blast_result_file ):
		os.popen( "makeblastdb -in " + peptide_file + " -out " + blast_db + " -dbtype prot" )
		os.popen( "blastp -query " + query_file + " -db " + blast_db + " -out " + blast_result_file + " -outfmt 6 -evalue 0.00001 -num_threads " + str( cpus ) )
	self_blast_result_file = output_dir + "self_blast_results.txt"
	self_blast_db = output_dir + "self_blastdb"
	if not os.path.isfile( self_blast_result_file ):
		os.popen( "makeblastdb -in " + query_file + " -out " + self_blast_db + " -dbtype prot" )
		os.popen( "blastp -query " + query_file + " -db " + self_blast_db + " -out " + self_blast_result_file + " -outfmt 6 -evalue 0.00001 -num_threads " + str( cpus ) )
	
	
	# --- load BLAST results --- #
	self_scores = load_self_BLAST_hit_scores( self_blast_result_file )
	blast_hits = load_BLAST_results( blast_result_file, self_scores, score_ratio_cutoff, similarity_cutoff )
	
	
	# --- analyse candidates in global alignments --- #
	tmp_dir = output_dir + "tmp/"
	if not os.path.exists( tmp_dir ):
		os.makedirs( tmp_dir )
	peps = load_sequences( peptide_file )
	ref_seqs = load_ref_seqs( query_file, name_mapping_table )
	
	
	alignment_per_candidate, sim_matrix_per_gene, candidates_by_gene = generate_global_alignments( mafft, peps, blast_hits, tmp_dir, ref_seqs )
	sim_matrix_folder = output_dir + "similarity_matrix/"
	if not os.path.exists( sim_matrix_folder ):
		os.makedirs( sim_matrix_folder )
	sim_per_pep = generate_sim_matrix_output_files( sim_matrix_folder, sim_matrix_per_gene, subject_name_mapping_table )
	
	
	# --- check conserved residues in alignment --- #
	pos_data_files = glob.glob( pos_data_dir + "*.txt" )
	pos_data_per_gene, regions_per_gene = load_pos_data_per_gene( pos_data_files )
	cons_res_matrix_folder = output_dir + "conserved_residues/"
	if not os.path.exists( cons_res_matrix_folder ):
		os.makedirs( cons_res_matrix_folder )
	cons_pos_per_pep = check_cons_res( cons_res_matrix_folder, pos_data_per_gene, alignment_per_candidate, candidates_by_gene, subject_name_mapping_table )
	
	
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
	summary = generate_final_pep_files( peps, final_pep_folder, candidates_by_gene, xsimcut, xconsrescut, xconsregcut, sim_per_pep, cons_pos_per_pep, cons_reg_per_pep, summary_file, subject_name_mapping_table )

#build phylogenetic tree with landmark sequences


if '--baits' in sys.argv and '--positions' in sys.argv and '--out' in sys.argv and '--subject' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
