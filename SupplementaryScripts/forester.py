### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

### based on KIPEs: https://doi.org/10.3390/plants9091103 ###

__usage__ = """
					python enrichment_check.py
					--in <INPUT_FOLDER>
					--out <OUTPUT_FOLDER>
					--ref <REFERENCE_SEQUENCE_FOLDER>
					
					optional:
					--mafft <PATH_TO_MAFFT>[mafft]
					--fasttree <PATH_TO_FASTTREE>[FastTree]
					--occ <MINIMAL_ALIGNMENT_OCCUPANCY>[0.1]
					"""


import re, os, sys, glob
from operator import itemgetter

# --- end of imports --- #


def load_sequences( multiple_fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( multiple_fasta_file ) as f:
		header = f.readline()[1:].strip().replace( " ", "" ).replace('(', '_').replace(')', '_').replace( "	", "" ).replace('.', '_').replace(':', '_').replace(';', '_')	#.replace( "-", "_" ).replace('[', '_').replace(']', '_')
		#if header[0] in "0123456789":
		#	header = "x" + header
		seq = ""
		line = f.readline()
		while line:
			if line[0] == '>':
				if len( seq ) > 0:
					sequences.update( { header: seq } )
				header = line.strip()[1:].replace( " ", "" ).replace('(', '_').replace(')', '_').replace( "	", "" ).replace('.', '_').replace(':', '_').replace(';', '_')	#.replace( "-", "_" ).replace('[', '_').replace(']', '_')
				#if header[0] in "0123456789":
				#	header = "x" + header
				seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		if len( seq ) > 0:
			sequences.update( { header: seq } )
	return sequences


def modify_FASTA_names( input_file, output_file, mapping_table_file ):
	"""! @brief modify names in FASTA file """
	
	seqs = load_sequences( input_file )
	mapping_table = {}
	with open( output_file, "w" ) as out:
		for idx, key in enumerate( seqs.keys() ):
			out.write( '>seq' + str( idx ).zfill(5) + '\n' + seqs[ key ] + '\n' )
			mapping_table.update( { key: 'seq' + str( idx ).zfill(5) } )
	
	with open( mapping_table_file, "w" ) as out:
		for key in mapping_table.keys():
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
	
	for key in mapping_table.keys():
		tree = tree.replace( key, mapping_table[ key ] )
	
	with open( output_tree_file, "w" ) as out:
		out.write( tree )


def construct_tree( input_file,  output_dir, mafft, fasttree, occupancy, name ):
	"""! @brief handle tree construction """
		
	if output_dir[-1] != '/':
		output_dir += "/"
	
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )
	
	mod_FASTA = output_dir + "names_modified.fasta"
	mapping_table_file = output_dir + "seq_names_mapping_table.txt"
	modify_FASTA_names( input_file, mod_FASTA, mapping_table_file )
	
	alignment_file = mod_FASTA + ".aln"
	os.popen( " ".join( [ mafft, mod_FASTA, ">", alignment_file, "2>", alignment_file+".log" ] ) )
	
	clean_alignment_file = alignment_file + ".cln"
	alignment_trimming( alignment_file, clean_alignment_file, occupancy )
	
	tree_file = clean_alignment_file + ".tre"
	os.popen( " ".join( [ fasttree, "-wag -nosupport <", clean_alignment_file, ">", tree_file, "2>", tree_file+".log" ] ) )
	
	output_tree_file = output_dir + name + "FINAL_TREE.tre"
	modify_names_in_tree( tree_file, output_tree_file, mapping_table_file )
	
	return output_tree_file


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


def alignment_trimming( aln_file, cln_aln_file, occupancy ):
	"""! @brief remove all alignment columns with insufficient occupancy """
	
	alignment = load_alignment( aln_file, {} )
	# --- if there is an alignment (expected case) 
	if len( alignment.keys() ) > 0:
		# --- identify valid residues in aligned sequences (columns with sufficient occupancy) --- #
		valid_index = []
		for idx, aa in enumerate( alignment.values()[0] ):
			counter = 0
			for key in alignment.keys():
				if alignment[ key ][ idx ] != "-":
					counter += 1
			if counter / float( len( alignment.keys() ) ) > occupancy:
				valid_index.append( idx )
		
		# --- generate new sequences --- #
		with open( cln_aln_file, "w" ) as out:
			for key in alignment.keys():
				seq = alignment[ key ]
				new_seq = []
				for idx in valid_index:
					new_seq.append( seq[ idx ] )
				out.write( ">" + key + '\n' + "".join( new_seq ) + '\n' )
	# --- just in case the alignment file is empyt (is this possible?) ---#
	else:
		with open( cln_aln_file, "w" ) as out:
			out.write( "" )


def main( arguments ):
	"""! @brief run everything """
	
	input_folder = arguments[ arguments.index('--in')+1 ]	#KIPEs output "final_pep_files"
	if not "final_pep_files" in input_folder:
		sys.stdout.write( "WARNING: the input folder does not contain 'final_pep_files', but the path to this KIPEs output folder should be provided.\n" )
		sys.stdout.flush()
	output_folder = arguments[ arguments.index('--out')+1 ]	#output folder
	refseq_folder = arguments[ arguments.index('--ref')+1 ]	#folder with reference sequences
	
	if input_folder[-1] != "/":
		input_folder += "/"
	
	if output_folder[-1] != "/":
		output_folder += "/"
	
	if '--mafft' in arguments:
		mafft = arguments[ arguments.index( '--mafft' )+1 ]
	else:
		mafft = "mafft"
	
	if '--fasttree' in arguments:
		fasttree = arguments[ arguments.index( '--fasttree' )+1 ]
	else:
		fasttree = "FastTree"
	
	if '--occ' in arguments:
		occupancy = float( arguments[ arguments.index('--occ') + 1 ] )
	else:
		occupancy = 0.1
	
	if '--clean' in arguments:
		clean = True
	else:
		clean = False
	
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
			os.popen( "cat " + filename + " " + ref_seq_file + " > " + tree_input_file )
			tree_tmp_folder = output_folder + ID + "_tmp/"

			tree = construct_tree( tree_input_file, tree_tmp_folder, mafft, fasttree, occupancy, ID )
			
			os.popen( "cp " + tree + " " + output_folder + ID + ".tre" )
			if clean:
				os.popen( "rm -r " + tree_tmp_folder )	#remove temp folders
	

if '--in' in sys.argv and '--out' in sys.argv and "--ref" in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
