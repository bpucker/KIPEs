### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python get_cons_pos.py
					--in <INPUT_FASTA_FILE>
					--out <OUTPUT_FOLDER>
					--ref <REF_SEQ_NAME>
					--name <PREFIX_OF_RESULT_FILES>
					
					optional:
					--mafft <PATH_TO_MAFFT>[mafft]
					--mincons <FLOAT, CONSERVATION_OF_RESIDUE_TO_BE_CONSIDERED>[1.0]
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import os, sys

# --- end of imports --- #

def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip()
		if " " in header:
			header = header.split(' ')[0]
		if "\t" in header:
			header = header.split('\t')[0]
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:]
					if " " in header:
						header = header.split(' ')[0]
					if "\t" in header:
						header = header.split('\t')[0]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def get_conserved_pos( alignment, mincons, ref_name ):
	"""! @brief identify conserved positions in reference sequence """
	
	cons_positions = []
	aln_idx = -1
	for idx, aa in enumerate( alignment[ ref_name ] ):
		if aa != "-":
			aln_idx += 1
		counter = 0
		for seq in alignment.values():
			if seq[ idx ] == aa:
				counter += 1
		if counter / len( alignment.values() ) >= mincons:
			cons_positions.append( aa + str( aln_idx+1 ) )
	return cons_positions


def main( arguments ):
	"""! @brief run everything """
	
	input_file = arguments[ arguments.index('--in')+1 ]
	output_folder = arguments[ arguments.index('--out')+1 ]
	ref_name = arguments[ arguments.index('--ref')+1 ]
	
	if '--mafft' in arguments:
		mafft = arguments[ arguments.index('--out')+1 ]
	else:
		mafft = "mafft"
	if '--mincons' in arguments:
		mincons = float( arguments[ arguments.index('--mincons')+1 ] )
	else:
		mincons = 1.0
	
	name = arguments[ arguments.index('--name')+1 ]
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	
	# --- check for presence of reference sequence in file --- #
	peps = load_sequences( input_file )
	try:
		peps[ ref_name ]
	except KeyError:
		sys.exit( "ERROR: Reference sequence name " + ref_name + " not detected." )
	
	# --- generate alignment --- #
	aln_file = output_folder + name + ".fasta.aln"
	os.popen( mafft + " " + input_file + " > " + aln_file )
	
	# --- load alignment and identify conserved residues --- #
	alignment = load_sequences( aln_file )
	cons_positions = get_conserved_pos( alignment, mincons, ref_name )
	print "number of identified conserved amino acid residues: " + str( len( cons_positions ) )
	
	output_file = output_folder + name + ".txt"
	with open( output_file, "w" ) as out:
		out.write( "\t".join( [ ref_name ] + cons_positions ) + '\n' )


if '--in' in sys.argv and '--out' in sys.argv and '--ref' in sys.argv and '--name' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
