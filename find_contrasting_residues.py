### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python find_contrasting_residues.py
					--in <INPUT_FASTA_FILE>
					--out <OUTPUT_FOLDER>
					--ref <REFERENCE_SEQUENCE_NAME>
					--pos <COMMA_SEPARATED_NAMES_OF_POSITIVE_SEQUENCES>
					--neg <COMMA_SEPARATED_NAMES_OF_NEGATIVE_SEQUENCES>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""


import os, sys

# --- end of imports --- #


def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip()
		if "\t" in header:
			header = header.split('\t')[0]
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:]
					if "\t" in header:
						header = header.split('\t')[0]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def main( arguments ):
	"""! @brief run everything """

	input_file = arguments[ arguments.index('--in')+1 ]
	output_folder = arguments[ arguments.index('--out')+1 ]
	pos_IDs = arguments[ arguments.index('--pos')+1 ].split(',')
	neg_IDs = arguments[ arguments.index('--neg')+1 ].split(',')
	ref_ID = arguments[ arguments.index('--ref')+1 ]
	
	
	if '--mafft' in arguments:
		mafft = arguments[ arguments.index('--mafft')+1 ]
	else:
		mafft = "mafft"
	
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	if output_folder[-1] != "/":
		output_folder += "/"
	
	aln_file = output_folder + "alignment.fasta.aln"
	os.popen( mafft + " " + input_file + " > " + aln_file )

	alignment = load_sequences( aln_file )
	
	try:
		alignment[ ref_ID ]
	except KeyError:
		sys.exit( "ERROR: reference ID is not present in alignment!" )

	results = []
	for idx, aa in enumerate( alignment[ ref_ID ] ):
		if aa != '-':
			pos_counter = 0
			for ID in pos_IDs:	#check for conservation in all sequences of this type
				if alignment[ ID ][ idx ] == aa:
					pos_counter += 1
			if pos_counter == len( pos_IDs ):	#perfect conservation
				neg_counter = 0
				for ID in neg_IDs:
					if alignment[ ID ][ idx ] == aa:
						neg_counter += 1
				if neg_counter == 0:	#valid position (only conserved in seq of interest)
					position = (idx+1) - alignment[ ref_ID ][ :idx].count( '-' )
					results.append( aa + str( position ) )
	
	print "number of detected informative residues: " + str( len( results ) )
	
	output_file = output_folder + "results.txt"
	with open( output_file, "w" ) as out:
		out.write( ref_ID + "\t" + "\t".join( results ) + '\n' )


if '--in' in sys.argv and '--ref' in sys.argv and '--pos' in sys.argv and '--neg' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
