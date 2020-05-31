### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python get_KIPEs_seqs.py
					--in <INPUT_FOLDER>
					--out <OUTPUT_FILE>
					--gene <GENE_NAME>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import glob, os, sys

# --- end of imports --- #


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


def main( arguments ):
	"""! @brief run everything """
	
	gene = arguments[ arguments.index('--gene')+1 ]
	input_folder = arguments[ arguments.index('--in')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]

	if input_folder[-1] != "/":
		input_folder += "/"

	filenames = glob.glob( input_folder + "*/final_pep_files/" + gene + ".fasta" )

	with open( output_file, "w" ) as out:
		for filename in filenames:
			ID = filename.split('/')[-3]
			seqs = load_sequences( filename )
			for idx, key in enumerate( seqs.keys() ):
				out.write( '>' + ID + "@" + str( idx+1 ) + "\n" + seqs[ key ] + '\n' )


if '--in' in sys.argv and '--out' in sys.argv and '--gene' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
