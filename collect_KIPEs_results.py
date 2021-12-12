### Boas Pucker ###
### b.pucker@tu-braunschweig.de ###

__version__ = "0.1"
__usage__ = """
					python3 collect_KIPEs_results.py
					--in <INPUT_FOLDER>
					--out <OUTPUT_FILE>
					--genes <GENES_OF_INTEREST>
					
					optional:
					--mapping <NAME_MAPPING_TABLE>
					--minsim <MINIMAL_SIMILARITY>
					--minmatches <MINIMAL_RESIDUE_MATCH_RATIO>
					"""

import os, sys, glob

# --- end of imports --- #

def load_kipes_results( filename, minsim, minmatches ):
	"""! @brief run everything """
	
	# --- find gene function of each gene --- #
	best_hit_per_gene = {}
	with open( filename, "r" ) as f:
		f.readline()
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				if float( parts[2] ) > best_hit_per_gene[ parts[0] ]['score']:
					best_hit_per_gene[ parts[0] ] = { 'id': parts[1].split('_')[0], 'score': float( parts[2] ) }
			except KeyError:
				best_hit_per_gene.update( { parts[0]: { 'id': parts[1].split('_')[0], 'score': float( parts[2] ) } } )
			line = f.readline()
	
	# --- group genes by gene function --- #
	data = {}
	with open( filename, "r" ) as f:
		f.readline()
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if float( parts[2] ) > minsim and float( parts[3] ) > minmatches:	#ensure minimal similarity
				if best_hit_per_gene[ parts[0] ]['id'] == parts[1].split('_')[0]:	#check if this annotated function is the best match
					try:
						data[ parts[1].split('_')[0] ].append( { 'id': parts[0], 'sim': parts[2], 'cons': parts[3] } )
					except KeyError:
						data.update( { parts[1].split('_')[0]: [ { 'id': parts[0], 'sim': parts[2], 'cons': parts[3] } ] } )
			line = f.readline()
	return data


def load_mapping_table( phytozome_to_species_mapping_file ):
	"""! @brief load phytozome to species mapping table """
	
	mapping_table = {}
	with open( phytozome_to_species_mapping_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			mapping_table.update( { parts[0]: parts[1] } )
			line = f.readline()
	return mapping_table


def main( arguments ):
	"""! @brief run everything """
	
	input_folder = arguments[ arguments.index('--in')+1 ]
	output_file  = arguments[ arguments.index('--out')+1 ]
	
	genes_of_interest = arguments[ arguments.index('--genes')+1 ]
	if "," in genes_of_interest:
		genes_of_interest = genes_of_interest.split(',')
	else:
		genes_of_interest = [ genes_of_interest ]
	
	if "--minsim" in arguments:
		minsim = float( arguments[ arguments.index('--minsim')+1 ] )
	else:
		minsim = 0.5	#overall sequence similarity
	if '--minsim' in arguments:
		minmatches = float( arguments[ arguments.index('--minsim')+1 ] )
	else:
		minmatches = 0.9	#matches concerning functionally relevant amino acids
	
	if "--mapping" in arguments:
		phytozome_to_species_mapping_file = arguments[ arguments.index('--mapping')+1 ]
		mapping_table = load_mapping_table( phytozome_to_species_mapping_file )
	else:
		mapping_table = {}
	
	kipes_result_files = glob.glob( input_folder + "*/summary.txt" )
	
	with open( output_file, "w" ) as out:
		out.write( "\t".join( [ "Species", "Gene", "ID", "Similarity", "ConservedResidues" ] ) + "\n" )
		for filename in kipes_result_files:
			kipes_results = load_kipes_results( filename, minsim, minmatches )
			try:
				spec = mapping_table[ filename.split('/')[-2] ]
			except KeyError:
				spec = filename.split('/')[-2]
			for gene in genes_of_interest:
				try:
					hits = kipes_results[ gene ]
					for hit in hits:
						out.write( "\t".join( [ spec, gene, hit['id'], hit['sim'], hit['cons'] ] ) + "\n" )
				except KeyError:
					out.write( "\t".join( [ spec, gene, ".", ".", "." ] ) + "\n" )


if '--in' in sys.argv and '--out' in sys.argv and '--genes' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
