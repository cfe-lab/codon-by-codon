#!/opt/local/Library/Frameworks/Python.framework/Versions/2.7/bin/python

# Start html output.
print ("Content-type: text/html")
print 
print ( "<html><head>" )
print ( "<title>Results</title>" ) 
print ( "</head><body>" )

import re
import sys  
sys.path.append("/Users/B_Team_iMac/Sites/cgi-bin/python_dependencies/")  # Add the path to openpyxl, (excel files.) (And other web dependencies.)
sys.path.append("/Users/B_Team_iMac/anaconda/pkgs/scipy-0.14.0-np19py27_0/lib/python2.7/site-packages/")  # Add the path to scipy.
sys.path.append("/Users/B_Team_iMac/anaconda/pkgs/numpy-1.9.0-py27_0/lib/python2.7/site-packages/")  # Add the path to numpy.

from scipy import stats  # For the p value calculation function.

import sequence_utils
##import format_utils

import math_utils
import mailer
import cgi


##### Get website input.


form = cgi.FieldStorage()  # Get form data from the website.

# Assign form data to variables.
protein_sequences = [ tuple(e.split('\t')) for e in form.getvalue("functionProtein").replace('\r', '\n').replace('\n\n', '\n').split('\n') ]  # Turn this into a list of tuples. -> (decimal_value, protein_sequence)
min_count = form.getvalue("minCount")
analysis_id = form.getvalue("analysisID")
email_address_string = form.getvalue("emailAddress")

##print ( "Got data.", min_count, analysis_id, email_address_string, protein_sequences )


##### Make sure data is acceptable.  ( Validate data. )


if not math_utils.is_string_int(min_count):
	print ( "<br><br><b><r style=\"color: red;\">Error:</r> Min count isn't an integer;</b> consider removing decimals or changing the value." )
	sys.exit(0)
else:
	min_count = int(min_count)

if not re.match(r"[^@]+@[^@]+\.[^@]+", email_address_string):
	print "<br><br><b><r style=\"color: red;\">Error:</r> Your email address (<em>{}</em>) is missing necessary characters,</b> please re-check its spelling.".format(email_address_string)
	sys.exit(0)

# TODO: Check the protein sequences. (length [DONE], proper characters, formatting, etc...)
# TODO: Add warning if sequence contains mixtures.

# Check if all sequences are the correct length and find said length.
sequence_length = len(protein_sequences[0][1])  # Init the length to be the length of the first protein sequence.
for tuple in protein_sequences:
	if len(tuple[1]) != sequence_length:
		print "<br><br><b><r style=\"color: red;\">Error:</r> All sequences are not the same length,</b> please re-check their formatting."
		sys.exit(0)

print ( "<br><br>Data has been checked.<br><br>" )


##### Class and function definitions.


class ProteinData:
	'''
	This class holds the data for a single protein character in its column column.
	'''

	def __init__(self, char, decimal_list):
		self.protein_char = char
		self.decimal_list = decimal_list

	# This function returns how many times this protein characters appears by getting the data
	# from the decimal value list's length.
	def get_occurrances(self):
		return len(self.decimal_list)

class ColumnData:
	'''
	This class holds the ProteinData for everything in this column.
	'''
	
	# nw_protein_list is a list of ProteinData classes, while with_protein_data is just one.
	def __init__(self, w_protein_dat, nw_protein_list):
		self.with_protein_data = w_protein_dat
		self.not_with_protein_data_list = nw_protein_list

	def get_with(self):
		pass

	def get_not_with(self, index):
		pass
	
	# This function gets all the decimal lists from not with and combines them.
	def get_not_with_all_decimal_lists(self):
		sum_decimal_list = []
		for data in self.not_with_protein_data_list:
			sum_decimal_list += data.decimal_list
		return sum_decimal_list

class ColumnOutput:
	
	def __init__(self):
		self.coord = -1
		self.w_amino = '!'
		self.w_median = -1.0
		self.nw_median = -1.0
		self.w_count = -1
		self.nw_count = -1
		self.p_value = -1.0


# This function extracts and calculates the needed data from the containers and gives it to the output_matrix.
# Note: These containers are of the type -> [[occurrence count, (protein char, [decimal, ...])], ...]
def _run_seq_test(position, w_container, nw_container):
	w_amino = w_container[0][1][0]
	
	w_decimal_list = w_container[0][1][1]  # Save the decimal list for later.
	w_median = math_utils.median( w_decimal_list )
	
	# Find the median of all decimal values.
	nw_decimal_list = []
	for item in nw_container:  # Iterate all decimal lists and add them together.
		nw_decimal_list += item[1][1]
	print( nw_decimal_list )
	nw_median = math_utils.median( nw_decimal_list )

	w_count = w_container[0][0]
	
	# Find the sum of all occurance counts.
	nw_count = 0
	for item in nw_container:  # Iterate all occurrence counts and add them together.
		nw_count += item[0]
	
	p_value = stats.kruskal(w_decimal_list, nw_decimal_list)[1]  # Use scipy to find the p value.
 	
	# Make sure this is the correct position (because 0 start -> 1 start) .
	output_matrix.append( { 'coord': position+1, 'w_amino': w_amino, 'w_median': w_median, 'nw_median': nw_median, 'w_count': w_count, 'nw_count': nw_count, 'p_value<br>': p_value } ) 


##### Run tests on the given sequences. 


#sequence_count = len(protein_sequences)  # How many total sequences there are.

output_matrix = []  # This matrix holds the data that will be added to the excel file. This is a list of dictionaries, so... matrix?

# Iterate over each item in the sequences.
for column_index in range(sequence_length):
	# Create and fill a list that contains all proteins at the current position and their decimal value.
	codon_list = [] # -> [ (char, decimal), ... ]
	for seq_index in range(len(protein_sequences)):
		codon_list.append( (protein_sequences[seq_index][1][column_index], protein_sequences[seq_index][0]) )  # list of -> (character, decimal)
	
	# Count how many of each protein there is.
	data_dict {}  # This is a dict of ProteinData classes with char as the key. -> char : ProteinData
	##codon_dict = {}
	sequence_count = 0  # This holds how many valid characters were counted.  ( for checking validity )
	for tuple in codon_list:
		char = tuple[0]
		decimal = tuple[1]
		if not char in sequence_utils.protein_mixture_list:  # Make sure that invalid characters (mixtures) are not counted in this step.
			sequence_count += 1
			if char in data_dict:
				data_dict[ char ].decimal_list += [ decimal ]
				##data_dict[tuple[0]][0] += 1
				##data_dict[tuple[0]] += [tuple[1]]  # dict with value of -> [occurrence count, decimal, ...]
			else:
				dat = ProteinData( char, [decimal] )
				data_dict[ char ] = dat
				##codon_dict[tuple[0]] = [1, tuple[1]]  # -> [occurrence count, decimal, ...]
	
	# Find the protein (or proteins) with the most ocurrences.	
	if len(codon_dict) > 1:
		most_occurrences = []  # -> [occurrence count, (protein char, decimal, ...), ...]
		for key, occurrences in codon_dict.iteritems():			
			if most_occurrences == []: # Case: First iteration, init list.
				most_occurrences = [ occurrences[0], (key, occurrences[1:]) ]
			else:
				# Compare the current occurrence value with the top values.
				if occurrences[0] > most_occurrences[0]:
					most_occurrences = [ occurrences[0], (key, occurrences[1:]) ]  # -> [occurrence count, (protein char, decimal, ...), ...]
				elif occurrences[0] == most_occurrences[0]:
					most_occurrences += [ (key, occurrences[1:]) ]
	else:
		##print ( 'IGNORE <br>' ) 
		continue  # Case: all characters equal, ignore.
	
	# Check if sequence is above min_count.
	if most_occurrences[0] < min_count or sequence_count - most_occurrences[0] < min_count:
		##print ( 'IGNORE <br>' ) 
		continue  # Case: with or without are smaller than min_count, ignore.
		
	#print( str(codon_list) + "<br>" + str(codon_dict) + "<br>" + str(most_occurrences) + "<br>----------<br>" )
	
	# If this is not 2 then with and not with are tied. (equal)
	if len(most_occurrences) <= 2:
		# Fit with and not-with into different containers that hold information about them.
		# Note: These containers are of the type -> [[occurrence count, (protein char, [decimal, ...])],...]
		w_container = [most_occurrences]  # With.
	
		codon_dict.pop(most_occurrences[1][0])
		nw_container = [[ v[0], (k, list( v[1:] )) ] for k, v in codon_dict.iteritems()]  # Not-with
	
		#print( str(w_container) + " /1/<br>" + str(nw_container) + " /2/<br><br>" )

		_run_seq_test(column_index, w_container, nw_container)  # Run test on containers.
	else:  #Case: do two analysis.
		#print( most_occurrences )
		for index in range( 1, len(most_occurrences) ):
			# Fit with and not-with into different containers that hold information about them.
			# Note: These containers are of the type -> [[occurrence count, (protein char, [decimal, ...])], ...]
			w_container = [[most_occurrences[0], most_occurrences[index]]]  # With.
			
			codon_dict.pop(w_container[0][1][0])  # Temp pop.
			nw_container = [[ v[0], (k, list( v[1:] )) ] for k, v in codon_dict.iteritems()]  # Not-with

			#print( str(w_container) + " /1/<br>" + str(nw_container) + " /2/<br><br>" )
			
			_run_seq_test(column_index, w_container, nw_container)  # Run test on containers.
			codon_dict[w_container[0][1][0]] = [w_container[0][0]] + w_container[0][1][1]  # Put back on.
		
	#print ( '<br>' )
	

##### Output the immediate results to the webpage.


print ( output_matrix )


##### Create an xlsx file.


'''
from openpyxl import Workbook
import openpyxl
'''


##### Send an email with the xlsx file in it.


print ( "</body></html>" )  # Complete the html output.


