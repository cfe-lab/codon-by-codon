#!/opt/local/Library/Frameworks/Python.framework/Versions/2.7/bin/python

# Start html output.
print ("Content-type: text/html \n")  # Note: This is 2.7 and 3.6 compatible.
print ( "<html><head>" )
print ( "<title>Results</title>" ) 
print ( "</head><body>" )

import re
import cgi
import sys

sys.path.append("/Users/B_Team_iMac/Sites/cgi-bin/python_dependencies/libraries/")  # Add the path to openpyxl, (excel files.)
from openpyxl import Workbook
import openpyxl

sys.path.append("/Users/B_Team_iMac/anaconda/pkgs/scipy-0.14.0-np19py27_0/lib/python2.7/site-packages/")  # Add the path to scipy.
sys.path.append("/Users/B_Team_iMac/anaconda/pkgs/numpy-1.9.0-py27_0/lib/python2.7/site-packages/")  # Add the path to numpy (scipy dependency).
from scipy import stats  # For the p value calculation function.

sys.path.append("/Users/B_Team_iMac/Sites/cgi-bin/python_dependencies/util_scripts/") 
import sequence_utils
##import format_utils
import math_utils
import mailer


##### Get website input.


form = cgi.FieldStorage()  # Get form data from the website.

# Assign form data to variables.
protein_sequences = [ tuple(e.split('\t')) for e in form.getvalue("functionProtein").replace('\r', '\n').replace("\n\n", '\n').replace(' ', '').split('\n') ]  # Turn this into a list of tuples. -> (decimal_value, protein_sequence)
min_count = form.getvalue("minCount")
analysis_id = form.getvalue("analysisID")
email_address_string = form.getvalue("emailAddress")
desc_string = form.getvalue("analysisID")

##print ( "Got data.", min_count, analysis_id, email_address_string, protein_sequences )


##### Make sure data is acceptable (validate data) and raise any warnings.


if not math_utils.is_string_int(min_count):
	print ( "<br><br><b><r style=\"color: red;\">Error:</r> Min count isn't an integer;</b> consider removing decimals or changing the value." )
	sys.exit(0)
else:
	min_count = int(min_count)

if email_address_string == "":
	print ( "<br><br><b><r style=\"color: red;\">Error:</r> Did not run analysis. You didn't enter an email address.</b> please enter an email address." )
	sys.exit(0)

elif not re.match(r"[^@]+@[^@]+\.[^@]+", email_address_string):
	print ( "<br><br><b><r style=\"color: red;\">Error:</r> Did not run analysis. Your email address (<em>{}</em>) is missing necessary characters,</b> please re-check its spelling.".format(email_address_string) )
	sys.exit(0)

# Check if all sequences are the correct length and find said length.
sequence_length = len(protein_sequences[0][1])  # Init the length to be the length of the first protein sequence.
for tuple in protein_sequences:
	if len(tuple[1]) != sequence_length:
		print ( "<br><br><b><r style=\"color: red;\">Error:</r> All sequences are not the same length,</b> please re-check their formatting." )
		sys.exit(0)

# Check if all sequences contain valid characters.
send_error = False
char_messages = ""
row_number = 0
for tuple in protein_sequences:
	row_number += 1	
	index = 0
	for char in tuple[1]:
		if (char in sequence_utils.valid_protein_character_list) == False:
			send_error = True
			char_messages += "<br><b>{}</b> was found at position {} of row {}.".format(char, index, row_number)  # Report any invalid characters.
		index += 1

# Print error message.
if send_error == True:
	print ( "<br><br><b><r style=\"color: red;\">Error:</r> Some invalid characters have been found,</b> please remove them to run the analysis." )	
	print ( char_messages )
	sys.exit(0)

print ( "Data raised no errors.<br><br>" )

# Gives a warning if the sequence contains mixture characters.
found_warning = False
for tuple in protein_sequences:
	if found_warning == True:  # This is the exit condition.
		break
	
	for char in tuple[1]:
		if (char in sequence_utils.protein_mixture_list) == True:
			found_warning = True
			print ( "<b><r style=\"color: orange;\">Warning:</r></b> Some mixture characters (X and/or -) have been found in this analysis and will be ignored.<br><br>" )
			break


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
	def __init__(self, column_pos, w_protein_dat, nw_protein_list):
		self.with_protein_data = w_protein_dat
		self.not_with_protein_data_list = nw_protein_list
		self.position = column_pos

	def get_with(self):
		return self.with_protein_data

	def get_not_with_list(self):
		return self.not_with_protein_data_list
	
	# This function gets all the decimal lists from not-with and combines them.
	def get_not_with_decimal_list(self):
		sum_decimal_list = []
		for data in self.not_with_protein_data_list:
			sum_decimal_list += data.decimal_list
		return sum_decimal_list

	# This function gets all the occurrance values from not-with and combines them into one number.
	def get_not_with_occurrances(self):
		sum_occurrance_count = 0
		for data in self.not_with_protein_data_list:
			sum_occurrance_count += data.get_occurrances()
		return sum_occurrance_count
	
	
class ColumnOutput:
	'''
	ColumnOutput contains the information that will be formatted into a .xlsx file.
	'''
	
	def __init__(self, coord, with_amino, with_median, notwith_median, with_count, notwith_count, p_value):
		self.coord = coord
		self.w_amino = with_amino
		self.w_median = with_median
		self.nw_median = notwith_median
		self.w_count = with_count
		self.nw_count = notwith_count
		self.p_value = p_value
	
	# This function formats all the information and returns it as a list. ( row )
	def get_formatted_row(self):
		return [ self.coord, self.w_amino, self.w_median, self.nw_median, self.w_count, self.nw_count, self.p_value ]
	
	
# This function extracts and calculates the needed data from the ColumnData class.
def _run_seq_test( column_data ):
	w_amino = column_data.get_with().protein_char  # Get Protein char for with.
	
	w_decimal_list = column_data.get_with().decimal_list  # Save the decimal list for later.
	print ( w_decimal_list )
	w_median = math_utils.median( w_decimal_list )
	
	# Find the median of all not-with decimal values.
	nw_decimal_list = column_data.get_not_with_decimal_list()
	nw_median = math_utils.median( nw_decimal_list )
	
	w_count = column_data.get_with().get_occurrances()
	
	# Find the sum of all not-with occurance counts.
	nw_count = column_data.get_not_with_occurrances()
	
	p_value = math_utils.round_to_sig_figs(stats.kruskal(w_decimal_list, nw_decimal_list)[1], 4)  # Use scipy to find the p value.
 	
	# Put all this information into the ColumnOutput class.
	output_column = ColumnOutput( column_data.position+1, w_amino, w_median, nw_median, w_count, nw_count, p_value )  # Make sure this is the correct position (because 0 start -> 1 start) .
	
	# Send the output column to the output matrix.
	output_matrix.append( output_column ) 
	
	
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
	data_dict = {}  # This is a dict of ProteinData classes with char as the key. -> { char : ProteinData }
	sequence_count = 0  # This holds how many valid characters were counted.  ( For checking validity. )
	for tuple in codon_list:
		char = tuple[0]
		decimal = float(tuple[1])

		if not char in sequence_utils.protein_mixture_list:  # Make sure that invalid characters (mixtures) are not counted in this step.
			sequence_count += 1

			if char in data_dict:
				data_dict[ char ].decimal_list += [ decimal ]  # Populate the decimal list.
			else:
				data_dict[ char ] = ProteinData( char, [decimal] )  # Init the ProteinData class.
	
	# Find the protein (or proteins) with the most ocurrences.	
	if len(data_dict) > 1:  # Case: there are at least two different proteins in this column.
		most_occurrences = [ data_dict.values()[0] ]  # Init the list with the first item. 

		for data in data_dict.values()[1:]:  # Iterate all data classes except the first one.			
			# Compare the current occurrence value with the top values.
			if data.get_occurrances() > most_occurrences[0].get_occurrances():
				most_occurrences = [ data ]  
			elif data.get_occurrances() == most_occurrences[0].get_occurrances():  # Case: both protein data classes have equal occurance values, include both. 
				most_occurrences += [ data ]
	else:
		continue  # Case: all characters equal, ignore.

	# Check if sequence is above min_count.
	if most_occurrences[0].get_occurrances() < min_count or sequence_count - most_occurrences[0].get_occurrances() < min_count:
		continue  # Case: with or without are smaller than min_count, ignore.
		
	#print( str(codon_list) + "<br>" + str(codon_dict) + "<br>" + str(most_occurrences) + "<br>----------<br>" )
	
	# If there are more than one protein data classes with the same occurance value, do the analysis for each protein. 
	if len(most_occurrences) <= 1:
		# Put with and not with into the same container.
		data_dict.pop( most_occurrences[0].protein_char )  # Data is not needed again, data is never re-assigned to the dict.
		
		column_data = ColumnData( column_index, most_occurrences[0], data_dict.values() )
		_run_seq_test( column_data )  # Run test on the data.
	
	else:  #Case: do two analysis.
		for index in range( 0, len(most_occurrences) ):
			# Put with and not with into the same container.
			data_dict.pop( most_occurrences[index].protein_char )  # Temp pop.

			column_data = ColumnData( column_index, most_occurrences[index], data_dict.values() )
			_run_seq_test( column_data )  # Run test on the data.
			
			data_dict[most_occurrences[index].protein_char] = most_occurrences[index]  # Add the item back to the dictionary. 
	
	
##### Create an xlsx file.


XLSX_FILENAME = "codon_by_codon_data"

wb = Workbook()  # Create a new workbook.
ws = wb.active  # Create a new page. (worksheet [ws])
ws.title = "Data"  # Page title

# Create the title row information (key).
ws.append( ["Coord", "Amino", "Median(With)", "Median(Without)", "N(With)", "N(Without)", "Kruskal-wallis p"] )

# Add rows to the document.
for item in sorted( output_matrix, key=lambda x: x.p_value ):
	ws.append( item.get_formatted_row() )

# Save a string version of the excel workbook and send it to the file builder.
file_text = openpyxl.writer.excel.save_virtual_workbook(wb)
xlsx_file = mailer.create_file( XLSX_FILENAME, 'xlsx', file_text )


##### Send an email with the xlsx file in it.


print ( "--"*35 )
print ( "<br><br>" )

# Add the body to the message and send it.
end_message = "This is an automatically generated email, please do not respond."
msg_body = "The included .xlsx file ({}.xlsx) contains the requested {}. \n\nAnalysis description: {} \n\n{}".format(XLSX_FILENAME, "codon analysis data", desc_string, end_message)

if mailer.send_sfu_email("codon_analysis", email_address_string, "Codon by codon analysis: {}".format( desc_string ), msg_body, [xlsx_file]) == 0:
	print ( "An email has been sent to <b>{}</b> with a full table of results. <br>Make sure <b>{}</b> is spelled correctly.".format(email_address_string, email_address_string) )

# Check if email is formatted correctly.
if not re.match(r"[^@]+@[^@]+\.[^@]+", email_address_string):
	print ( "Your email address (<b>{}</b>) is likely spelled incorrectly, please re-check its spelling.".format(email_address_string) )
	
print ( "<br><br>" )
print ( "--"*35 )
	
print ( "<br><br> python version: " + sys.version )  # Print version number.

print ( "</body></html>" )  # Complete the html output.


