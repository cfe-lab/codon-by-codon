#!/usr/bin/perl

    local ($buffer, @pairs, $pair, $name, $value, %FORM);
    # Read in text
    $ENV{'REQUEST_METHOD'} =~ tr/a-z/A-Z/;
    if ($ENV{'REQUEST_METHOD'} eq "POST")
    {
        read(STDIN, $buffer, $ENV{'CONTENT_LENGTH'});
    }else {
        $buffer = $ENV{'QUERY_STRING'};
    }
    # Split information into name/value pairs
    @pairs = split(/&/, $buffer);
    foreach $pair (@pairs)
    {
        ($name, $value) = split(/=/, $pair);
        $value =~ tr/+/ /;
        $value =~ s/%(..)/pack("C", hex($1))/eg;
        $FORM{$name} = $value;
    }
    $first_name = $FORM{first_name};
    $last_name  = $FORM{last_name};

$| = 1;

print "Content-type:text/html\r\n\r\n";

use Scalar::Util qw(looks_like_number);
use Email::Valid;
use DBI;

print "<html>\n";
print "<body>\n";

my $textareainput = $FORM{'functionProtein'};
my $analysisID = $FORM{'analysisID'};		$analysisID =~ s/\n//g;		$analysisID =~ s/\r//g;
my $email_address = $FORM{'emailAddress'};	$email_address =~ s/\n//g;	$email_address =~ s/\r//g;
my $minCount = $FORM{'minCount'};

print "<h2>Performing codon by codon analysis</h2>\n";

my $proteinLength = -99;

my $emailValidator = (Email::Valid->address($email_address) ? 'yes' : 'no');
if ($emailValidator ne 'yes') { print "<br/>Error detected in user input - the email <b>$email_address</b> is invalid.</body></html>"; exit; }

if (! ($minCount =~ /^[0-9]+$/)) { print "<br/>Error detected in user input - the minimum count <b>$minCount</b> is not a positive integer.</body></html>"; exit; }

print "<p>Analysis description: <b>$analysisID</b>. Min count: <b>$minCount</b></p>\n";

my @rows = split(/\n/, $textareainput);
my $row;

my $lineNumber = 1;


my $db = DBI->connect('DBI:mysql:scratch_space', 'scratch', 'MABimac') || die "Couldn't connect to database: " . DBI->errstr;

# Clear any old contents of the temporary table
my $query = "DELETE FROM codonByCodonTemp";
my $dbHandle = $db->prepare($query);
$dbHandle->execute || die "Couldn't execute SQL";

# Populate the table with the data
foreach $row (@rows) {


	my @numDelimiters = ($row =~ /\t/g);
	my $numDelimiters = @numDelimiters;
	if ($numDelimiters != 1) { print "<br/>Error on input line $lineNumber: Make sure you have 2 fields, separated by a tab!</body></html>"; exit; }

	my ($function, $aminoSeq) = split(/\t/, $row);
	$aminoSeq =~ s/\n//g;
	$aminoSeq =~ s/\r//g;

	# Make sure the function is a number
	if (!looks_like_number($function) ) { print "<br/>Error on input line $lineNumber: <b>$function</b> is not a number!"; exit; }

	# Make sure the amino sequence is actually an amino sequence
	if (! ($aminoSeq =~ /^[ACDEFGHIKLMNPQRSTVWYX\-_]+$/) ) {
		print "<br/>Error on input line $lineNumber: <b>'$aminoSeq'</b> contains invalid characters!";
		exit;
		}

	$query = "INSERT INTO codonByCodonTemp VALUES ('$function', '$aminoSeq')";

	$dbHandle = $db->prepare($query);
	$dbHandle->execute || die "Couldn't execute SQL";

	if ( ($proteinLength != -99) && ($proteinLength != length($aminoSeq))) { print "<br/>Error on input line $lineNumber: sequence lengths not consistent"; exit; }
	$proteinLength = length($aminoSeq);

	$lineNumber++;
	}

# START PROCEDURE FOR CODON BY CODON ANALYSIS HERE
print "<p>The analysis will take a few minutes to execute. You will recieve an email with the final results once this analysis is complete.</p>\n";
print "<p>Assumed protein length: $proteinLength</p>\n";

print "<p><table border='0'><tr>";

my $counter = 0;

my $count = 1;
while ($count <= $proteinLength) {
	my $query =	"SELECT COUNT(*) FROM (SELECT substring(sequence,$count,1) aminos, COUNT(*) " .
			"FROM codonByCodonTemp GROUP BY substring(sequence,$count,1) HAVING aminos REGEXP '[ACDEFGHIKLMNPQRSTVWY]') AS numUniqueAminos";

	my $dbHandle = $db->prepare($query);                    $dbHandle->execute || die "Couldn't execute SQL";
	my $numUniqueAminos = $dbHandle->fetchrow_array();      if($numUniqueAminos < 2) { $count++; next; }

	my $query =	"SELECT substring(sequence,$count,1) aminos, COUNT(*) FROM codonByCodonTemp " .
			"GROUP BY substring(sequence,$count,1) HAVING aminos REGEXP '[ACDEFGHIKLMNPQRSTVWY]'";

	my $dbHandle = $db->prepare($query);
	$dbHandle->execute || die "Couldn't execute SQL";

	while (my @aminoAcids = $dbHandle->fetchrow_array()) {
		my $aminoAcid = $aminoAcids[0];

		# Perform min count filtering first
		my $queryWith = 	"SELECT COUNT(*) FROM codonByCodonTemp " .
					"WHERE substring(sequence,$count,1) LIKE '$aminoAcid' AND substring(sequence,$count,1) REGEXP '[ACDEFGHIKLMNPQRSTVWY]'";

		my $dbHandle2 = $db->prepare($queryWith);	$dbHandle2->execute || die "Couldn't execute SQL";
		my @fields = $dbHandle2->fetchrow_array();	my ($countWith) = @fields;
		if ($countWith < $minCount) { next; }

		my $queryWithout = 	"SELECT COUNT(*) FROM codonByCodonTemp " .
					"WHERE substring(sequence,$count,1) NOT LIKE '$aminoAcid' AND substring(sequence,$count,1) REGEXP '[ACDEFGHIKLMNPQRSTVWY]'";

		$dbHandle2 = $db->prepare($queryWithout);	$dbHandle2->execute || die "Couldn't execute SQL";
		my @fields = $dbHandle2->fetchrow_array();	my ($countWithout) = @fields;
		if ($countWithout < $minCount) { next; }

		if ($counter % 15 == 0) { print "</tr><tr>"; }
		print "<td>$count $aminoAcid</td>\n";
		$counter++;

		my $query =     "(SELECT 0, function FROM codonByCodonTemp " .
				"WHERE substring(sequence,$count,1) NOT LIKE '$aminoAcid' AND substring(sequence,$count,1) REGEXP '[ACDEFGHIKLMNPQRSTVWY]') " .
				"UNION ALL " .
				"(SELECT 1, function  FROM codonByCodonTemp " .
				"WHERE substring(sequence,$count,1) LIKE '$aminoAcid' AND substring(sequence,$count,1) REGEXP '[ACDEFGHIKLMNPQRSTVWY]')";

		$dbHandle2 = $db->prepare($query);
		$dbHandle2->execute || die "Couldn't execute SQL";

		### Generate CSV files for R ###
		my $fileName = $count . "_" . $aminoAcid;
		open(OUTPUT, ">scriptOutput/$fileName.csv");
		print OUTPUT "$aminoAcid,geneFunction\n";
		while (my @fields = $dbHandle2->fetchrow_array()) { print OUTPUT "$fields[0],$fields[1]\n"; }
		close(OUTPUT);

		### Generate R script for each csv file ###
		open(R_SCRIPT, ">scriptOutput/$fileName.rcScript");
		print R_SCRIPT "rcData<-read.csv(\"$fileName.csv\");\n";
		print R_SCRIPT "kruskal.test(geneFunction ~ $aminoAcid, data = rcData);\n";
		print R_SCRIPT "median(rcData\$geneFunction[rcData\$$aminoAcid == 1]);\n";
		print R_SCRIPT "median(rcData\$geneFunction[rcData\$$aminoAcid == 0]);\n";
		print R_SCRIPT "length(rcData\$geneFunction[rcData\$$aminoAcid == 1]);\n";
		print R_SCRIPT "length(rcData\$geneFunction[rcData\$$aminoAcid == 0]);\n";
		close(R_SCRIPT);		

		### Execute R script for each csv file (Results automatically saved) ###
		chdir("scriptOutput");
		system("R CMD BATCH $fileName.rcScript");
		chdir("..");
		
		}

	$count++;
	}

print "</tr></table></p>\n";

print "<p>Compiling results into summary file...\n";

# Parse the R output
open(OUTPUT, ">scriptOutput/summary.out");

chdir("scriptOutput");
my @R_OutputFiles = <*.Rout>;
my $file;
foreach $file (@R_OutputFiles) {
	open(DATA, $file) || die "Couldn't open .Rout file for final summary";

	my @fields = split(/[.]/, $file);
	my $info = $fields[0];
	my @realFields = split(/_/,$info);
	my $coordinate = $realFields[0];
	my $amino = $realFields[1];

	while(my $line = <DATA>) {
		chomp $line;

		if ($line =~ m/p-value/) {
			print OUTPUT "$coordinate\t$amino\t$line";
			$line = <DATA>; $line = <DATA>; $line = <DATA>; chomp $line; print OUTPUT "$line\t";
			$line = <DATA>; $line = <DATA>; chomp $line; print OUTPUT "$line\t";
			$line = <DATA>; $line = <DATA>; chomp $line; print OUTPUT "$line\t";
			$line = <DATA>; $line = <DATA>; chomp $line; print OUTPUT "$line\n";
			}
		}
	}
close(OUTPUT);

# Make a final output file
my $fileName = "summary.out";

open (OUTPUT, ">codon_by_codon_analysis.out");
print OUTPUT "Coord\tAmino\tMedian(With)\tMedian(Without)\tN(With)\tN(Without)\tKruskal-wallis p\n";


open(DATA, $fileName) || die "Couldn't open $fileName";
while(my $line = <DATA>) {
	chomp $line;
	$line =~ s/\[1\]//g;
	$line =~ s/Kruskal-Wallis chi-squared = //g;
	$line =~ s/, df = 1, p-value = /\t/g;
	$line =~ s/[ ]+/\t/g;
	$line =~ s/\t+/\t/g;
	my @fields = split(/\t/,$line);
	my ($coord, $amino, $chi_squared, $p, $medianWith, $medianWithout, $N_with, $N_without) = @fields;
	print OUTPUT "$coord\t$amino\t$medianWith\t$medianWithout\t$N_with\t$N_without\t$p\n";
	}
close DATA;

print " done!</p>\n";

system("rm -rf *.Rout");
system("rm -rf *.rcScript");
system("rm -rf *.csv");
system("rm -rf summary.out");

# Clear the database when done
$query = "DELETE FROM codonByCodonTemp";
$dbHandle = $db->prepare($query);
$dbHandle->execute || die "Couldn't execute SQL";

print "<p>Emailing <b>$email_address</b> ...";

use MIME::Lite;
use Net::SMTP;

my $destination = $email_address;
my $cc = 'zbrumme@sfu.ca';
my $subject = "Codon by codon analysis: $analysisID";
my $mainText = "Your codon by codon analysis job is complete!\n\nAnalysis description: $analysisID\n\nTo read this file, make sure you use the IMPORT function in excel. Do not directly try to open this file by double clicking it!\n\nThis is an automatically generated email, please do not respond.";
my $filePath = "codon_by_codon_analysis.out";
my $fileName = "codon_analysis.csv";
$msg = MIME::Lite->new( From => 'codon.analysis@sfu.ca', To => $destination, Cc => $cc, Subject => $subject, Type => 'multipart/mixed', );
$msg->attach( Type => 'TEXT', Data => $mainText );
$msg->attach( Type => 'application/zip', Path => $filePath, Filename => $fileName, Disposition => 'attachment' );
$msg->send('smtp','mailhost.sfu.ca', Debug=>0 );

print " done! (If you are not in North America it can take mora than an hour for this email to be recieved by your mail server)</p>\n";

print "<h2>Use the import function in excel to open the tab-delimited results.</h2>\n";

print "</body></html>";
