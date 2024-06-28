#!/usr/bin/perl -X
#Description: Converts standard catalog files from JPL and CDMS to tab-delimited files (dropping the tag and QNs)
#Author: Jacob C. Laas
#Date: April 24, 2013

# Converts standard catalog files from JPL and CDMS to tab-delimited files (dropping the tag and QNs)

use Tie::File;
use File::Copy;

# copies all *.cat to *.txt
@files = glob("*.cat");					# glob creates a list of files and this becomes an array named files
foreach $filename (@files) {					# loops through the files array
	while ($filename =~ m/(.+)\.cat/g) {
		copy("$1.cat","$1.txt") or die "Copy failed: $!";
	}
}

# uses Tie::File to do in-place editing of all *.txt files
@files = glob("*.txt");
foreach $file (@files) {
	tie @fullfile, 'Tie::File', $file;		# tie allows for in-place editing by linking an array with the contents of the currently working file
	foreach (@fullfile) { # at each line...					# loops through the file, line-by-line
		s/(.{13})(.{8})(.{8})(.{2})(.{10})(.{3})(.{7})(.+)/$1\t$2\t$3\t$4\t$5\t$6/g; # ...tab-delimits all columns and drops tag and QNs
		s/\s\.(\d)/0\.$1/g;					# ensures all decimals have a leading digit
		s/\s\-\.(\d)/\-0\.$1/g;				# ensures all negative decimals have a leading digit
	}
	untie @fullfile;						# closes the open file and lets the loop move on to the next file
}

