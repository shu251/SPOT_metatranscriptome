#! /usr/bin/perl -w

$usage = "samtools view [BAM file] | perl $0 [Read fastq] [output read]\n";
die $usage unless (@ARGV == 2);

while (<STDIN>)
{
  unless (/^\@/)
  {
    @fields = split /\t/, $_;
    unless ($fields[1] == 4)
    {
      $remove{$fields[0]} = 1;
    }
  }
}
close INFILE;

print $ARGV[0], "\t"; 
print scalar(keys %remove), " reads removed\n";
open INFILE, $ARGV[0] or die;
open OUTFILE, ">$ARGV[1]" or die;
while (<INFILE>)
{
  if (/^\@((\S+).*)/)
  {
    if ((exists $remove{$1}) or (exists $remove{$2})) { $flag = 0; }
    else { $flag = 1; }
  }
  if ($flag)
  { print OUTFILE $_; }
}
close INFILE;
close OUTFILE;
