#! /usr/bin/perl -w

$usage = "samtools view [BAM file] | perl $0 [Read1 fastq] [Read2 fastq] [output read1] [output read2]\n";
die $usage unless (@ARGV == 4);

while (<STDIN>)
{
  unless (/^\@/)
  {
    @fields = split /\t/, $_;
    unless (($fields[1] == 77) or ($fields[1] == 141))
    {
      $remove{$fields[0]} = 1;
    }
  }
}
close INFILE;

print $ARGV[0], "\t", $ARGV[1], "\t"; 
print scalar(keys %remove), " read pairs removed\n";
for $i(0..1)
{
  open INFILE, $ARGV[$i] or die;
  open OUTFILE, ">$ARGV[$i+2]" or die;
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
}
