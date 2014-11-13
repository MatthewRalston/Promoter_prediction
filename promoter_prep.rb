#!/usr/bin/env ruby


require 'bio'

#INDIR="motifs"
#OUTDIR="pspm"
i,o=ARGV[0..1]
INDIR,OUTDIR=i.chomp("/"),o.chomp("/")

# This function takes a hashmap of fasta sequences and returns a PSPM
def fasta_to_pspm(hashmap)
  temp = []; array = []
  hashmap.each {|key,value| temp.push(value.split(''))}; num = temp.size
  temp[0].size.times do |x|
    c=0;a=0;g=0;t=0
    temp.each do |liszt|
      case
      when (liszt[x] == "A" or liszt[x] == "a")
        a+=1
      when (liszt[x] == "C" or liszt[x] == "c")
        c+=1
      when (liszt[x] == "G" or liszt[x] == "g")
        g+=1
      when (liszt[x] == "T" or liszt[x] == "t")
        t+=1
      end
    end
    array.push([(a/num.to_f).round(12),(c/num.to_f).round(12),(g/num.to_f).round(12),(t/num.to_f).round(12)])
  end
  return array
end


# This function prints a pspm to mast format
def pspmFile(outdir,name,pspm)
  File.open(outdir+"/"+name+".pspm",'w') do |file|
    file.puts("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.345 C 0.155 G 0.155 T 0.345  \n")
    file.puts("MOTIF #{name}\nletter-probability matrix: alength= 4 w= #{pspm.size} nsites= #{pspm.size} \n")
    pspm.each {|liszt| file.puts(" #{liszt.join("  ")}")}
  end
end


def main
  Dir.glob(INDIR+"/*.txt") do |infile|
    motif=infile.split("/")[1].split(".")[0]
    fastas={}
    Bio::FastaFormat.open(infile).each_entry{|f| fastas[f.entry_id]=f.seq}
    pspm=fasta_to_pspm(fastas)
    pspmFile(OUTDIR,motif,pspm)
  end
end

main
