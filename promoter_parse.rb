#!/usr/bin/env ruby
=begin
------------------------------------------------------------------------------
--                                                                          --
--                                MATTHEW RALSTON                           --
--                                                                          --
--                       P R O M O T E R _ P A R S E . R B                  --
--                                                                          --
------------------------------------------------------------------------------
-- TITLE                                                                    --
--                                                                          --
--  Fall 2014                                                               --
--                                                                          --
------------------------------------------------------------------------------
-- This file is designed to parse output from Mast promoter predictions.    --
-- The pairing of bipartite motifs requires the paired transcription factor --
-- to have two motif pspms used for mast, and that the results are stored as--
-- Factor_motif1 and Factor_motif2 or Factor_1 and Factor_2
--                                                                          --
------------------------------------------------------------------------------
=end

################################################
#
#               R E Q U I R E
#
################################################

require 'parallel'



################################################
#
#               U S E R    V A R I A B L E S
#
################################################
i,o,cores,mindist,maxdist=ARGV
MINDIST,MAXDIST,CORES=mindist.to_i,maxdist.to_i,cores.to_i
P=p.to_f
INDIR,OUTDIR=i.chomp("/"),o.chomp("/")
RESULTS=Dir[INDIR+"/*.txt"].map{|x| x.split("/")[1]}
pairs=RESULTS.map{|x| x.split("_")}.reject{|x| x.size == 1}
PAIRED_FACTORS=pairs.map{|x| x[0]}.uniq
MOTIFS=pairs.map{|x| x[1]}.uniq # e.g. _upstream | _downstream
SINGLES=RESULTS-pairs.map{|x| x.join("_")}



################################################
#
#               S U B - R O U T I N E S
#
################################################

################################################

################################################

def process_pairs
  Parallel.each(PAIRED_FACTORS, :in_processes => CORES) do |factor|
    motif1=[]; motif2=[]
    File.open(INDIR+"/"+factor+"_"+MOTIFS[0],'r').each {|line| motif1 << line.chomp.split unless line[0]=="#"}
    File.open(INDIR+"/"+factor+"_"+MOTIFS[1],'r').each {|line| motif2 << line.chomp.split unless line[0]=="#"}
    motif1.map!{|x| x[1]=x[1][0]; x[2]=x[2].to_i; x[3]=x[3].to_i;x[-1]=x[-1].to_f;x}.sort!{|x,y| x[2]<=>y[2]}
    motif2.map!{|x| x[1]=x[1][0]; x[2]=x[2].to_i; x[3]=x[3].to_i;x[-1]=x[-1].to_f;x}.sort!{|x,y| x[2]<=>y[2]}
    matches=find_matches(motif1,motif2,factor)
    File.open(OUTDIR+"/"+factor+".results.txt",'w') {|file| matches.each {|x| file.puts(x.join("\t"))}}
  end
end

def find_matches(motif1,motif2,name)
  matches=[]
  motif1.each do |d|
    s,e=d[2..3]
    strand=d[1]
    chrom=d[0]
    results=[]
    if strand == "+"
      i=(0...motif2.size).to_a.bsearch{|x| motif2[x][2] - e >= MINDIST}
      while i && motif2[i][2] - e <= MAXDIST && motif2[i][2] - e >= MINDIST
        results << motif2[i]
        i+=1
      end
    else
      i=(0...motif2.size).to_a.bsearch{|x| s - motif2[x][3] <= MAXDIST}
      while i && s - motif2[i][3] >= MINDIST && s - motif2[i][3] <= MAXDIST
        results << motif2[i]
        i+=1
      end
    end
    results.select!{|x| x[1]==strand && x[0]==chrom}
    results.sort!{|x,y| x[-1]<=>y[-1]} if results[0]
    if results[0]
      matches << [chrom,"MAST\tpromoter",[s,results[0][2]].min,[e,results[0][3]].max,".",strand,".\tID #{name}; Upstream_pval #{d[-1]}; Downstream_pval #{results[0][-1]};"]
    end
  end
  matches
end

def process_single
  Parallel.each(SINGLES, :in_processes => CORES) do |factor|
    motif=[]
    File.open(INDIR+"/"+factor,'r').each {|line| motif << line.chomp.split unless line[0]=="#"}
    motif.map!{|x| x[1]=x[1][0];[x[0],"MAST\tPromoter",x[2],x[3],".",x[1],".\tID #{factor.split(".")[0]}_motif; pval #{x[-1]};"]}
    File.open(OUTDIR+"/"+factor,'w') {|file| motif.each {|rec| file.puts(rec.join("\t"))}}
  end
end

def main
  process_pairs
  process_single
end
#*****************************************************************************#
################################################
#
#-----------------------------------------------
#
#                   M A I N
#-----------------------------------------------
#
################################################

main





##########################  E O F   ############################################
