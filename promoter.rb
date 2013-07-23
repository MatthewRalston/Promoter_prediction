=begin
------------------------------------------------------------------------------
--                                                                          --
--                                 MATTEW RALSTON                           --
--                                                                          --
--                               P R O M O T E R . R B                      --
--                                                                          --
------------------------------------------------------------------------------
-- TITLE                                                                    --
--                                                                          --
--  Summer 2013                                                             --
--                                                                          --
------------------------------------------------------------------------------
-- This file is designed to be a pipeline for the conversion of a fasta     --
-- format file of individual sequences thought or known to represent 				--
-- a promoter consensus sequence and converts them first into a position    --
-- specific probability matrix (PSPM). Then the PSPM can be used by the MAST--
-- program to scan input sequences for the promoter motif. These results    --
-- are then parsed to yield the p-value and position of the matches. 				--
-- The parser can take inputs from both bipartite and single motifs.				--
-- WARNING: The input fasta format sequences must be of 
--																																					--
-- VARIABLES																																--
-- 1. The first variables that must be set is the number and name of the 		--
-- single and bipartite motif containing fasta format files.								--
-- This is handled by environment variables named 'BIPARTITE' and 'SINGLE'  --
-- 2. The background frequencies of the nucleotide resiudes is handled by   --
-- an environment variable named 'BACKGROUND_FREQUENCIES.'
--                                                                          --
------------------------------------------------------------------------------
=end

################################################
#
#               R E Q U I R E
#
################################################

require 'bio'



################################################
#
#               U S E R    V A R I A B L E S
#
################################################

SINGLE = ['spo0A.txt']
# The variable below describes pairs of motifs representing a bipartite consensus sequence. The first motif in the pair should be 'upstream'.
# The motifs should be named like 'motifname_anythingatall'
BIPARTITE = [ ['siga_one.txt', 'siga_two.txt'] ]
# The variable below describes the background frequences of nucleotides in order.
# => 										   A			C				G			T
BACKGROUND_FREQUENCIES = [0.345, 0.155, 0.155, 0.345]
TARGETS = infile('targets.txt')
MAXDIST = 30
MINDIST = 10
CUTOFF = 0.05


################################################
#
#               S U B - R O U T I N E S
#
################################################

################################################

################################################

def infile(filein)
  hashmap = {}; c = 0; keys = []
  File.open(filein, 'r').each(sep="\n") do |line|
    if c == 0
      keys = line.split("\t")
    else
      line=line.split("\t")
      (puts "duplicate entry: #{line[0]}"; exit) if hashmap[line[0]]
      hashmap[line[0]] = {}
      begin
        for i in 1..keys.size()
          hashmap[line[0]][keys[i]] = line[i]
        end
      end
    end
    c += 1
  end
  return hashmap
end

# This function reads a fasta format file, checking if the length of each sequence is uniform.
#
def fastaread(filein)
  hashmap = {}; len = 0; c = 0
  ffile = Bio::FastaFormat.open(filein)
  ffile.each_entry do |f|
		len = f.seq.size if c==0
		f.seq
		(puts "one or more of the sequences in #{filein} contained Uracil residues."; exit) if (f.seq.include?("U") or f.seq.include?('u'))
    f.seq.size == len ? hashmap[f.entry_id] = f.seq.upcase : (puts "#{filein} contained sequence(s) that had non-uniform length. Please reformat!"; exit)
		c=1
  end
  ffile.close()
  return hashmap
end


#This function takes a hashmap of fasta sequences and returns a position specific probability matrix of the nucleotides per position.
# The 2D array is arranged such that rows represent the position in the MSA, the columns represent the four letters of the DNA alphabet,
# and the values in the ith row and jth column represent the probability of the 'j' residue in the 'i' position based on the MSA.
#
def fasta_to_pspm(hashmap)
	temp = []; array = []
	hashmap.each {|key,value| temp.push(value.split(''))}; num = temp.size
	temp[0].size.times do |x|
		c=0;a=0;g=0;t=0
		temp.each do |liszt|
			case liszt[x]
			when "A"
				a+=1
			when "C"
				c+=1
			when "G"
				g+=1
			when "T"
				t+=1
			end
		end
		array.push([a/num,c/num,g/num,t/num])
	end
	return array
end


# This function takes a hash of 2D arrays (containing PSPMs) and prints each in MAST format.
# Then each PSPM is used to run the 'mast' program from the MEME suite. The raw results in both
# HTML and flat text format can be found in directories according to the name of the motif.
# The results are then gathered by a call to the mastparser function, the relative coordinates are converted to
# strand specific genomic coordinates and returned to this function in a hash.
# If the motif is part of a bipartite motif, we add the hash to a special hash 'bipartite' and continue.
# Otherwise, the results are printed to a file named 'motifname.res.txt' through a call to the 'outfile' funciton.
# The bipartite motifs are handled by the 'bipartite_merger' function.
#
def mast(hashmap)
	bipartite = {}
	hashmap.each do |key,value|
		`mkdir #{key}`
		File.open('pspm.txt', 'w') do |file|
			file.puts("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA #{BACKGROUND_FREQUENCIES[0]} C #{BACKGROUND_FREQUENCIES[1]} G #{BACKGROUND_FREQUENCIES[2]} T #{BACKGROUND_FREQUENCIES[3]}  \n")
			file.puts("MOTIF #{key}\nletter-probability matrix: alength= 4 w= #{value.size} nsites= #{value.size} \n")
			value.each {|liszt| file.puts(" #{liszt.join("  ")}")}
		end
		`mast pspm.txt targets.txt -o #{key}/results -hit_list -norc -mt 0.01 -w -compe -ev 1000 > #{key}/results.txt`
		`mast pspm.txt targets.txt -o #{key}/results -mt 0.01 -w -compe -ev 1000`
		temp = mastparser("#{key}/results.txt")
		bipart = False
		BIPARTITE.each {|liszt| bipart = True if liszt.include?(key)}
		bipart ? bipartite[key] = temp : outfile(temp, "#{key}.res.txt")
	end
	bipartite_merger(bipartite)
end


# This function parses output from the mast program from the meme suite. Results meeting the significance cutoff are passed to the 'adder' function.
# The resulting hashmap is loaded into the 2D hashmap to be returned. The results are organized such that sequences matching the motif are loaded into the hashmap
# with an ID of the form 'sequencename.n' where n is a zero based natural number, such that no IDs are duplicate.
#
def mastparser(infile)
	c=0; hashmap = {}; name = ''; n = 0
	File.open(infile).each(sep="\n") do |line|
		c < 2 ? c+=1 : (line = line.split; (line[0] == name ? (hashmap["#{line[0]}.#{n.to_s}"] = adder(line); n+=1) : (n=0; hashmap["#{line[0]}.#{n.to_s}"] = adder(line); n=1)) if line[-1].to_f <= CUTOFF)
	end
	return hashmap
end


# This function is the workhorse of the mastparser function. It detects which strand the target is on and converts the relative coordinates to genomic coordinates.
#
def adder(liszt)
	hash = {}
	if TARGETS[liszt[0]]["Strand"] == '+'
		hash["Start"] = TARGETS[line[0]]["Start"].to_i + liszt[2].to_i - 1; hash["End"] = TARGETS[line[0]]["Start"].to_i + liszt[3].to_i - 1; hash["p-value"] = liszt[-1]
	else
		hash["Start"] = TARGETS[line[0]]["End"].to_i - liszt[3].to_i - 1; hash["End"] = TARGETS[line[0]]["End"].to_i + liszt[2].to_i - 1; hash["p-value"] = liszt[-1]
	return hash
end


# This function takes a hashmap for a single motif and prints the information.
#
def outfile(hashmap, filename)
	File.open(filename, 'w') do |file|
		file.puts("ID\tStart\tEnd\tp-value")
		hashmap.each {|key,value|  file.puts("#{key}\t#{value["Start"]}\t#{value["End"]}\t#{value["p-value"]}")}
	end
end

# 
# This function takes a bipartite hash of hashes and merges the two motifs into a single bipartite motif.
# First, we iterate over the pairs of bipartite motifs in the constant list 'BIPARTITE,' initializing the pair of hashes (one, two).
# Next, we sort the keys of these hashes according to their genomic coordinates and initialize two lists (pairs, names)
# Then we begin iterating through one hash, looking for a matching downstream motif within the distance restrictions first with binary search.
# More than one downstream motif may match the distance criteria for a single upstream motif, so these are found by considering the nearest neighbors in the array.
# This function is strand specific.
#
def bipartite_merger(bipartite)
	BIPARTITE.each do |first, second|
		one = bipartite[first]; two = bipartite[second]
		onekeys = one.keys.sort {|x,y| one[x]["Start"] <=> one[y]["Start"]}; plus = []; minus = []
		two.each {|key,value| value["Strand"] == "+" ? plus << key : minus << key}
		plus = plus.sort{|x,y| two[x]["Start"] <=> two[y]["Start"]}
		minus = minus.sort{|x,y| two[y]["Start"] <=> two[y]["Start"]}
		
		pairs = []; names = []
		one.each do |key, value|
			l = -1; r = 1; genekey = key.split(".")[0]; 
			if TARGETS[geneeky]["Strand"] == '+'
				twokeys = plus
				start = Array(0..twokeys.size).bsearch {|x| two[twokeys[x]]["Start"] - value["End"] < MINDIST ? 1 : (two[twokeys[x]]["Start"] - value["End"] > MAXDIST ? -1 : (pairs.push([value, two[twokeys[x]]]); names.push(genekey); return 0))}
				(pairs.push([value, two[twokeys[l]]]); names.push(genekey); l-=1) while (two[twokeys[l]]["Start"] - value["End"] >= MINDIST and two[twokeys[l]]["Start"] - value["End"] <= MAXDIST)
				(pairs.push([value, two[twokeys[r]]]); names.push(genekey); r+=1) while (two[twokeys[r]]["Start"] - value["End"] >= MINDIST and two[twokeys[r]]["Start"] - value["End"] <= MAXDIST)	
			else
				twokeys = minus
				start = Array(0..twokeys.size).bsearch {|x| value["Start"] - two[twokeys[x]]["End"] > MAXDIST ? 1 : (value["Start"] - two[twokeys[x]]["End"] < MINDIST ? -1 : (pairs.push([value, two[twokeys[x]]]); names.push(genekey); return 0))}
				(pairs.push([value, two[twokeys[l]]]); names.push(genekey); l-=1) while (value["Start"] - two[twokeys[l]]["End"] >= MINDIST and value["Start"] - two[twokeys[l]]["End"] <= MAXDIST)
				(pairs.push([value, two[twokeys[r]]]); names.push(genekey); r+=1) while (value["Start"] - two[twokeys[r]]["End"] >= MINDIST and value["Start"] - two[twokeys[r]]["End"] <= MAXDIST)	
		end
		File.open("#{first.split("_")[0]}.res.txt", 'w') do |file|
			file.puts("ID\tStart\tEnd\tLeft_p-value\tRight_p-value")
			name = ''; c = 0
			pairs.size.times do |x|
				name == names[x] ? (file.puts("#{name}.#{c}\t#{pairs[x][0]["Start"]}\t#{pairs[x][1]["End"]}\t#{pairs[x][0]["p-value"]}\t#{pairs[x][1]["p-value"]}"); c+=1) : (name = names[x]; c=0; file.puts("#{name}.#{c}\t#{pairs[x][0]["Start"]}\t#{pairs[x][1]["End"]}\t#{pairs[x][0]["p-value"]}\t#{pairs[x][1]["p-value"]}"); c=1)
			end
		end
	end
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







##########################  E O F   ############################################