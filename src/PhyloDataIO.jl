module PhyloDataIO

using BioSequences
using FASTX

export read_PhyloData, read_phylip, read_fasta, write_phylip, write_fasta, read_taxa_char_table, read_nexus, write_PhyloData, write_nexus

function detect_sequence_type(seq::String)::String
    aa_only = Set(['E','F','I','L','P','Q','Z','J','O','U'])
    for char in uppercase(seq)
        if char in aa_only
            return  "protein"
        end
        if isdigit(char)  || char == '+'
        return "standard"
        end
    end
    return "dna"
end

function read_PhyloData(filename::String)::Dict{String,String}
    header =nothing
    sequences = Dict{String,String}()

    open(filename, "r") do fh
        header = strip(first(eachline(fh)))
    end
    if header === nothing
        @error "File $filename is empty"
        exit(1)
    end
    if startswith(header, ">")
         sequences= read_fasta(filename)
    elseif startswith(header, "#")
         sequences= read_nexus(filename)
    elseif occursin(r"^\d+\s+\d+", header)
         sequences= read_phylip(filename)
    else
         sequences= read_taxa_char_table(filename)
    end
    return sequences
end    

function write_PhyloData(seqs::Dict{String,String},filename::String,format::String)
    if format == "n"
        write_nexus(filename, seqs)
    elseif format == "p"
        write_phylip(filename, seqs)
    else
        write_fasta(filename, seqs)
    end
end



function read_taxa_char_table(filename::String)::Dict{String,String}
    seqs= Dict{String,String}()
    open(filename, "r") do fh
        for line in eachline(fh)
            line = strip(line)
            if isempty(line)
                continue
            end
            current_line = split(line, r"\s+", limit=2)
            seqs[String(current_line[1])] = String(current_line[2])
        end
    end
    return seqs
end

function read_phylip(filename::String)::Dict{String,String}    
    seqs = Dict{String,String}()
  
    open(filename, "r") do fh
        for line in eachline(fh)
            line = strip(line)
            if isempty(line)
                continue
            end
            if occursin(r"^\d+\s+\d+$", line)
                continue
            end
            current_line = split(line, r"\s+", limit=2)
            seqs[String(current_line[1])] = String(current_line[2])
        end
    end
    return seqs
end

function read_fasta(filename::String)::Dict{String,String}
    sequences = Dict{String, String}()  
    FASTAReader(open(filename, "r")) do reader
        for record in reader
            sequences[identifier(record)] = sequence(String, record)
        end
    end
    return sequences
end

function read_nexus(filename::String)::Dict{String,String}
    seqs = Dict{String,String}()
    open(filename, "r") do fh
        line_counter=0
        for line in eachline(fh)
             line = strip(line)
            if isempty(line)
                continue
            end
            if occursin(r"^matrix"i, line)
                line_counter += 1
                continue
            end
            if line_counter == 1 && startswith(line, ";")
                break
            end
            if line_counter == 0
                continue
            end
            current_line = split(line, r"\s+", limit=2)
            seqs[String(current_line[1])] = String(current_line[2])
        end
    end
    return seqs
end

function write_phylip(filename::String, seqs::Dict{String,String})
    taxa= sort(collect(keys(seqs)))
    ntaxa= length(taxa)
    nsites= length(seqs[taxa[1]])
    open(filename, "w") do fh
        println(fh, "$ntaxa $nsites")
        for taxon in taxa
            print(fh, "$taxon  ")
            println(fh, seqs[taxon])
        end
    end
end

function write_nexus(filename::String, seqs::Dict{String,String})
    taxa= sort(collect(keys(seqs)))
    ntaxa= length(taxa)
    nsites= length(seqs[taxa[1]])
    seqtype= detect_sequence_type(seqs[taxa[1]])
    open(filename, "w") do fh
        println(fh, "#Nexus")
        println(fh, "Begin data;")
        println(fh, "Dimension ntax=$ntaxa nchar=$nsites;")
        println(fh, "format datatype=$seqtype missing=? gap=-;")
        println(fh, "Matrix")
        for taxon in taxa
            print(fh, "$taxon  ")
            println(fh, seqs[taxon])
        end
        println(fh, ";")
        println(fh, "end;")
    end
end

function write_fasta(filename::String, seqs::Dict{String,String})
    taxa= sort(collect(keys(seqs)))
    open(filename, "w") do fh
        for taxon in taxa
            println(fh, ">$taxon")
            println(fh, seqs[taxon])
        end
    end
end


end