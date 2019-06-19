#! /usr/bin/env julia
# import Pkg; Pkg.add("ArgParse"); Pkg.add("GZip")
import GZip
using Printf
using ArgParse

### A Fastq structure consisting of read entry items,
mutable struct Fastq
    header::String
    sequence::String
    scores::String
end

###Reads an entry of a fastq file and returns a Fastq Object with corresponding entry items: header, sequence, and quality scores
function read_entry(file_ptr::GZip.GZipStream)
    he = String(chomp(readline(file_ptr)))
    se = String(chomp(readline(file_ptr)))
    _  = String(chomp(readline(file_ptr))) # Ignore this line
    sc = String(chomp(readline(file_ptr)))
    return Fastq(he,se,sc)
end

### Takes a DNA sequence and encodes it as a BitArray
function encode_sequence(sequence::String)
    encoder = Dict('A'=> [false,false], 'C'=> [false, true], 'G'=> [true, false], 'T'=> [true,true], 'N' => [false,false])
    encoded_sequence = trues(length(sequence)*2)
    for i = 1:length(sequence)
        base = sequence[i]
        encoded_base = encoder[base]
        encoded_sequence[2*i-1]   = encoded_base[1]
        encoded_sequence[2*i] = encoded_base[2]
    end
    return encoded_sequence::BitArray{1}
end


### Takes a encoded BitArray and returns a DNA sequence as a String
function decode_sequence(encoded_sequence::BitArray{1})
    decoder = Dict([false,false] => 'A', [false, true] => 'C', [true, false] => 'G', [true,true] => 'T')
    decoded_sequence = ""
    for i in 1:2:length(encoded_sequence)
        decoded_sequence = string(decoded_sequence, decoder[encoded_sequence[i:i+1]])
    end
    return decoded_sequence::String
end


### Convert a string of ASCII characters representing quality values into a list of integer scores
function string2int(string::String)
    result = Int[]
    for i in eachindex(string)
        push!(result, Int(string[i]))
    end
    return result
end


### Convert a string of quality scores to an encoded form that matches encoded_sequences
function encode_scores(scores::String)
    encoded_scores = zeros(length(scores)*2)
    for i = 1:length(scores)
        score = scores[i]
        encoded_score = Int(score)/2
        encoded_scores[2*i-1] = encoded_score
        encoded_scores[2*i] = encoded_score
    end
    return encoded_scores
end


function decode_scores(encoded_scores::Vector{Float64})
    decoded_scores = ""
    for i in 1:2:length(encoded_scores)
        decoded_scores = string(decoded_scores, Char(encoded_scores[i]*2))
    end
    return decoded_scores
end


### A quick and dirty function to reverse complement a DNA sequence
function reverse_complement(sequence::String)
    sequence = String(reverse(sequence))
    encoded_sequence = encode_sequence(sequence)
    reversed_encoded = map(!, encoded_sequence)
    return decode_sequence(reversed_encoded)
end


### Converts a string of phred character scores to a vector error probabilities
function phred2prob(phred::String)
    probs = Vector{Float64}(undef, length(phred))
    for i in eachindex(phred)
        probs[i] = phred2prob(phred[i])
    end
    return probs
end
function phred2prob(phred::Char)
    return 10^(-(Int(phred)-33)/10)
end

### Converts a vector of probability scores to a string of phred character scores
function prob2phred(probs::Vector{Float64})
    phreds = Vector{Char}(undef, length(probs))
    for i in eachindex(probs)
        phreds[i] = prob2phred(probs[i])
    end
    return join(phreds)
end
function prob2phred(prob::Float64)
    return Char(Int(round(-10*log(10,prob)+33)))
end

### Converts a string of phred character scores to graphical characters representing quality
###   phred_low   - The character representing the lowest score in a range of scores
###   phred_high  - The character representing the highest score in a range of scores
###   offset      - In an UTF16 character table, the offset from the phred_low character that maps to the 
###                 first character in a new character set (integers encoding '!' => '▃' for example)
###   steps       - # of scores to map to a single character 
function phred2sym(phred::String; phred_low::Char='!', phred_high::Char='K', offset::Int=9601, steps::Int=6)
    syms = Vector{Char}(undef, length(phred))
    low  = Int(phred_low)
    high = Int(phred_high)
    for i in eachindex(phred)
        pcnt = (Int(phred[i])-low)/(high-low)
        pcnt = pcnt >= 1.0 ? 1.0 : pcnt
        new = Int(floor(pcnt*steps))+offset
        syms[i] = Char(new)
        #println(Int(phred[i]), "-", low, "/", high, "*", steps, "+", offset,"=",new,"(",Char(new),")")
    end
    return join(syms)
end


### Align the query sequence against the anchor sequence and return the position and score of the best alignment
### NOTE: Passing sequences in reverse order will work, but will be slower
### Input:
###        anchor - [String] of ACTGs that is the anchor sequence to align against
###         query - [String] of ACTGs that will be aligned to the anchor sequence
###         start - [Int] The position (1-indexed) relative to the anchor string to begin the alignment (default -query)
###          stop - [Int] The position (1-indexed) relative to the anchor string to terminate the alignment (default length(anchor)+length(query))
### Output:
###    best_start - The nucleotide position, relative to anchor where the best query match occurs (can be negative, zero indexed)
###    best_score - The number of matching nucleotides at the best_start position (no threshold is applied, can be a horrible score!)
function align(anchor::String, query::String; start::Int=-length(query), stop::Int=length(anchor))
    if start < -length(query) || start > length(anchor)
        error("The alignment start position was misspecified")
    end
    if stop > length(anchor) || stop < 1
        error("The alignment stop position was misspecified")
    end
    
    # Process and encode input variables
    # Encoded indices
    start = (2*start)
    stop  = (2*stop)

    # Encode sequences
    #scores = Vector{Float64}(undef, cld(stop-start,2))
    scores = fill(0.0, length(query)+length(anchor)+1)
    anchor  =  encode_sequence(anchor)
    query =  encode_sequence(query)

    ## For the desired search range, score the alignments 
    full_bitmask = BitArray(x%2 == 0 for x = 1:length(query))
    for i in start:2:stop
        if i < 0
            overlap = min(i + length(query), length(anchor))
            
            anchor_start = 1
            anchor_stop  = overlap
            
            query_start = -i + 1
            query_stop  = min(query_start+overlap-1, length(query))
            
            bitmask = BitArray(x%2 == 0 for x = 1:(query_stop-query_start)+1)
            
        elseif (i + length(query)) <= length(anchor)
            anchor_start = i + 1
            anchor_stop  = i + length(query)
            
            query_start = 1
            query_stop  = length(query)
            
            bitmask = full_bitmask

        else
            anchor_start = i + 1
            anchor_stop  = length(anchor)
            
            query_start = 1
            query_stop  = length(anchor) - i
            
            bitmask = BitArray(x%2 == 0 for x = 1:(query_stop-query_start)+1)
            
        end
        @inbounds anchor_sub = anchor[anchor_start:anchor_stop]
        @inbounds query_sub = query[query_start:query_stop]

        ## Find matches
        ##                         G   C   T                       A   G   T  
        ## e.g. if anchor_sub = [ 1 0 0 1 1 1 ] and query_sub = [ 0 0 1 0 1 1 ]  (negate)
        ##                      [ 1 0 0 1 1 1 ]      xor        [ 1 1 0 1 0 0 ]  =>  [ 0 1 0 0 1 1 ]
        ##                                                                           [ 0 0 1 0 0 1 ] (right 1 bitshift)
        ##                                                                           [ 0 0 0 0 0 1 ] (anded)
        ##                                                                           [ X 0 X 0 X 1 ] (bitmask)
        ##                                                                                 1         (sum) == total matches
        ## \xor+[tab] for the ⊻ symbol
        partial = (anchor_sub .⊻ .~query_sub)
        align_sub = (partial .& (partial >> 1)) .& bitmask

        # Essentially the percent matched
        if length(query_sub)>0
            score = 2*sum(align_sub) / length(query_sub)
        else
            score = 0
        end
        @inbounds scores[cld((i+length(query)),2)+1] = score
    end

    best_score, best_start = findmax(scores)
    best_start = best_start - cld(length(query),2) - 1
    return (best_start, best_score)
end 


function print_alignment(top::String, bottom::String, offset::Int; width::Int=150)
    flip = false

    if (offset < 1)
        flip = true
        offset = abs(offset)
        temp   = top
        top    = bottom
        bottom = temp
    end


    total_length = offset + length(bottom)

    topover = top[(offset+1):end]
    top     = top[1:offset]
    botover = bottom[1:(end-offset)]
    bottom  = bottom[length(botover):end]


    if (total_length > width)
        s_width = cld(width,3)
        if (length(top) > s_width)
            top = "..." * top[(end-s_width+3):end]
        end
        if (length(topover) > s_width)
            midpoint = cld(s_width,2)
            odd      = length(topover)%2
            println(midpoint)
            topover = topover[1:(midpoint-1-odd)] * "..." * topover[(end-midpoint+1):end]
            botover = botover[1:(midpoint-1-odd)] * "..." * botover[(end-midpoint+1):end]
        end
        if (length(bottom) > s_width)
            bottom = bottom[1:(s_width-3)] * "..."
        end
    end

    if flip
        println(" "^length(top),botover,bottom)
        println(top,topover)
    else
        println(top,topover)
        println(" "^length(top),botover,bottom)
    end
end


## Only called if there is overlap!! (yes, the if/elseifs are a bit redundant), otherwise there are two molecules 
function make_molecule(read1::Fastq, read2::Fastq, offset::Int)
    final_size   = offset+length(read2.sequence)
    final        = Vector{Char}(undef, final_size)
    final_scores = Vector{Char}(undef, final_size)
    for i in 1:final_size
        j = i-offset
        # Before the offset, just take what read1 provides
        if i<=offset
            final[i]        = read1.sequence[i]
            final_scores[i] = read1.scores[i]
        # For the overlap region
        elseif i>=offset && i<=length(read1.sequence)
            r1 = read1.sequence[i]
            s1 = phred2prob(read1.scores[i])
            r2 = read2.sequence[j]
            s2 = phred2prob(read2.scores[j])
            fb = ""
            fs = 0
            # They agree on the base
            if r1==r2
                fb = r1
                fs = prob2phred(s1) # Could also consider multiplying scores or something similar
            # No agreement, read1's score is better
            elseif s1<s2
                fb = r1
                fs = prob2phred(s1)
            # No agreement and read2's score is better
            else
                fb = r2
                fs = prob2phred(s2)
            end
            final[i]        = fb
            final_scores[i] = fs
        # After the end of read1, just take what read2 provides
        else
            final[i]        = read2.sequence[j]
            final_scores[i] = read2.scores[j]
        end
    end
    return (join(final), join(final_scores))
end


### Function to parse command line arguments
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--mode", "-m"
            help = "The display mode: 'none' to just write an output file, 'summary' to write a summary and output file, 'stream' to stream to STDOUT, and 'show' to show parsed alignments"
            arg_type = String  
            default = "summary"
        "--format", "-f"
            help = "The output format ('tsv', 'fastq', 'fastq.gz', 'fasta', 'fasta.gz') of written/streamed results"
            arg_type = String
            default = "fastq"
        "--read1"
            help = "The read 1 fastq(.gz) file"
            arg_type = String
            required = true
        "--read2"
            help = "The read 2 fastq(.gz) file"
            arg_type = String
            required = true
        # CGCAATTCCTTTAGTGGTACCTTTCTATTCTCACTCT for example files
        "--leader", "-L"
            help = "The sequence leading the region of interest"
            arg_type = String
            default = ""
        # CTTTCAACAGTTTCGGCCGAACCTCCACC for example files 
        "--follower", "-F"
            help = "The sequence following the region of interest"
            arg_type = String
            default = ""
        "--align-threshold", "-a"
            help = "Alignment threshold before a read is discarded"
            arg_type = Float64
            default = 0.70
        "--expected-length", "-e"
            help = "The insert length between leader and follower sequences (will discard all others)"
            arg_type = Int
            default = -1
        "--output", "-o"
            help = "The name of the file to write output"
            arg_type = String
            default = ""
    end

    return parse_args(s)
end




## The main function 
function main()
    ### Parse arguments
    parsed_args  = parse_commandline()
    mode         = parsed_args["mode"]
    format       = parsed_args["format"]
    
    # Detect and open file type
    read1_ifh = try 
        GZip.open(parsed_args["read1"])
    catch err
        if isa(err, LoadError)
            open(parsed_args["read1"])
        end
    end
    read2_ifh = try 
        GZip.open(parsed_args["read2"])
    catch err
        if isa(err, LoadError)
            open(parsed_args["read2"])
        end
    end

    leader       = parsed_args["leader"]
    follower_    = reverse_complement(parsed_args["follower"])
    alignment_threshold = parsed_args["align-threshold"]
    expected_length     = parsed_args["expected-length"]
    output       = parsed_args["output"]
    output_ofh   = nothing

    # Check some of the arguments
    if !(mode in ["none", "summary", "stream", "show"])
        print("Please choose 'none', 'summary', 'stream' or 'show' for '--mode'")
        exit()
    end
    if !(format in ["tsv", "fastq", "fastq.gz", "fasta", "fasta.gz"])
        print("Please choose 'tsv', 'fastq', 'fastq.gz', 'fasta', or 'fasta.gz' for '--format'")
        exit()
    end
    # Make sure a file is specified if the mode is 'none' or 'summary'
    if ((mode=="none" || mode=="summary") && isempty(output))
        print("Please provide a file to write results with '--output' if either 'none', 'summary', or 'show' is specified as the mode")
        exit()
    end

    # Open the output file if needed and write a tsv output header if tsv format
    if (mode in ["none", "summary", "show"]) 
        if (format=="fastq.gz" ||format=="fasta.gz")
            output_ofh = GZip.open(output, "w")
        else
            output_ofh = open(output, "w")
            if format=="tsv"
                write(output_ofh, "Insert\tHeader\tRead_Alignment_Score\tAlignment_Qualities")
                !isempty(leader) ? write(output_ofh, "\tLeader_Alignment_Score") : nothing
                !isempty(follower_) ? write(output_ofh, "\tFollower_Alignment_Score") : nothing
                write(output_ofh, "\n")
            end
        end
    end


    # Summary statistics
    entry = 0
    successful = 0
    uil_err = 0
    pra_err = 0
    pla_err = 0
    pfa_err = 0
    zil_err = 0

    while !eof(read1_ifh)
        entry += 1
        mode=="show" ? @printf("Entry: %-6d ", entry) : nothing
        # Show a progress indicator if nothing else is requested
        if (mode=="summary" && (entry%1000==0))
            prog = ['⣀','⡄','⠆','⠃','⠉','⠘','⠰','⢠']
            print(stderr, "Processing ", prog[Int(entry/1000)%(length(prog))+1], "\r")
            flush(stderr)
        end

        # Get the next read entry and parse
        read1 = read_entry(read1_ifh)
        read2 = read_entry(read2_ifh)

        # Check for empty strings (gzip doesn't handle eof properly...)
        if isempty(read1.sequence) || isempty(read2.sequence)
            mode=="show" ? println("\e[1m\e[38;2;255;0;0;249m!\033[0m empty entry") : nothing
            continue
        end
        # Reverse complement read2's sequence and scores
        read2.sequence = reverse_complement(read2.sequence)
        read2.scores = reverse(read2.scores)

        # Align the two reads (quality align is almost 2x slower)
        #(alignment_offset, alignment_score) = quality_align(read1.sequence, read2.sequence)
        (alignment_offset, alignment_score) = align(read1.sequence, read2.sequence)
        
        # If the read alignment doesn't goes well, continue
        if alignment_score <= alignment_threshold
            mode=="show" ? @printf("\e[1m\e[38;2;255;0;0;249m!\033[0m poor read alignment     (%.2f < %.2f)\n", alignment_score, alignment_threshold) : nothing
            pra_err += 1
            continue
        end

        if leader != ""
            # Get alignment of leader sequence
            (leader_start, leader_score) = align(read1.sequence, leader, start=0, stop=length(read1.sequence))

            if leader_score <= alignment_threshold
                mode=="show" ? @printf("\e[1m\e[38;2;255;0;0;249m!\033[0m poor leader alignment   (%.2f < %.2f)\n", leader_score, alignment_threshold) : nothing
                pla_err += 1
                continue
            end
        else
            leader_start = 0
            leader_score = 0
        end

        if follower_ != ""
            # Get alignment of follower sequence
            (follower_start, follower_score) = align(read2.sequence, follower_, start=0, stop=length(read2.sequence))
            if follower_score <= alignment_threshold
                mode=="show" ? @printf("\e[1m\e[38;2;255;0;0;249m!\033[0m poor follower alignment (%.2f < %.2f)\n", follower_score, alignment_threshold) : nothing
                pfa_err += 1
                continue
            end
        else
            follower_start = length(read2.sequence)
            follower_score = 0
        end

        ## Get the insert coordinates
        insert_start = (leader_start+length(leader))
        insert_end   = (alignment_offset+follower_start)

        # make the final molecule
        molecule1 = ""
        molecule1_scores = ""
        molecule2 = ""
        molecule2_scores = ""
        if alignment_score >= alignment_threshold
            (molecule1, molecule1_scores) = make_molecule(read1, read2, alignment_offset)
        else
            molecule1        = read1.sequence
            molecule1_scores = read1.scores
            molecule2        = read2.sequence
            molecule2_scores = read2.scores
        end
        final_molecule = molecule1 * molecule2
        final_scores   = molecule1_scores * molecule2_scores
        insert = final_molecule[insert_start+1:insert_end]
        scores = final_scores[insert_start+1:insert_end]

        if (expected_length!=-1 && expected_length!=length(insert))
            mode=="show" ? @printf("\e[1m\e[38;2;255;0;0;249m!\033[0m unexpected insert size  (%4d)\n", length(insert)) : nothing
            uil_err += 1
            continue
        end
        if (length(insert)==0)
            mode=="show" ? @printf("\e[1m\e[38;2;255;0;0;249m!\033[0m 0 bp insert length\n") : nothing
            zil_err += 1
            continue
        end

        if mode=="show"
            println("\n\e[1m\e[38;2;0;255;0;249m>>>\033[0m")
            if leader_start >= 0 && leader != ""
                println(repeat("~", leader_start), leader)
            end
            println(phred2sym(read1.scores))
            println(read1.sequence)
            ## Bottom strand extends before top strand, need to truncate
            if alignment_offset < 0
                ss = phred2sym(read2.scores)
                println("...", read2.sequence[(abs(alignment_offset)+4):end])
                print("...")
                ## Cant subset unicode character strings, need to manually subset! xP
                for (i,j) in enumerate(eachindex(ss))
                    if i>=(abs(alignment_offset)+4)
                        print(ss[j])
                    end
                end
                println()
            else
                println(repeat(" ", alignment_offset), read2.sequence)
                println(repeat(" ", alignment_offset), phred2sym(read2.scores))
            end
            if follower_ != ""
                ## Bottom strand extends before top strand, need to truncate
                if alignment_offset + follower_start > 0
                    print(repeat(" ", alignment_offset+follower_start), follower_)
                end
                if length(read2.sequence)-follower_start-length(follower_) >= 0
                    print(repeat("~", length(read2.sequence)-follower_start-length(follower_)))
                end
            end
            print("\n\n")
            println(repeat(" ",insert_start), phred2sym(scores))
            if (insert_start>17)
                println("Extracted insert:", repeat(" ", insert_start-17), insert)
            else
                println(repeat(" ", insert_start), insert)
            end
            @printf("Length: %4d",length(insert))
            @printf("\nOffset: %4d",insert_start)
            @printf("\n\nAlignment score: %1.2f\n", alignment_score)
            if leader != ""
                @printf("Leader score:    %1.2f\n", leader_score)
            end
            if follower_ != ""
                @printf("Follower score:  %1.2f\n", follower_score)
            end

            println("\e[1m\e[38;2;0;255;0;249m<<<\033[0m")
        end
        
        # Write output
        line = ""
        if format=="tsv"
            line = @sprintf("%s\t%s\t%1.3f\t%s", insert, read1.header, alignment_score, join(scores))
            !isempty(leader) ? line *= @sprintf("\t%1.3f", leader_score) : nothing
            !isempty(follower_) ? line *= @sprintf("\t%1.3f", follower_score) : nothing
            line *= "\n"
        end
        if format=="fastq" || format=="fastq.gz"
            line = @sprintf("%s  %1.3f\n", read1.header, alignment_score)
            line *= insert
            line *= "\n+\n"
            line *= @sprintf("%s\n", join(scores))
        end
        if format=="fasta" || format=="fasta.gz"
            line = @sprintf("\n>%s  %1.3f\n", read1.header, alignment_score)
            line *= insert
        end
        if mode=="stream"
            print(line)
        elseif mode=="none" || mode=="summary" || mode=="show"
            write(output_ofh, line)
        end

        successful += 1
    end

    ### Display summary results if requested
    if (mode in ["none", "summary", "show"])
        println()
        println("\n🧬  Summary")
        @printf("Total processed:     \t%8d\n", entry)
        @printf("# errors:            \t%8d\n", (entry-successful))
        @printf("# successful:        \t%8d\n", successful)
        @printf("Successfully parsed: \t%8.2f%%\n", 100*Float64(successful)/Float64(entry))

        println("\n🧬  Errors")
        @printf("Read alignment:     \t%8d\n", pra_err)
        if leader != ""
            @printf("Leader alignment:   \t%8d\n", pla_err)
        end
        if follower_ != ""
            @printf("Follower alignment: \t%8d\n", pfa_err)
        end
        @printf("Insert size:        \t%8d\n", uil_err)
        @printf("No insert:          \t%8d\n\n", zil_err)
        @printf(stderr, "\rWrote \"%s\"\n", parsed_args["output"])
    end

    close(read1_ifh)
    close(read2_ifh)
    (mode in  ["none", "summary", "show"]) ? close(output_ofh) : nothing
end

main()

