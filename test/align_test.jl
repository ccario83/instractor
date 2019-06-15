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


#         "GGAGCTGTCGTATTCCAATCAGGTGTGATGCTCGGGGATCCGAATTCCTGGAGCCATCCGCAGTTCGAAAAGCCCCAGAGGCCAGAC"
#                                                                          "CGAAAAG"




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
function align(anchor::String, query::String; start::Int=-length(query), stop::Int=length(anchor), verbose::Bool=false)
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
    scores = fill(0, length(query)+length(anchor)+1)
    anchor  =  encode_sequence(anchor)
    query =  encode_sequence(query)

    verbose && println("Query  length: ", cld(length(query),2), " (",length(query)," encoded)")
    verbose && println("Anchor length: ", cld(length(anchor),2), " (",length(anchor)," encoded)")
    verbose && println("Scores length: ", length(scores))
    verbose && println("Start: ", start)
    verbose && println("Stop:  ", stop)

    ## For the desired search range, score the alignments 
    full_bitmask = BitArray(x%2 == 0 for x = 1:length(query))
    for i in start:2:stop
        verbose &&  println("\n>>>> Index: ",i," (Anchor Nucleotide: ",cld(i,2),") <<<<")
        if i < 0
            verbose &&  println("BEFORE")
            overlap = min(i + length(query), length(anchor))
            
            anchor_start = 1
            anchor_stop  = overlap
            
            query_start = -i + 1
            query_stop  = min(query_start+overlap-1, length(query))
            
            bitmask = BitArray(x%2 == 0 for x = 1:(query_stop-query_start)+1)
            
        elseif (i + length(query)) <= length(anchor)
            verbose &&  println("WITHIN")
            anchor_start = i + 1
            anchor_stop  = i + length(query)
            
            query_start = 1
            query_stop  = length(query)
            
            bitmask = full_bitmask

        else
            verbose &&  println("PAST")
            anchor_start = i + 1
            anchor_stop  = length(anchor)
            
            query_start = 1
            query_stop  = length(anchor) - i
            
            bitmask = BitArray(x%2 == 0 for x = 1:(query_stop-query_start)+1)
            
        end
        @inbounds anchor_sub = anchor[anchor_start:anchor_stop]
        @inbounds query_sub = query[query_start:query_stop]

        verbose &&  println("Query:  ", decode_sequence(query_sub))
        verbose &&  println("Anchor: ", decode_sequence(anchor_sub))
        verbose &&  println("Bitmask:", bitmask)

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

        score = sum(align_sub)
        @inbounds scores[cld((i+length(query)),2)+1] = score
    end

    best_score, best_start = findmax(scores)
    best_start = best_start - cld(length(query),2) - 1


    anchor = decode_sequence(anchor)
    query  = decode_sequence(query)
    println("Best score: ", best_score, " @ index ", best_start, " (Max score: ",length(query),")\n")
    println(repeat(" ", min(length(query)+best_start,length(query))), repeat(" ", max(0,best_start)), query)
    println(repeat(" ", max(0,length(query))), anchor)
    lookup = anchor[max(1,best_start+1):min(length(anchor),best_start+length(query))]
    verbose && println("\nIndexed from anchor: ", lookup)
    verbose && println("Matched query: ", lookup==query)
    println("="^50)
    return (best_start, best_score)
end 


println("="^50)
println("Shorter anchor, left overhang")
align("GGAAAAA", "AAAAGCCCCGATAGA")
println("Shorter anchor, no overhang")
align("CCCC", "AAAAGCCCCGATAGA")
println("Shorter anchor, right overhang")
align("AAAAGCCCCGATAGA", "GACAGATTT")

println("\n\n")

println("="^50)
println("Equal length, left overhang")
align("GATAGAGGGTTACCA", "AAAAGCCCCGATAGA")
println("Equal length, no overhang")
align("GGAAGCCCCGATAGT", "AAAAGCCCCGATAGA")
println("Equal length, rigth overhang")
align("GGTTAAAAAGCGTTT", "AAAAGCCCCGATAGA")

println("\n\n")

println("="^50)
println("Shorter query, left overhang")
align("GATAGAGGGTTACCA", "CGATGATA")
println("Shorter length, no overhang")
align("GGAAGCCCCGATAGT", "CCCC")
println("Shorter query, rigth overhang")
align("GGTTAAAAAGCGTTT", "TTTCCCC")



println("\n\n")

println("="^50)
println("No match")
align("GACAGAGGGCCACCA", "TTTTTTT")
println("Poor match")
align("GACAGAGGGCCACCA", "TTGTTT")