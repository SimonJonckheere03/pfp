//
//  pfp_algo.cpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//

#include <pfp_algo.hpp>
#include <cctype>
#include <random>

namespace
{

bool
is_structural_sentinel(char c)
{
    return c == vcfbwt::pfp::DOLLAR || c == vcfbwt::pfp::DOLLAR_SEQUENCE || c == vcfbwt::pfp::DOLLAR_PRIME;
}

char
replace_non_acgt(char original, std::minstd_rand& gen)
{
    static constexpr char valid_chars[] = {'A', 'C', 'G', 'T'};
    if (original != 'A' && original != 'C' && original != 'G' && original != 'T')
    {
        std::uniform_int_distribution<std::size_t> distribution(0, 3);
        return valid_chars[distribution(gen)];
    }
    return original;
}

char
replace_non_acgt_with_seed(char original, const std::string& seed, std::size_t& seed_index)
{
    if (original != 'A' && original != 'C' && original != 'G' && original != 'T')
    {
        char replacement = seed[seed_index];
        seed_index = (seed_index + 1) % seed.size();
        return replacement;
    }
    seed_index = 0;
    return original;
}

std::string
normalize_biological_sequence(const std::string& sequence, const vcfbwt::pfp::Params& params)
{
    if (not params.acgt_only)
    {
        return sequence;
    }

    std::string normalized;
    normalized.reserve(sequence.size());

    std::minstd_rand gen(42);
    std::string seed;
    if (params.seed_length > 0)
    {
        seed.reserve(params.seed_length);
        for (std::size_t i = 0; i < params.seed_length; ++i)
        {
            seed.push_back(replace_non_acgt('N', gen));
        }
        gen.seed(42);
    }

    std::size_t seed_index = 0;
    for (char c : sequence)
    {
        char normalized_char = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        if (params.seed_length == 0)
        {
            normalized_char = replace_non_acgt(normalized_char, gen);
        }
        else
        {
            normalized_char = replace_non_acgt_with_seed(normalized_char, seed, seed_index);
        }
        normalized.push_back(normalized_char);
    }

    return normalized;
}

void
write_bigbwt_payload_and_map(const std::string& full_text, const std::string& payload_path, const std::string& map_path)
{
    std::string payload;
    payload.reserve(full_text.size());

    std::vector<unsigned char> packed_bits((full_text.size() + 7) / 8, 0);
    vcfbwt::long_type kept_count = 0;

    for (std::size_t i = 0; i < full_text.size(); ++i)
    {
        if (is_structural_sentinel(full_text[i]))
        {
            continue;
        }
        payload.push_back(full_text[i]);
        packed_bits[i / 8] |= static_cast<unsigned char>(1u << (i % 8));
        ++kept_count;
    }

    std::ofstream payload_out(payload_path, std::ios::binary);
    payload_out.write(payload.data(), payload.size());
    vcfbwt::DiskWrites::update(payload_out.tellp());
    payload_out.close();

    std::ofstream map_out(map_path, std::ios::binary);
    const vcfbwt::long_type full_length = full_text.size();
    map_out.write(reinterpret_cast<const char*>(&full_length), sizeof(full_length));
    map_out.write(reinterpret_cast<const char*>(&kept_count), sizeof(kept_count));
    if (not packed_bits.empty())
    {
        map_out.write(reinterpret_cast<const char*>(packed_bits.data()), packed_bits.size());
    }
    vcfbwt::DiskWrites::update(map_out.tellp());
    map_out.close();
}

std::string
build_reference_segment_text(const std::string& reference, const vcfbwt::pfp::Params& params, bool first)
{
    std::string segment;
    segment.reserve(reference.size() + params.w + 1);

    if (first) { segment.push_back(vcfbwt::pfp::DOLLAR); }
    else
    {
        segment.append(params.w - 1, vcfbwt::pfp::DOLLAR_PRIME);
        segment.push_back(vcfbwt::pfp::DOLLAR_SEQUENCE);
    }

    segment.append(reference);

    segment.append(params.w - 1, vcfbwt::pfp::DOLLAR_PRIME);
    segment.push_back(vcfbwt::pfp::DOLLAR_SEQUENCE);
    return segment;
}

std::string
build_sample_segment_text(const std::string& haplotype_sequence, const vcfbwt::pfp::Params& params, bool is_last_contig)
{
    std::string segment;
    segment.reserve(haplotype_sequence.size() + (2 * params.w));
    segment.append(params.w - 1, vcfbwt::pfp::DOLLAR_PRIME);
    segment.push_back(vcfbwt::pfp::DOLLAR_SEQUENCE);
    segment.append(haplotype_sequence);
    segment.append(params.w - 1, vcfbwt::pfp::DOLLAR_PRIME);
    if (is_last_contig) { segment.append(params.w, vcfbwt::pfp::DOLLAR); }
    else { segment.push_back(vcfbwt::pfp::DOLLAR_SEQUENCE); }
    return segment;
}

size_t
reference_payload_offset(const std::vector<vcfbwt::pfp::ReferenceParse>& references_parse, size_t ref_index)
{
    size_t offset = 0;
    for (size_t i = 0; i < ref_index; ++i)
    {
        offset += references_parse[i].length();
    }
    return offset;
}

void
parse_materialized_segment(const std::string& text, vcfbwt::pfp::Dictionary& dictionary, const vcfbwt::pfp::Params& params,
                           const std::set<vcfbwt::hash_type>& ignored_trigger_strings, std::vector<vcfbwt::hash_type>& parse)
{
    std::string phrase;
    phrase.reserve(text.size());
    vcfbwt::KarpRabinHash kr_hash(params.w);

    for (char c : text)
    {
        phrase.push_back(c);
        if (phrase.size() == params.w) { kr_hash.initialize(phrase); }
        else if (phrase.size() > params.w) { kr_hash.update(phrase[phrase.size() - params.w - 1], phrase.back()); }

        if ((phrase.size() > params.w) and ((kr_hash.get_hash() % params.p) == 0))
        {
            std::string_view ts(&(phrase[phrase.size() - params.w]), params.w);
            if (ts.find(vcfbwt::pfp::DOLLAR) != std::string_view::npos) { continue; }
            vcfbwt::hash_type ts_hash = vcfbwt::KarpRabinHash::string_hash(ts);
            if (ignored_trigger_strings.contains(ts_hash)) { continue; }

            parse.push_back(dictionary.check_and_add(phrase));
            phrase.erase(phrase.begin(), phrase.end() - params.w);
            kr_hash.reset();
            kr_hash.initialize(phrase);
        }
    }

    if (phrase.size() > params.w) { parse.push_back(dictionary.check_and_add(phrase)); }
    else if (phrase.size() < params.w)
    {
        spdlog::error("Materialized segment shorter than the parsing window at finalization");
        std::exit(EXIT_FAILURE);
    }
}

} // namespace

//------------------------------------------------------------------------------

bool
vcfbwt::pfp::Dictionary::contains(const std::string& phrase)
{
    // lock the dictionary
    std::lock_guard<std::mutex> guard(dictionary_mutex);

    hash_type phrase_hash = string_hash(&(phrase[0]), phrase.size());
    const auto& ptr = hash_string_map.find(phrase_hash);
    
    return ((ptr != hash_string_map.end()) and (ptr->second.phrase == phrase));
}

vcfbwt::hash_type
vcfbwt::pfp::Dictionary::add(const std::string& phrase)
{
    // lock the dictionary
    std::lock_guard<std::mutex> guard(dictionary_mutex);

    this->sorted = false;
    
    hash_type phrase_hash = string_hash(&(phrase[0]), phrase.size());
    if (hash_string_map.contains(phrase_hash))
    {
        spdlog::error("Hash collision! Hash already in the dictionary");
        std::exit(EXIT_FAILURE);
    }
    
    DictionaryEntry entry(phrase);
    hash_string_map.insert(std::make_pair(phrase_hash, entry));

    // std::cerr << phrase_hash << " " << entry.phrase << std::endl;

    if (this->size() >= (std::numeric_limits<size_type>::max() - insertions_safe_guard))
    { spdlog::error("Dictionary too big for type {}", typeid(size_type).name()); std::exit(EXIT_FAILURE); }
    
    return phrase_hash;
}

vcfbwt::hash_type
vcfbwt::pfp::Dictionary::check_and_add(const std::string &phrase)
{
    // lock the dictionary
    std::lock_guard<std::mutex> guard(dictionary_mutex);

    // Check if present
    hash_type phrase_hash = string_hash(&(phrase[0]), phrase.size());
    const auto& ptr = hash_string_map.find(phrase_hash);

    if ((ptr != hash_string_map.end()) and (ptr->second.phrase != phrase))
    {
        spdlog::error("Hash collision! Hash already in the dictionary for a different phrase");
        std::exit(EXIT_FAILURE);
    }
    else if (ptr != hash_string_map.end()) { return phrase_hash; }

    this->sorted = false;

    DictionaryEntry entry(phrase);
    hash_string_map.insert(std::make_pair(phrase_hash, entry));

    // std::cerr << phrase_hash << " " << entry.phrase << std::endl;

    if (this->size() >= (std::numeric_limits<size_type>::max() - insertions_safe_guard))
    { spdlog::error("Dictionary too big for type {}", typeid(size_type).name()); std::exit(EXIT_FAILURE); }

    return phrase_hash;
}

vcfbwt::hash_type
vcfbwt::pfp::Dictionary::get(const std::string& phrase) const
{
    return string_hash(&(phrase[0]), phrase.size());
}

void
vcfbwt::pfp::Dictionary::sort()
{
    // lock the dictionary
    std::lock_guard<std::mutex> guard(dictionary_mutex);

    // sort the dictionary
    for (auto& entry : this->hash_string_map)
    {
        this->sorted_phrases.emplace_back(std::ref(entry.second.phrase), entry.first);
    }
    std::sort(sorted_phrases.begin(), sorted_phrases.end(), ref_smaller);
    
    // insert in hashmap
    this->hash_to_ranks.clear();
    for (size_type i = 0; i < sorted_phrases.size(); i++)
    {
        hash_to_ranks.insert(std::make_pair(sorted_phrases[i].second, i + 1)); // 1 based
    }
    
    this->sorted = true;
}

const std::string&
vcfbwt::pfp::Dictionary::sorted_entry_at(std::size_t i)
{
    std::lock_guard<std::mutex> guard(dictionary_mutex);
    if (not this->sorted) { sort(); } return sorted_phrases[i].first.get();
}

vcfbwt::size_type
vcfbwt::pfp::Dictionary::hash_to_rank(hash_type hash)
{
    if (not this->sorted) { sort(); }
    auto it = hash_to_ranks.find(hash);
    if (it != hash_to_ranks.end()) { return it->second; }
    else { spdlog::error("Something went wrong"); std::exit(EXIT_FAILURE); }
}

//------------------------------------------------------------------------------

void
vcfbwt::pfp::ReferenceParse::init(const std::string& reference, bool first)
{
    if (params.ignore_ts_file.size() > 0)
    {
        spdlog::info("Reading trigger string to be ignored from: {}", params.ignore_ts_file);
        std::ifstream ts_file(params.ignore_ts_file);
        if (not ts_file.is_open()) { spdlog::error("ERROR!"); std::exit(EXIT_FAILURE); }
    
        char c = '\5';
        while (c != ENDOFDICT)
        {
            std::string ts = "";
            while ((c = ts_file.get()) != ENDOFWORD and c != ENDOFDICT)
            {
                ts.push_back(c);
            }
            if (ts.size() != 0) { this->to_ignore_ts_hash.insert(KarpRabinHash::string_hash(ts)); }
        }
        ts_file.close();
        spdlog::info("To be ingored trigger strings: {}", this->to_ignore_ts_hash.size());
    }

    const std::string biological_reference = normalize_biological_sequence(reference, this->params);
    this->segment_text = build_reference_segment_text(biological_reference, this->params, first);
    
    std::string phrase;
    spdlog::info("Parsing reference contig " + this->ref_id);
    
    // Karp Robin Hash Function for sliding window
    KarpRabinHash kr_hash(this->params.w);
    
    // Reference as first sample, just one dollar to be compatible with Giovanni's pscan.cpp
    if (first) phrase.append(1, DOLLAR);
    else 
    {
         // Append w-1 dollar prime and 1 dollar sequence at the end, reference as the first sample
        phrase.append(this->params.w - 1, DOLLAR_PRIME);
        phrase.append(1, DOLLAR_SEQUENCE);   
        kr_hash.initialize(phrase);    
    }
    
    for (std::size_t ref_it = 0; ref_it < biological_reference.size(); ref_it++)
    {
        char c = biological_reference[ref_it];
        phrase.push_back(c);
        if (phrase.size() == params.w) { kr_hash.initialize(phrase); }
        else if (phrase.size() > params.w) { kr_hash.update(phrase[phrase.size() - params.w - 1], phrase[phrase.size() - 1]); }
        
        if ((phrase.size() > this->params.w) and ((kr_hash.get_hash() % this->params.p) == 0))
        {
            std::string_view ts(&(phrase[phrase.size() - params.w]), params.w);
            hash_type ts_hash = KarpRabinHash::string_hash(ts);
            if (to_ignore_ts_hash.contains(ts_hash)) { continue; }
            
            hash_type hash = dictionary.check_and_add(phrase);
            
            this->parse.push_back(hash);
            this->trigger_strings_position.push_back(ref_it - this->params.w + 1);
    
            phrase.erase(phrase.begin(), phrase.end() - this->params.w); // Keep the last w chars
            
            kr_hash.reset(); kr_hash.initialize(phrase);
        }
    }
    
    // Last phrase
    if (phrase.size() >= this->params.w)
    {
        // Append w-1 dollar prime and 1 dollar sequence at the end, reference as the first sample
        phrase.append(this->params.w - 1, DOLLAR_PRIME);
        phrase.append(1, DOLLAR_SEQUENCE);
    
        hash_type hash = dictionary.check_and_add(phrase);
    
        this->parse.push_back(hash);
        this->trigger_strings_position.push_back(biological_reference.size() - 1);
    }
    else { spdlog::error("The reference doesn't have w dollar prime at the end!"); std::exit(EXIT_FAILURE); }
}


//------------------------------------------------------------------------------

void
vcfbwt::pfp::ParserVCF::init(const Params& params, const std::string& prefix, std::vector<ReferenceParse>& rp, Dictionary& dict, Dictionary& dict_rev, std::size_t t)
{
    this->w = params.w; this->out_file_prefix = prefix; this->p = params.p; this->tags = t; this->parse_size = 0; this->parse_size_rev = 0;
    if (not ((tags & MAIN) or (tags & WORKER))) { spdlog::error("A parser must be either the main parser or a worker"); std::exit(EXIT_FAILURE); }
    
    if (tags & MAIN)
    {
        this->out_file_name = out_file_prefix + EXT::PARSE;
        this->out_file_name_rev = out_file_prefix + ".rev" + EXT::PARSE;
    }
    if ((tags & MAIN) and params.compute_lifting) { this->out_lift_name = out_file_prefix + EXT::LIFTING; }
    if ((tags & MAIN) and ( params.report_lengths or params.compute_lifting ) ) { this->out_len_name = out_file_prefix + EXT::LENGTHS; }
    if ((tags & MAIN) and params.compress_dictionary) { tags = tags | COMPRESSED; }
    this->tmp_out_file_name = TempFile::getName("parse");
    this->out_file.open(tmp_out_file_name, std::ios::binary);
    if(params.compute_lifting)
    {
        this->tmp_out_lift_name = TempFile::getName("lift");
        this->out_lift.open(tmp_out_lift_name, std::ios::binary);
    }
    if(params.report_lengths or params.compute_lifting)
    {
        this->tmp_out_len_name = TempFile::getName("lengths");
        this->out_len.open(tmp_out_len_name, std::ios::binary);
    }

    this->references_parse = &rp;
    this->dictionary = &dict;
    this->dictionary_rev = &dict_rev;
    
    this->params = params;
}


void
vcfbwt::pfp::ParserVCF::operator()(const vcfbwt::Sample& sample)
{
    this->samples_processed.push_back(sample.id());

    for(auto& contig: sample.contigs)
    {
        if (contig.get_ploidy() <= this->working_genotype)
            continue;
        this->contigs_processed.push_back(std::make_pair(sample.id(), contig.id()));

        const size_t ref_index = contig.get_reference_index();
        ReferenceParse& reference_parse = (*references_parse)[ref_index];
        std::vector<long long int>& tsp = reference_parse.trigger_strings_position;

        Contig::iterator segment_iterator(contig, this->working_genotype);
        std::string haplotype_sequence;
        while (not segment_iterator.end())
        {
            haplotype_sequence.push_back(*segment_iterator);
            ++segment_iterator;
        }
        if (this->params.acgt_only)
        {
            std::string normalized_haplotype_sequence = normalize_biological_sequence(haplotype_sequence, this->params);
            std::string segment_text = build_sample_segment_text(normalized_haplotype_sequence, this->params,
                                                                 contig.last(this->working_genotype));

            std::vector<hash_type> materialized_parse;
            parse_materialized_segment(segment_text, *this->dictionary, this->params,
                                       reference_parse.to_ignore_ts_hash, materialized_parse);

            if (not materialized_parse.empty())
            {
                out_file.write(reinterpret_cast<char*>(materialized_parse.data()),
                               sizeof(hash_type) * materialized_parse.size());
                this->parse_size += materialized_parse.size();
            }

            this->reverse_segments.push_back({std::move(segment_text), &reference_parse.to_ignore_ts_hash});
        }
        else
        {
            this->reverse_segments.push_back({build_sample_segment_text(haplotype_sequence, this->params,
                                                                        contig.last(this->working_genotype)),
                                              &reference_parse.to_ignore_ts_hash});

            // Karp Robin Hash Function for sliding window
            KarpRabinHash kr_hash(this->params.w);
            std::string phrase;
            // Every contig starts with w-1 dollar prime and one dollar seq
            phrase.append(this->w - 1, DOLLAR_PRIME);
            phrase.append(1, DOLLAR_SEQUENCE);
            kr_hash.initialize(phrase);

            std::size_t start_window = 0, end_window = 0;

            Contig::iterator contig_iterator(contig, this->working_genotype);

            while (not contig_iterator.end())
            {
                // Compute where we are on the reference
                std::size_t pos_on_reference = contig_iterator.get_ref_it();
                
                if ( not ((contig_iterator.get_var_it() > 0) and (contig_iterator.prev_variation() > (pos_on_reference - (8 * this->w)))))
                {
                    // Set start postion to the position in the reference parse after the last computed phrase
                    if (params.use_acceleration and ((phrase.size() == this->w) and ((pos_on_reference != 0) and (phrase[0] != DOLLAR_PRIME))))
                    {
                        start_window = end_window;
                        while ((tsp[start_window] + this->w) <= pos_on_reference and (start_window < tsp.size() - 2))
                        { start_window++; }
            
                        // Iterate over the parse up to the next variation
                        while (tsp[end_window + 1] < (long long int)(contig_iterator.next_variation() - (this->w + 1))) { end_window++; }
                        
                        // If the window is not empty
                        if ((start_window < end_window - 1) and (tsp[end_window] > pos_on_reference))
                        {
                            spdlog::debug("------------------------------------------------------------");
                            spdlog::debug("copied from {} to {}", tsp[start_window], tsp[end_window] + this->w);
                            spdlog::debug("next variation: {}", contig_iterator.next_variation());
                            spdlog::debug("skipped phrases: {}", end_window - start_window);
                            
                            // copy from parse[start_window : end_window]
                            out_file.write((char*) &(reference_parse.parse[start_window]), sizeof(hash_type) * (end_window - start_window + 1));
                            this->parse_size += end_window - start_window + 1;
                    
                            // move iterators and re initialize phrase
                            contig_iterator.go_to(tsp[end_window]);
                            phrase.clear();
                            for (std::size_t i = 0; i < this->w; i++) { ++contig_iterator; phrase.push_back(*contig_iterator);}
                            
                            kr_hash.reset(); kr_hash.initialize(phrase);
                            
                            ++contig_iterator;
                            spdlog::debug("New phrase [{}]: {}", phrase.size(), phrase);
                            spdlog::debug("------------------------------------------------------------");
                        }
                    }
                }

                // Next phrase should contain a variation so parse as normal, also if we don't
                // want to use the acceleration we should always end up here
                phrase.push_back(*contig_iterator);
                kr_hash.update(phrase[phrase.size() - params.w - 1], phrase[phrase.size() - 1]);
                ++contig_iterator;
            
                if ((phrase.size() > this->params.w) and ((kr_hash.get_hash() % this->params.p) == 0))
                {
                    std::string_view ts(&(phrase[phrase.size() - params.w]), params.w);
                    hash_type ts_hash = KarpRabinHash::string_hash(ts);
                    if (reference_parse.to_ignore_ts_hash.contains(ts_hash)) { continue; }
                    
                    hash_type hash = this->dictionary->check_and_add(phrase);
                
                    out_file.write((char*) (&hash), sizeof(hash_type)); this->parse_size += 1;
            
                    if (phrase[0] != DOLLAR_PRIME)
                    {
                        spdlog::debug("------------------------------------------------------------");
                        spdlog::debug("Parsed phrase [{}] {}", phrase.size(), phrase);
                        spdlog::debug("------------------------------------------------------------");
                    }
                    
                    phrase.erase(phrase.begin(), phrase.end() - this->w); // Keep the last w chars
            
                    kr_hash.reset(); kr_hash.initialize(phrase);
                }
            }

            assert(contig_iterator.length() == (contig_iterator.get_sam_it() - 1)); // -1 because of the last voi itration

            // Last phrase
            if (phrase.size() >= this->w)
            {
                // Append w dollar prime at the end of each sample, also w DOLLAR if it's the last sample
                phrase.append(this->w - 1, DOLLAR_PRIME);
                if (contig.last(this->working_genotype)) { phrase.append(this->w, DOLLAR); }
                else { phrase.append(1, DOLLAR_SEQUENCE); }

                hash_type hash = this->dictionary->check_and_add(phrase);
                
                out_file.write((char*) (&hash), sizeof(hash_type));   this->parse_size += 1;
            }
            else { spdlog::error("A sample doesn't have w dollar prime at the end!"); std::exit(EXIT_FAILURE); }
        }

        // Build the lifting
        if(params.compute_lifting)
        {
            const size_t genotype = this->working_genotype;
            // Initialize the Lift builder
            lift::Lift_builder lvs_builder(contig.get_reference().size());
            // Iterate throgh all the variations
            for(size_t i = 0; i < contig.variations.size(); ++i)
            {
                const vcfbwt::Variation& variation = contig.get_variation(i);

                const size_t var_genotype = contig.genotypes[i][genotype];
                // Get only variation in the current genotype
                if( var_genotype == 0)
                    continue;

                const int rlen = variation.alt[0].size(); // Length of the reference allele
                const int alen = variation.alt[var_genotype].size(); // Length of the alternate allele
                lvs_builder.set(variation.pos, variation.types[var_genotype], rlen, alen );  
            }

            lift::Lift lift(lvs_builder);
            size_t offset = reference_payload_offset(*references_parse, ref_index);
            out_lift.write((char *)&offset, sizeof(offset)); 
            lift.serialize(out_lift);
        }

        // Reporting contig lengths
        if(params.report_lengths or params.compute_lifting)
        {
            const size_t length = haplotype_sequence.size();
            const std::string contig_name = params.haplotype_name_prefix + sample.id() + "_H" + std::to_string(this->working_genotype + 1) + "_" + contig.id();
            out_len << contig_name << " " << length << std::endl;
        }
    }
}



void
vcfbwt::pfp::ParserVCF::close()
{
    if (closed) return; closed = true;
    
    if ((tags & MAIN) or (tags & WORKER))
    {
        vcfbwt::DiskWrites::update(out_file.tellp()); // Disk Stats
        this->out_file.close();
        if(params.compute_lifting) this->out_lift.close();
        if(params.report_lengths or params.compute_lifting) this->out_len.close();
    }
    
    // Output parse, substitute hash with rank
    if (tags & MAIN)
    {
        // close all the registered workers and merge their dictionaries
        spdlog::info("Main parser: closing all registered workers");
        for (auto worker : registered_workers) { worker.get().close(); }
        
        // Occurrences
        std::vector<size_type> occurrences(this->dictionary->size(), 0);
        
        spdlog::info("Main parser: Replacing hash values with ranks in MAIN, WORKERS and reference, wirting .last ans .sai");
        
        std::string last_file_name = out_file_prefix + EXT::LAST;
        std::ofstream last_file(last_file_name);
    
        std::string sai_file_name = out_file_prefix + EXT::SAI;
        std::ofstream sai_file(sai_file_name);
        
        std::size_t pos_for_sai = 0;
    
        // MAIN mmap file and substitute
        if (this->parse_size != 0)
        {
            std::error_code error;
            mio::mmap_sink rw_mmap = mio::make_mmap_sink(tmp_out_file_name, 0, mio::map_entire_file, error);
            if (error) { spdlog::error(error.message()); std::exit(EXIT_FAILURE); }
    
            for (size_type i = 0; i < (rw_mmap.size() / sizeof(hash_type)); i++)
            {
                hash_type hash;
                std::memcpy(&hash, rw_mmap.data() + (i * sizeof(hash_type)), sizeof(hash_type));
                size_type rank = this->dictionary->hash_to_rank(hash);
                std::memcpy(rw_mmap.data() + (i * sizeof(size_type)), &rank, sizeof(size_type));
                occurrences[rank - 1] += 1;
    
                const std::string& dict_string = this->dictionary->sorted_entry_at(rank - 1);
                last_file.put(dict_string[(dict_string.size() - this->params.w) - 1]);
    
                if (pos_for_sai == 0) { pos_for_sai = dict_string.size() - 1; } // -1 is for the initial $ of the first word
                else { pos_for_sai += dict_string.size() - this->params.w; }
                sai_file.write((char*) &pos_for_sai, IBYTES);
            }
            rw_mmap.unmap();
            truncate_file(tmp_out_file_name, this->parse_size * sizeof(size_type));
        }
        
        // repeat for reference
        for(auto& reference_parse : *(this->references_parse))
        {
            if (not reference_parse.parse.empty())
            {
                for (size_type i = 0; i < reference_parse.parse.size(); i++)
                {
                    hash_type rank = this->dictionary->hash_to_rank(reference_parse.parse[i]);
                    reference_parse.parse[i] = rank;
                    occurrences[rank - 1] += 1;
        
                    const std::string& dict_string = this->dictionary->sorted_entry_at(rank - 1);
                    last_file.put(dict_string[(dict_string.size() - this->params.w) - 1]);
        
                    if (pos_for_sai == 0) { pos_for_sai = dict_string.size() - 1; } // -1 is for the initial $ of the first word
                    else { pos_for_sai += dict_string.size() - this->params.w; }
                    sai_file.write((char*) &pos_for_sai, IBYTES);
                }
            }
        }
        
        // repeat for every worker
        for (auto worker : registered_workers)
        {
            if (worker.get().parse_size != 0)
            {
                std::error_code error;
                mio::mmap_sink rw_mmap = mio::make_mmap_sink(worker.get().tmp_out_file_name, 0, mio::map_entire_file, error);
                if (error) { spdlog::error(error.message()); std::exit(EXIT_FAILURE); }
    
                for (size_type i = 0; i < (rw_mmap.size() / sizeof(hash_type)); i++)
                {
                    hash_type hash;
                    std::memcpy(&hash, rw_mmap.data() + (i * sizeof(hash_type)), sizeof(hash_type));
                    size_type rank = this->dictionary->hash_to_rank(hash);
                    std::memcpy(rw_mmap.data() + (i * sizeof(size_type)), &rank, sizeof(size_type));
                    occurrences[rank - 1] += 1;
    
                    const std::string& dict_string = this->dictionary->sorted_entry_at(rank - 1);
                    last_file.put(dict_string[(dict_string.size() - this->params.w) - 1]);
    
                    if (pos_for_sai == 0) { pos_for_sai = dict_string.size() - 1; } // -1 is for the initial $ of the first word
                    else { pos_for_sai += dict_string.size() - this->params.w; }
                    sai_file.write((char*) &pos_for_sai, IBYTES);
                }
                rw_mmap.unmap();
                truncate_file(worker.get().tmp_out_file_name, worker.get().parse_size * sizeof(size_type));
            }
        }
    
        vcfbwt::DiskWrites::update(last_file.tellp());
        last_file.close();
    
        vcfbwt::DiskWrites::update(sai_file.tellp());
        sai_file.close();
        
        // Merging files together
        spdlog::info("Main parser: concatenating parsings from workers and reference, reference as first");
        std::ofstream merged(out_file_name, std::ios_base::binary);
        
        // Reference
        size_type out_parse_size = 0;
        for(auto& reference_parse : *(this->references_parse))
        {
            out_parse_size += reference_parse.parse.size();
            for (auto& e : reference_parse.parse)
            {
                size_type out_e = e;
                merged.write((char*) &out_e, sizeof(size_type));
            }
        }
        // Main
        if ((this->parse_size) != 0)
        {
            std::ifstream main_parse(this->tmp_out_file_name);
            merged << main_parse.rdbuf();
            out_parse_size += this->parse_size;
        }
        
        // Workers
        for (auto worker : registered_workers)
        {
            if (worker.get().parse_size != 0)
            {
                out_parse_size += worker.get().parse_size;
                std::ifstream worker_parse(worker.get().tmp_out_file_name);
                merged << worker_parse.rdbuf();
            }
        }
        vcfbwt::DiskWrites::update(merged.tellp()); // Disk Stats
        merged.close();

        // Merging the length components
        if (params.report_lengths or params.compute_lifting)
        {
            std::ofstream merged(out_len_name);

            // Reference
            for(auto& reference_parse : *(this->references_parse))
            {
                const size_t length = reference_parse.length();
                const std::string contig_name = reference_parse.id();
                merged << contig_name << " " << length << std::endl;
            }

            // Main
            if ((this->parse_size) != 0)
            {
                std::ifstream main_parse(this->tmp_out_len_name);
                if(not main_parse.is_open())
                    spdlog::error("Cannot open file ", this->tmp_out_len_name);
                merged << main_parse.rdbuf();
                main_parse.close();
            }
            // Workers
            for (auto worker : registered_workers)
            {
                if (worker.get().parse_size != 0)
                {
                    std::ifstream worker_parse(worker.get().tmp_out_len_name);
                    if(not worker_parse.is_open())
                        spdlog::error("Cannot open file ", worker.get().tmp_out_len_name);
                    merged << worker_parse.rdbuf();
                    worker_parse.close();
                }
            }
            vcfbwt::DiskWrites::update(merged.tellp()); // Disk Stats
            merged.close();
        }

        // Merging the lifting components
        if (params.compute_lifting)
        {
            std::ofstream merged(out_lift_name, std::ios_base::binary);

            // if compute_lifting => report_lengths
            std::vector<size_t> onset(1,0);
            std::vector<size_t> lengths;
            size_t u = 0;            
            size_t w = this->params.w;            
            std::vector<std::string> names;

            // Reading the lengths
            std::string tmp_name;
            std::size_t tmp_length;
            std::ifstream in_lidx(out_len_name);
            while (in_lidx >> tmp_name >> tmp_length )
            { 
                if (tmp_name != "")
                {
                    u += tmp_length;
                    names.push_back(tmp_name);
                    lengths.push_back(tmp_length);
                    onset.push_back(u);
                }
            }
            ++u;
            in_lidx.close();
            // Build the seqidx structure
            sdsl::sd_vector_builder builder(u, onset.size());
            for (auto idx : onset)
                builder.set(idx);

            sdsl::sd_vector<> starts(builder);
            // Writing the seqidx on disk
            merged.write((char *)&u, sizeof(u));
            merged.write((char *)&w, sizeof(w));

            starts.serialize(merged);
            sdsl::serialize(names.size(), merged);
            for(size_t i = 0; i < names.size(); ++i)
            {
                sdsl::serialize(names[i].size(), merged);
                merged.write((char *)names[i].data(), names[i].size());
            }

            const size_t reference_count = this->references_parse->size();
            merged.write((char *)&reference_count, sizeof(reference_count));

            // Build the empty liftings for the references
            size_t clen = 0;
            for(size_t i = 0; i < this->references_parse->size(); ++i)
            {
                const size_t len = lengths[i];

                sdsl::bit_vector ibv(len);
                sdsl::bit_vector dbv(len);
                sdsl::bit_vector sbv(len);

                lift::Lift lift(ibv, dbv, sbv); 

                sdsl::serialize(clen, merged);
                lift.serialize(merged);  

                clen += len;       
            }

            // Main
            if ((this->parse_size) != 0)
            {
                std::ifstream main_parse(this->tmp_out_lift_name);
                merged << main_parse.rdbuf();
                main_parse.close();
            }
            // Workers
            for (auto worker : registered_workers)
            {
                if (worker.get().parse_size != 0)
                {
                    out_parse_size += worker.get().parse_size;
                    std::ifstream worker_parse(worker.get().tmp_out_lift_name);
                    merged << worker_parse.rdbuf();
                    worker_parse.close();
                }
            }
            vcfbwt::DiskWrites::update(merged.tellp()); // Disk Stats
            merged.close();
        }


        this->parse_size = out_parse_size;
        
        // Print dicitionary on disk
        if (tags & UNCOMPRESSED)
        {
            spdlog::info("Main parser: writing dictionary to disk NOT COMPRESSED");
            std::string dict_file_name = out_file_prefix + EXT::DICT;
            std::ofstream dict(dict_file_name);
        
            for (size_type i = 0; i < this->dictionary->size(); i++)
            {
                dict.write(this->dictionary->sorted_entry_at(i).c_str(), this->dictionary->sorted_entry_at(i).size());
                dict.put(ENDOFWORD);
            }
        
            dict.put(ENDOFDICT);
            
            vcfbwt::DiskWrites::update(dict.tellp()); // Disk Stats
            dict.close();
        }
        
        if (tags & COMPRESSED)
        {
            spdlog::info("Main parser: writing dictionary to disk COMPRESSED");
            std::ofstream dicz(out_file_prefix + EXT::DICT_COMPRESSED);
            std::ofstream lengths(out_file_prefix + EXT::DICT_COMPRESSED_LENGTHS);

            for (size_type i = 0; i < this->dictionary->size(); i++)
            {
                std::size_t shift = 1; // skip dollar on first phrase
                if (i != 0) { shift = this->w; }
                dicz.write(this->dictionary->sorted_entry_at(i).c_str() + shift,
                           this->dictionary->sorted_entry_at(i).size() - shift);
                int32_t len = this->dictionary->sorted_entry_at(i).size() - shift;
                lengths.write((char*) &len, sizeof(int32_t));
            }
    
            vcfbwt::DiskWrites::update(dicz.tellp()); // Disk Stats
            dicz.close();
    
            vcfbwt::DiskWrites::update(lengths.tellp()); // Disk Stats
            lengths.close();
        }
    
        // Outoput Occurrencies
        if(this->params.compute_occurrences)
        {
            spdlog::info("Main parser: writing occurrences to file");
            std::string occ_file_name = out_file_prefix + EXT::OCC;
            std::ofstream occ(occ_file_name, std::ios::out | std::ios::binary);
            occ.write((char*)&occurrences[0], occurrences.size() * sizeof(size_type));
            
            vcfbwt::DiskWrites::update(occ.tellp()); // Disk Stats
            occ.close();
        }

        spdlog::info("Main parser: building reverse parse from the exact global reverse text");
        std::vector<hash_type> reverse_parse_hashes;
        std::set<hash_type> ignored_trigger_strings_rev;

        for (const auto& reference_parse : *(this->references_parse))
        {
            ignored_trigger_strings_rev.insert(reference_parse.to_ignore_ts_hash.begin(), reference_parse.to_ignore_ts_hash.end());
        }

        for (const auto& segment : this->reverse_segments)
        {
            if (segment.ignored_trigger_strings != nullptr)
            {
                ignored_trigger_strings_rev.insert(segment.ignored_trigger_strings->begin(), segment.ignored_trigger_strings->end());
            }
        }

        for (const auto& worker_ref : registered_workers)
        {
            const auto& worker = worker_ref.get();
            for (const auto& segment : worker.reverse_segments)
            {
                if (segment.ignored_trigger_strings != nullptr)
                {
                    ignored_trigger_strings_rev.insert(segment.ignored_trigger_strings->begin(), segment.ignored_trigger_strings->end());
                }
            }
        }

        std::string forward_text;
        std::ifstream merged_parse_stream(out_file_name, std::ios::binary);
        if (not merged_parse_stream.is_open())
        {
            spdlog::error("Unable to reopen merged forward parse {} while building the reverse parse", out_file_name);
            std::exit(EXIT_FAILURE);
        }

        size_type rank = 0;
        while (merged_parse_stream.read((char*) &rank, sizeof(rank)))
        {
            if ((rank == 0) or (rank > this->dictionary->size()))
            {
                spdlog::error("Invalid rank {} found in merged forward parse while building the reverse parse", rank);
                std::exit(EXIT_FAILURE);
            }

            const std::string& dict_string = this->dictionary->sorted_entry_at(rank - 1);
            forward_text.append(dict_string.data(), dict_string.size() - this->params.w);
        }
        merged_parse_stream.close();
        forward_text.append(this->params.w, DOLLAR);
        if (forward_text.size() < (this->params.w + 1))
        {
            spdlog::error("Forward text is too short to build the reverse parse");
            std::exit(EXIT_FAILURE);
        }
        if (forward_text.front() != DOLLAR)
        {
            spdlog::error("Forward text does not start with the expected leading dollar");
            std::exit(EXIT_FAILURE);
        }
        for (std::size_t i = 0; i < this->params.w; ++i)
        {
            if (forward_text[forward_text.size() - this->params.w + i] != DOLLAR)
            {
                spdlog::error("Forward text does not end with {} trailing dollars", this->params.w);
                std::exit(EXIT_FAILURE);
            }
        }

        std::string reverse_text;
        reverse_text.reserve(forward_text.size());
        reverse_text.push_back(DOLLAR);
        reverse_text.append(forward_text.rbegin() + this->params.w, forward_text.rend() - 1);
        reverse_text.append(this->params.w, DOLLAR);

        write_bigbwt_payload_and_map(
            forward_text,
            out_file_prefix + EXT::BIG_BWT_PAYLOAD,
            out_file_prefix + EXT::BIG_BWT_MAP);
        write_bigbwt_payload_and_map(
            reverse_text,
            out_file_prefix + ".rev" + EXT::BIG_BWT_PAYLOAD,
            out_file_prefix + ".rev" + EXT::BIG_BWT_MAP);

        parse_materialized_segment(reverse_text, *this->dictionary_rev, this->params, ignored_trigger_strings_rev, reverse_parse_hashes);
        this->parse_size_rev = reverse_parse_hashes.size();

        std::vector<size_type> occurrences_rev(this->dictionary_rev->size(), 0);
        std::string last_file_name_rev = out_file_prefix + ".rev" + EXT::LAST;
        std::ofstream last_file_rev(last_file_name_rev);
        std::string sai_file_name_rev = out_file_prefix + ".rev" + EXT::SAI;
        std::ofstream sai_file_rev(sai_file_name_rev);
        std::ofstream merged_rev(out_file_name_rev, std::ios_base::binary);
        std::size_t pos_for_sai_rev = 0;

        for (const hash_type hash : reverse_parse_hashes)
        {
            const size_type rank = this->dictionary_rev->hash_to_rank(hash);
            merged_rev.write((char*) &rank, sizeof(rank));
            occurrences_rev[rank - 1] += 1;

            const std::string& dict_string = this->dictionary_rev->sorted_entry_at(rank - 1);
            last_file_rev.put(dict_string[(dict_string.size() - this->params.w) - 1]);

            if (pos_for_sai_rev == 0) { pos_for_sai_rev = dict_string.size() - 1; }
            else { pos_for_sai_rev += dict_string.size() - this->params.w; }
            sai_file_rev.write((char*) &pos_for_sai_rev, IBYTES);
        }

        vcfbwt::DiskWrites::update(last_file_rev.tellp());
        last_file_rev.close();
        vcfbwt::DiskWrites::update(sai_file_rev.tellp());
        sai_file_rev.close();
        vcfbwt::DiskWrites::update(merged_rev.tellp());
        merged_rev.close();

        if (tags & UNCOMPRESSED)
        {
            spdlog::info("Main parser: writing reverse dictionary to disk NOT COMPRESSED");
            std::string dict_file_name_rev = out_file_prefix + ".rev" + EXT::DICT;
            std::ofstream dict_rev(dict_file_name_rev);

            for (size_type i = 0; i < this->dictionary_rev->size(); i++)
            {
                dict_rev.write(this->dictionary_rev->sorted_entry_at(i).c_str(), this->dictionary_rev->sorted_entry_at(i).size());
                dict_rev.put(ENDOFWORD);
            }

            dict_rev.put(ENDOFDICT);
            vcfbwt::DiskWrites::update(dict_rev.tellp());
            dict_rev.close();
        }

        if(this->params.compute_occurrences)
        {
            spdlog::info("Main parser: writing reverse occurrences to file");
            std::string occ_file_name_rev = out_file_prefix + ".rev" + EXT::OCC;
            std::ofstream occ_rev(occ_file_name_rev, std::ios::out | std::ios::binary);
            occ_rev.write((char*)&occurrences_rev[0], occurrences_rev.size() * sizeof(size_type));
            vcfbwt::DiskWrites::update(occ_rev.tellp());
            occ_rev.close();
        }
    }
}


//------------------------------------------------------------------------------

void
vcfbwt::pfp::ParserFasta::init(const Params& params, const std::string& prefix)
{
    this->params = params;
    this->w = this->params.w;
    this->p = this->params.p;
    this->parse_size = 0;
    this->out_file_prefix = prefix;
    this->out_file_name = prefix + EXT::PARSE;
    this->out_file.open(out_file_name);
    this->out_lift_name = out_file_prefix + EXT::LIFTING;
    this->out_len_name = out_file_prefix + EXT::LENGTHS;

    if(params.compute_lifting)
    {
        init_ref(this->in_file_path, this->w, true);
    }
}

void
vcfbwt::pfp::ParserFasta::init_ref(const std::string& ref_path, const size_t w, bool last)
{
    spdlog::info("Reading reference file: {}", ref_path);
    std::ifstream in_stream(ref_path);
    if (not vcfbwt::is_gzipped(in_stream))
    {
        spdlog::error("Reference file expected gzip compressed");
        std::exit(EXIT_FAILURE);
    }
    
    zstr::istream is(in_stream);
    std::string line;
    
    while (getline(is, line)) 
    { 
        // TODO: Remove the second reference append when we can use multiple ReferenceParses.
        if ( not (line.empty() or line[0] == '>') ) { this->references.back().append(line); this->reference.append(line); } 
        else 
        { 
            // Update lengths
            if( this->ref_sum_lengths.size() > 1 ) this->ref_sum_lengths.push_back(this->ref_sum_lengths.back());
            else this->ref_sum_lengths.push_back(0);
            if(this->references.size()>0) this->ref_sum_lengths.back() += this->references.back().size() + w; // The +w is for the separator that has to be counted in the offset
            // Create the new reference
            std::string name = "";
            for (size_t i = 1; i < line.size(); name.push_back(line[i++]))
                if (isspace(line[i])) break;
            // line.substr(1, line.find(' ') - 1); // -1 to remove the space
            spdlog::info("Read contig {}.", name);
            // TODO: store also the description
            this->references_id.insert(std::make_pair(name, this->references.size())); 
            this->references_name.push_back(name); 
            this->references.push_back(""); 
        }
    }
    if (not last) reference.push_back(pfp::SPECIAL_TYPES::DOLLAR_PRIME);
    spdlog::info("Done reading {}", ref_path);
}

void
vcfbwt::pfp::ParserFasta::operator()()
{
    // Open input file with kseq
    gzFile fp; kseq_t *record;
    fp = gzopen(this->in_file_path.c_str(), "r");
    if (fp == 0)
    {
        spdlog::error("Failed to open input file {}", in_file_path);
        exit(EXIT_FAILURE);
    }
    
    std::string phrase;
    spdlog::info("Parsing sequence");
    
    // Karp Robin Hash Function for sliding window
    KarpRabinHash kr_hash(this->params.w);
    
    // First sequence start with one dollar
    phrase.append(1, DOLLAR);
    
    record = kseq_init(fp);
    while(kseq_read(record) >= 0)
    {
        std::string sequence_name("<error reading sequence name>"), sequence_comment;
        if (record->name.s != NULL) { sequence_name = record->name.s; }
        if (record->comment.s != NULL) { sequence_comment = record->comment.s; }
        this->sequences_processed.push_back(sequence_name + " " + sequence_comment);
        spdlog::debug("Parsed:\t{}", sequence_name + " " + sequence_comment);
        
        // Previous last phrase
        if (phrase[0] != DOLLAR and phrase.size() >= this->params.w)
        {
            // Append w-1 dollar prime, and one dollar seq at the end of each sequence
            phrase.append(this->params.w - 1, DOLLAR_PRIME);
            phrase.append(1, DOLLAR_SEQUENCE);
    
            hash_type hash = this->dictionary.check_and_add(phrase);
     
            out_file.write((char*) (&hash), sizeof(hash_type)); this->parse_size += 1;
    
            // Reset phrase
            phrase.erase();
            phrase.append(this->params.w - 1, DOLLAR_PRIME);
            phrase.append(1, DOLLAR_SEQUENCE);
            kr_hash.reset(); kr_hash.initialize(phrase);
        }
        
        for (std::size_t seq_it = 0; seq_it < record->seq.l; seq_it++)
        {
            char c = record->seq.s[seq_it];
        
            phrase.push_back(c);
            if (phrase.size() == params.w) { kr_hash.initialize(phrase); }
            else if (phrase.size() > params.w) { kr_hash.update(phrase[phrase.size() - params.w - 1], phrase[phrase.size() - 1]); }
        
            if ((phrase.size() > this->params.w) and ((kr_hash.get_hash() % this->params.p) == 0))
            {
                std::string_view ts(&(phrase[phrase.size() - params.w]), params.w);
                hash_type ts_hash = KarpRabinHash::string_hash(ts);
            
                hash_type hash = this->dictionary.check_and_add(phrase);
    
                out_file.write((char*) (&hash), sizeof(hash_type)); this->parse_size += 1;
                
                phrase.erase(phrase.begin(), phrase.end() - this->params.w); // Keep the last w chars
                
                kr_hash.reset(); kr_hash.initialize(phrase);
            }
        }
    }
    
    // Last phrase
    if (phrase.size() > this->params.w)
    {
        // Append w-1 dollar prime, and one dollar seq at the end of each sequence
        phrase.append(this->params.w - 1, DOLLAR_PRIME);
        
        // Append w dollar at the end
        phrase.append(this->params.w, DOLLAR);
        
        hash_type hash = this->dictionary.check_and_add(phrase);
        
        out_file.write((char*) (&hash), sizeof(hash_type)); this->parse_size += 1;
    }
    else { spdlog::error("Missing w DOLLAR at the end!"); std::exit(EXIT_FAILURE); }
    
    kseq_destroy(record);
    gzclose(fp);
}


void
vcfbwt::pfp::ParserFasta::close()
{
    if (closed) return; closed = true;
    
    vcfbwt::DiskWrites::update(out_file.tellp()); // Disk Stats
    this->out_file.close();
    
    // Occurrences
    std::vector<size_type> occurrences(this->dictionary.size(), 0);
    
    spdlog::info("Main parser: Replacing hash values with ranks, writing .last and .sai");
    
    std::string last_file_name = out_file_prefix + EXT::LAST;
    std::ofstream last_file(last_file_name);
    
    std::string sai_file_name = out_file_prefix + EXT::SAI;
    std::ofstream sai_file(sai_file_name);
    
    std::size_t pos_for_sai = 0;
    
    // mmap file and substitute
    if (this->parse_size != 0)
    {
        std::error_code error;
        mio::mmap_sink rw_mmap = mio::make_mmap_sink(out_file_name, 0, mio::map_entire_file, error);
        if (error) { spdlog::error(error.message()); std::exit(EXIT_FAILURE); }
        
        for (size_type i = 0; i < (rw_mmap.size() / sizeof(hash_type)); i++)
        {
            hash_type hash;
            std::memcpy(&hash, rw_mmap.data() + (i * sizeof(hash_type)), sizeof(hash_type));
            size_type rank = this->dictionary.hash_to_rank(hash);
            std::memcpy(rw_mmap.data() + (i * sizeof(size_type)), &rank, sizeof(size_type));
            occurrences[rank - 1] += 1;
    
            const std::string& dict_string = this->dictionary.sorted_entry_at(rank - 1);
            last_file.put(dict_string[(dict_string.size() - this->params.w) - 1]);
    
            if (pos_for_sai == 0) { pos_for_sai = dict_string.size() - 1; } // -1 is for the initial $ of the first word
            else { pos_for_sai += dict_string.size() - this->params.w; }
            sai_file.write((char*) &pos_for_sai, IBYTES);
        }
        rw_mmap.unmap();
        truncate_file(out_file_name, this->parse_size * sizeof(size_type));
    }
    
    vcfbwt::DiskWrites::update(last_file.tellp());
    last_file.close();
    
    vcfbwt::DiskWrites::update(sai_file.tellp());
    sai_file.close();
    
    // Print dicitionary on disk
    spdlog::info("Main parser: writing dictionary on disk NOT COMPRESSED");
    std::string dict_file_name = out_file_prefix + EXT::DICT;
    std::ofstream dict(dict_file_name);
    
    for (size_type i = 0; i < this->dictionary.size(); i++)
    {
        dict.write(this->dictionary.sorted_entry_at(i).c_str(), this->dictionary.sorted_entry_at(i).size());
        dict.put(ENDOFWORD);
    }
    
    dict.put(ENDOFDICT);
    
    vcfbwt::DiskWrites::update(dict.tellp()); // Disk Stats
    dict.close();
    
    if (this->params.compress_dictionary)
    {
        spdlog::info("Main parser: writing dictionary on disk COMPRESSED");
        std::ofstream dicz(out_file_prefix + EXT::DICT_COMPRESSED);
        std::ofstream lengths(out_file_prefix + EXT::DICT_COMPRESSED_LENGTHS);

        for (size_type i = 0; i < this->dictionary.size(); i++)
        {
            std::size_t shift = 1; // skip dollar on first phrase
            if (i != 0) { shift = this->w; }
            dicz.write(this->dictionary.sorted_entry_at(i).c_str() + shift,
                       this->dictionary.sorted_entry_at(i).size() - shift);
            int32_t len = this->dictionary.sorted_entry_at(i).size() - shift;
            lengths.write((char*) &len, sizeof(int32_t));
        }
        
        vcfbwt::DiskWrites::update(dicz.tellp()); // Disk Stats
        dicz.close();
        
        vcfbwt::DiskWrites::update(lengths.tellp()); // Disk Stats
        lengths.close();
    }
    
    // Outoput Occurrencies
    if(this->params.compute_occurrences)
    {
        spdlog::info("Main parser: writing occurrences to file");
        std::string occ_file_name = out_file_prefix + EXT::OCC;
        std::ofstream occ(occ_file_name, std::ios::out | std::ios::binary);
        occ.write((char*)&occurrences[0], occurrences.size() * sizeof(size_type));
        
        vcfbwt::DiskWrites::update(occ.tellp()); // Disk Stats
        occ.close();
    }

    // Merging the length components
    if (params.report_lengths or params.compute_lifting)
    {
        std::ofstream merged(out_len_name);

        // Reference
        for(int i = 0; i < this->references.size(); i++)
        {
            const size_t length = this->references[i].size();
            const std::string contig_name = this->references_name[i];
            merged << contig_name << " " << length << std::endl;
        }

        vcfbwt::DiskWrites::update(merged.tellp()); // Disk Stats
        merged.close();
    }

    // Merging the lifting components
    if (params.compute_lifting)
    {
        //Opens output lift file name as merged 
        std::ofstream merged(out_lift_name, std::ios_base::binary);

        // if compute_lifting => report_lengths
        std::vector<size_t> onset(1,0);
        std::vector<size_t> lengths;
        size_t u = 0;            
        size_t w = this->params.w;            
        std::vector<std::string> names;

        // Reading the lengths
        std::string tmp_name;
        std::size_t tmp_length;
        // while (not in_lidx.eof()) 
        std::ifstream in_lidx(out_len_name);
        while (in_lidx >> tmp_name >> tmp_length )
        { 
            if (tmp_name != "")
            {
                u += tmp_length;
                names.push_back(tmp_name);
                lengths.push_back(tmp_length);
                onset.push_back(u);
            }
        }
        ++u;
        in_lidx.close();
        // Build the seqidx structure
        sdsl::sd_vector_builder builder(u, onset.size());
        for (auto idx : onset)
            builder.set(idx);

        sdsl::sd_vector<> starts(builder);
        // Writing the seqidx on disk
        merged.write((char *)&u, sizeof(u));
        merged.write((char *)&w, sizeof(w));

        starts.serialize(merged);
        sdsl::serialize(names.size(), merged);
        for(size_t i = 0; i < names.size(); ++i)
        {
            sdsl::serialize(names[i].size(), merged);
            merged.write((char *)names[i].data(), names[i].size());
        }

        const size_t reference_count = this->references.size();
        merged.write((char *)&reference_count, sizeof(reference_count));

        // Build the empty liftings for the references
        size_t clen = 0;
        for(size_t i = 0; i < this->references.size(); ++i)
        {
            const size_t len = references[i].size();

            sdsl::bit_vector ibv(len);
            sdsl::bit_vector dbv(len);
            sdsl::bit_vector sbv(len);

            lift::Lift lift(ibv, dbv, sbv); 

            sdsl::serialize(clen, merged);
            lift.serialize(merged);  

            clen += len;       
        }

        vcfbwt::DiskWrites::update(merged.tellp()); // Disk Stats
        merged.close();
    }
}

//------------------------------------------------------------------------------

void
vcfbwt::pfp::ParserText::init(const Params& params, const std::string& prefix)
{
    this->params = params;
    this->w = this->params.w;
    this->p = this->params.p;
    this->parse_size = 0;
    this->out_file_prefix = prefix;
    this->out_file_name = prefix + EXT::PARSE;
    this->out_file.open(out_file_name);
}

void
vcfbwt::pfp::ParserText::operator()()
{
    // Open input file with kseq
    gzFile fp;
    fp = gzopen(this->in_file_path.c_str(), "r");
    if (fp == 0)
    {
        spdlog::error("Failed to open input file {}", in_file_path);
        exit(EXIT_FAILURE);
    }
    
    std::string phrase;
    spdlog::info("Parsing {}", in_file_path);
    
    // Karp Robin Hash Function for sliding window
    KarpRabinHash kr_hash(this->params.w);
    
    // First sequence start with one dollar
    phrase.append(1, DOLLAR);
    
    char c;
    while(gzread(fp, &c, 1) > 0)
    {
        phrase.push_back(c);
        if (phrase.size() == params.w) { kr_hash.initialize(phrase); }
        else if (phrase.size() > params.w) { kr_hash.update(phrase[phrase.size() - params.w - 1], phrase[phrase.size() - 1]); }
        
        if ((phrase.size() > this->params.w) and ((kr_hash.get_hash() % this->params.p) == 0))
        {
            std::string_view ts(&(phrase[phrase.size() - params.w]), params.w);
            hash_type ts_hash = KarpRabinHash::string_hash(ts);
            
            hash_type hash = this->dictionary.check_and_add(phrase);
            
            out_file.write((char*) (&hash), sizeof(hash_type)); this->parse_size += 1;
            
            phrase.erase(phrase.begin(), phrase.end() - this->params.w); // Keep the last w chars
            
            kr_hash.reset(); kr_hash.initialize(phrase);
        }
    }
    
    // Last phrase
    if (phrase.size() > this->params.w)
    {
        // Append w dollar at the end
        phrase.append(this->params.w, DOLLAR);
        
        hash_type hash = this->dictionary.check_and_add(phrase);
        
        out_file.write((char*) (&hash), sizeof(hash_type)); this->parse_size += 1;
    }
    else { spdlog::error("A sequence doesn't have w DOLLAR at the end!"); std::exit(EXIT_FAILURE); }
    
    gzclose(fp);
}


void
vcfbwt::pfp::ParserText::close()
{
    if (closed) return; closed = true;
    
    vcfbwt::DiskWrites::update(out_file.tellp()); // Disk Stats
    this->out_file.close();
    
    // Occurrences
    std::vector<size_type> occurrences(this->dictionary.size(), 0);
    
    spdlog::info("Main parser: Replacing hash values with ranks, writing .last and .sai");
  
    std::string last_file_name = out_file_prefix + EXT::LAST;
    std::ofstream last_file(last_file_name);
    
    std::string sai_file_name = out_file_prefix + EXT::SAI;
    std::ofstream sai_file(sai_file_name);
    
    std::size_t pos_for_sai = 0;
    
    // mmap file and substitute
    if (this->parse_size != 0)
    {
        std::error_code error;
        mio::mmap_sink rw_mmap = mio::make_mmap_sink(out_file_name, 0, mio::map_entire_file, error);
        if (error) { spdlog::error(error.message()); std::exit(EXIT_FAILURE); }
        
        for (size_type i = 0; i < (rw_mmap.size() / sizeof(hash_type)); i++)
        {
            hash_type hash;
            std::memcpy(&hash, rw_mmap.data() + (i * sizeof(hash_type)), sizeof(hash_type));
            size_type rank = this->dictionary.hash_to_rank(hash);
            std::memcpy(rw_mmap.data() + (i * sizeof(size_type)), &rank, sizeof(size_type));
            occurrences[rank - 1] += 1;
            
            const std::string& dict_string = this->dictionary.sorted_entry_at(rank - 1);
            last_file.put(dict_string[(dict_string.size() - this->params.w) - 1]);
    
            if (pos_for_sai == 0) { pos_for_sai = dict_string.size() - 1; } // -1 is for the initial $ of the first word
            else { pos_for_sai += dict_string.size() - this->params.w; }
            sai_file.write((char*) &pos_for_sai, IBYTES);
        }
        rw_mmap.unmap();
        truncate_file(out_file_name, this->parse_size * sizeof(size_type));
    }
    
    vcfbwt::DiskWrites::update(last_file.tellp());
    last_file.close();
    
    vcfbwt::DiskWrites::update(sai_file.tellp());
    sai_file.close();
    
    // Print dicitionary on disk
    spdlog::info("Main parser: writing dictionary on disk NOT COMPRESSED");
    std::string dict_file_name = out_file_prefix + EXT::DICT;
    std::ofstream dict(dict_file_name);
    
    for (size_type i = 0; i < this->dictionary.size(); i++)
    {
        dict.write(this->dictionary.sorted_entry_at(i).c_str(), this->dictionary.sorted_entry_at(i).size());
        dict.put(ENDOFWORD);
    }
    
    dict.put(ENDOFDICT);
    
    vcfbwt::DiskWrites::update(dict.tellp()); // Disk Stats
    dict.close();
    
    if (this->params.compress_dictionary)
    {
        spdlog::info("Main parser: writing dictionary on disk COMPRESSED");
        std::ofstream dicz(out_file_prefix + EXT::DICT_COMPRESSED);
        std::ofstream lengths(out_file_prefix + EXT::DICT_COMPRESSED_LENGTHS);
        
        for (size_type i = 0; i < this->dictionary.size(); i++)
        {
            std::size_t shift = 1; // skip dollar on first phrase
            if (i != 0) { shift = this->w; }
            dicz.write(this->dictionary.sorted_entry_at(i).c_str() + shift,
                       this->dictionary.sorted_entry_at(i).size() - shift);
            int32_t len = this->dictionary.sorted_entry_at(i).size() - shift;
            lengths.write((char*) &len, sizeof(int32_t));
        }
        
        vcfbwt::DiskWrites::update(dicz.tellp()); // Disk Stats
        dicz.close();
        
        vcfbwt::DiskWrites::update(lengths.tellp()); // Disk Stats
        lengths.close();
    }
    
    // Output Occurrencies
    if(this->params.compute_occurrences)
    {
        spdlog::info("Main parser: writing occurrences to file");
        std::string occ_file_name = out_file_prefix + EXT::OCC;
        std::ofstream occ(occ_file_name, std::ios::out | std::ios::binary);
        occ.write((char*)&occurrences[0], occurrences.size() * sizeof(size_type));
        
        vcfbwt::DiskWrites::update(occ.tellp()); // Disk Stats
        occ.close();
    }
    
}

//------------------------------------------------------------------------------

void
vcfbwt::pfp::ParserUtils::read_parse(std::string parse_file_name, std::vector<size_type>& parse)
{
    std::ifstream parse_file(parse_file_name, std::ios::binary);
    if (not parse_file.is_open()) { spdlog::error("Error opening file: {}", parse_file_name); std::exit(EXIT_FAILURE); }
    
    while (not parse_file.eof()) { size_type i; parse_file.read((char*) &i, sizeof(size_type)); parse.push_back(i); }
    
    // remove last entry, apparently it reads the last twice, strano
    parse.pop_back();
}

void
vcfbwt::pfp::ParserUtils::read_dictionary(std::string dic_file_name, std::vector<std::string>& dictionary_vector)
{
    std::ifstream dic_file(dic_file_name);
    if (not dic_file.is_open()) { spdlog::error("Error opening file: {}", dic_file_name); std::exit(EXIT_FAILURE); }
    
    char c = '\5';
    while (c != ENDOFDICT)
    {
        std::string phrase = "";
        while ((c = dic_file.get()) != ENDOFWORD and c != ENDOFDICT)
        {
            phrase.push_back(c);
        }
        if (phrase.size() != 0) { dictionary_vector.push_back(phrase); }
    }
}

void vcfbwt::pfp::ParserUtils::read_compressed_dictionary(std::string dic_file_name, std::string len_file_name, std::vector<std::string> &dictionary_vector)
{
    std::ifstream dic_file(dic_file_name);
    if (not dic_file.is_open()) { spdlog::error("Error opening file: {}", dic_file_name); std::exit(EXIT_FAILURE); }
    std::ifstream len_file(len_file_name);
    if (not len_file.is_open()) { spdlog::error("Error opening file: {}", len_file_name); std::exit(EXIT_FAILURE); }

    std::streampos fsize = len_file.tellg();
    len_file.seekg(0, std::ios::end);
    fsize = len_file.tellg() - fsize;
    len_file.seekg(0, std::ios::beg);

    size_t n_phrases = fsize/4;
    if ((fsize % 4) != 0) { spdlog::error("Invalid file: {}", len_file_name); std::exit(EXIT_FAILURE); }

    dictionary_vector.reserve(n_phrases);

    for( size_t i = 0; i < n_phrases; ++i)
    {
        int32_t len = 0;
        std::string phrase = "";
        len_file.read((char*)&len, sizeof(int32_t));
        phrase.resize(len);
        dic_file.read(phrase.data(), len);
        dictionary_vector.push_back(phrase);
    }
}

void
vcfbwt::pfp::ParserUtils::merge(std::string left_prefix, std::string right_prefix, std::string out_prefix, const Params& params)
{
    // Merged Dictionary and Parse
    Dictionary dictionary;
    
    // Read Left and Right Dictionary
    spdlog::info("Loading dictionaries from disk");
    std::vector<std::string> left_dictionary;   read_dictionary(left_prefix + EXT::DICT, left_dictionary);
    std::vector<std::string> right_dictionary; read_dictionary(right_prefix + EXT::DICT, left_dictionary);
    
    // Read Left parse and substitute ranks with hash again storing in the merged dictionary
    spdlog::info("Adjusting ranks");
    std::error_code error;
    mio::mmap_sink rw_mmap = mio::make_mmap_sink(left_prefix + EXT::PARSE, 0, mio::map_entire_file, error);
    if (error) { spdlog::error(error.message()); std::exit(EXIT_FAILURE); }
    
    for (size_type i = 0; i < (rw_mmap.size() / sizeof(hash_type)); i++)
    {
        hash_type rank;
        std::memcpy(&rank, rw_mmap.data() + (i * sizeof(hash_type)), sizeof(hash_type));
        
        std::string& phrase = left_dictionary[rank];
        hash_type hash = 0;
        if (dictionary.contains(phrase))    { hash = dictionary.get(phrase); }
        else                                { hash = dictionary.add(phrase); }
        
        std::memcpy(rw_mmap.data() + (i * sizeof(hash_type)), &hash, sizeof(hash_type));
    }
    rw_mmap.unmap(); left_dictionary.resize(0);
    
    // Read right parse, changing first phrase;
    rw_mmap = mio::make_mmap_sink(left_prefix + EXT::PARSE, 0, mio::map_entire_file, error);
    if (error) { spdlog::error(error.message()); std::exit(EXIT_FAILURE); }
    
    for (size_type i = 0; i < (rw_mmap.size() / sizeof(hash_type)); i++)
    {
        hash_type rank;
        std::memcpy(&rank, rw_mmap.data() + (i * sizeof(hash_type)), sizeof(hash_type));
        
        std::string& phrase = right_dictionary[rank];
        
        if (rank == 1) { phrase[0] = DOLLAR_PRIME; phrase.insert( 0, params.w - 1, DOLLAR_PRIME); }
        
        hash_type hash = 0;
        if (dictionary.contains(phrase))    { hash = dictionary.get(phrase); }
        else                                { hash = dictionary.add(phrase); }
        
        std::memcpy(rw_mmap.data() + (i * sizeof(hash_type)), &hash, sizeof(hash_type));
    }
    rw_mmap.unmap(); right_dictionary.resize(0);
    
    // Sort the Dictionary
    spdlog::info("Sorting merged dictionary");
    std::vector<std::pair<std::string, hash_type>> entries;
    for (auto& entry : dictionary.hash_string_map)
    {
        entries.push_back(std::make_pair(entry.second.phrase, entry.first));
    }
    std::sort(entries.begin(), entries.end());
    
    // Insert in hashmap
    std::unordered_map<hash_type, hash_type> hash_to_ranks;
    for (hash_type i = 0; i < entries.size(); i++)
    {
        hash_to_ranks.insert(std::make_pair(entries[i].second, i + 1)); // 1 based
    }
    
    // Merge parsings
    spdlog::info("Creating output parse");
    std::ofstream out_parse(out_prefix + EXT::PARSE, std::ios::binary);
    
    // Left
    mio::mmap_source rl_mmap; rl_mmap.map(left_prefix + EXT::PARSE, error);
    if (error) { spdlog::error(error.message()); std::exit(EXIT_FAILURE); }
    for (size_type i = 0; i < (rl_mmap.size() / sizeof(hash_type)); i++)
    {
        hash_type hash;
        std::memcpy(&hash, rl_mmap.data() + (i * sizeof(hash_type)), sizeof(hash_type));
        
        auto it = hash_to_ranks.find(hash);
        if (it != hash_to_ranks.end()) { out_parse.write((char*) &(it->second), sizeof(hash_type)); }
        else { spdlog::error("Something wrong happend"); std::exit(EXIT_FAILURE); }
    }
    rl_mmap.unmap();
    
    // Right
    mio::mmap_source rr_mmap; rr_mmap.map(left_prefix + EXT::PARSE, error);
    if (error) { spdlog::error(error.message()); std::exit(EXIT_FAILURE); }
    for (size_type i = 0; i < (rr_mmap.size() / sizeof(hash_type)); i++)
    {
        hash_type hash;
        std::memcpy(&hash, rr_mmap.data() + (i * sizeof(hash_type)), sizeof(hash_type));
        
        auto it = hash_to_ranks.find(hash);
        if (it != hash_to_ranks.end()) { out_parse.write((char*) &(it->second), sizeof(hash_type)); }
        else { spdlog::error("Something wrong happend"); std::exit(EXIT_FAILURE); }
    }
    rr_mmap.unmap();
    out_parse.close();
    
    // Print dicitionary on disk
    spdlog::info("Writing dictionary to disk NOT COMPRESSED");
    std::string dict_file_name = out_prefix + EXT::DICT;
    std::ofstream dict(dict_file_name);
    
    for (auto& entry : entries)
    {
        dict.write(entry.first.c_str(), entry.first.size());
        dict.put(ENDOFWORD);
    }
    
    dict.put(ENDOFDICT);
    dict.close();
}
