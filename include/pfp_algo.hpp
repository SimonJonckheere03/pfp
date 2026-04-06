//
//  pfp_algo.hpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//

#ifndef pfp_algo_hpp
#define pfp_algo_hpp

#include <vector>
#include <unordered_map>
#include <set>
#include <iostream>
#include <fstream>
#include <cstddef>
#include <vcf.hpp>
#include <utils.hpp>
#include <internals.hpp>
#include <leviosam.hpp>

namespace vcfbwt
{
namespace pfp
{

//------------------------------------------------------------------------------

enum SPECIAL_TYPES
{
    ENDOFDICT = 0,
    ENDOFWORD = 1,
    DOLLAR = 2,
    DOLLAR_SEQUENCE = 4,
    DOLLAR_PRIME = 5
};

static const unsigned char acgt_only_table[256] = {
      0,  1,  2,  3,  4,  5,'N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','A','N','C','N','N','N','G','N','N','N','N','N','N','N','N',
    'N','N','N','N','T','N','N','N','N','N','N','N','N','N','N','N',
    'N','A','N','C','N','N','N','G','N','N','N','N','N','N','N','N',
    'N','N','N','N','T','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    };

//------------------------------------------------------------------------------

class Dictionary
{
public:

    vcfbwt::size_type insertions_safe_guard = 1000;

    std::mutex dictionary_mutex;

    bool sorted = false;
    
    struct DictionaryEntry
    {
        std::string phrase;
    
        DictionaryEntry() = default;
        DictionaryEntry(const std::string& s) : phrase(s) {}
    };
    
    std::unordered_map<hash_type, DictionaryEntry> hash_string_map;
    std::vector<std::pair<std::reference_wrapper<std::string>, hash_type>> sorted_phrases;
    std::unordered_map<hash_type, size_type> hash_to_ranks;
    
    void sort();
    
    Dictionary() = default;
    
    hash_type add(const std::string& phrase);
    hash_type check_and_add(const std::string& phrase);
    hash_type get(const std::string& phrase) const;
    bool contains(const std::string& phrase);

    size_type hash_to_rank(hash_type hash);
    
    size_type size() const { return hash_string_map.size(); }
    
    const std::string& sorted_entry_at(std::size_t i);
    
    static void merge(Dictionary& destination, const Dictionary& source);
    
    friend class ParserVCF;
};

//------------------------------------------------------------------------------

struct Params
{
    hash_type p = 100;
    hash_type w =  10;
    std::size_t seed_length = 100;
    bool compress_dictionary = false;
    bool use_acceleration = false;
    bool print_out_statistics_csv = false;
    bool compute_occurrences = true;
    bool auPair = false;
    bool compute_lifting = false;
    bool report_lengths = false;
    bool acgt_only = false;
    std::string ignore_ts_file;
    std::string haplotype_name_prefix;
};

struct Statistics
{
    std::size_t parse_length = 0;
    std::size_t total_dictionary_length = 0;
    std::size_t num_of_phrases_dictionary = 0;
};

class ReferenceParse
{

    const std::size_t ref_length;
    const std::string ref_id;
public :
    Dictionary& dictionary;
    Dictionary& dictionary_rev;
    std::vector<hash_type> parse;
    std::string segment_text;

    std::vector<long long int> trigger_strings_position; // position of first char of each trigger string
    // std::vector<size_type> trigger_strings_position; // position of first char of each trigger string
    std::set<hash_type> to_ignore_ts_hash;
    
    const Params& params;

    const std::string& id() const { return this->ref_id; }
    const std::size_t length() const { return this->ref_length; }
    
    void init(const std::string& reference, bool first = true);
    
    ReferenceParse(const std::string& reference, const std::string& id, Dictionary& dict, Dictionary& dict_rev, const Params& pms, bool first = true) : dictionary(dict), dictionary_rev(dict_rev), params(pms), ref_length(reference.size()), ref_id(id) { this->init(reference, first); }
    const hash_type& operator[](size_type i) const { return this->parse[i]; }
};

// Create a parse on disk
class ParserVCF
{

private:
    
    std::ofstream out_file;
    std::string out_file_prefix;
    std::string out_file_name;
    std::string tmp_out_file_name;

    // Reverse file variables
    std::string out_file_name_rev;
    size_type parse_size_rev = 0;

    std::ofstream out_lift;
    std::string out_lift_prefix;
    std::string out_lift_name;
    std::string tmp_out_lift_name;

    std::ofstream out_len;
    std::string out_len_prefix;
    std::string out_len_name;
    std::string tmp_out_len_name;
    
    struct ReverseSegment
    {
        std::string text;
        const std::set<hash_type>* ignored_trigger_strings = nullptr;
    };

    std::vector<std::string> samples_processed;
    std::vector<std::pair<std::string, std::string>> contigs_processed;
    std::vector<ReverseSegment> reverse_segments;
    
    Params params;
    Statistics statistics;
    
    std::vector<ReferenceParse>* references_parse = nullptr;
    Dictionary* dictionary = nullptr;
    Dictionary* dictionary_rev = nullptr; // NEW: The missing reverse dictionary pointer

    // Shorthands
    hash_type w, p;
    size_type parse_size, tags;
    
    bool closed = false;
    
    std::vector<std::reference_wrapper<ParserVCF>> registered_workers;
    
    std::size_t working_genotype = 0;
    
public:
    
    enum tags
    {
        MAIN = 1,
        WORKER = 2,
        COMPRESSED = 4,
        UNCOMPRESSED = 8
    };
    
    void init(const Params& params, const std::string& out_prefix, std::vector<ReferenceParse>& rp, Dictionary& dict, Dictionary& dict_rev, std::size_t t = MAIN | UNCOMPRESSED);
    
    ParserVCF(const Params& params, const std::string& out_prefix, std::vector<ReferenceParse>& rp, Dictionary& dict, Dictionary& dict_rev, std::size_t t = MAIN | UNCOMPRESSED)
    {
        if (out_prefix.empty()) { this->init(params, "out", rp, dict, dict_rev, t); }
        else { this->init(params, out_prefix, rp, dict, dict_rev, t); }
    }
    
    ParserVCF() = default;
    
    ~ParserVCF()
    {
        close();
        size_type total_length = 0;
        for (auto& entry : dictionary->hash_string_map) { total_length += entry.second.phrase.size(); }
    
        // Fill out statistics
        this->statistics.parse_length = this->parse_size;
        this->statistics.num_of_phrases_dictionary = this->dictionary->sorted_phrases.size();
        this->statistics.total_dictionary_length = total_length;
        std::string name;
        if (tags & WORKER) { name = "Worker"; } else { name = "Main"; }
        spdlog::info("{} -\tParse size: {}\tDic Size: {} Dic Total Length: {}", name ,parse_size, dictionary->size(), total_length);
        
        if (params.print_out_statistics_csv and (tags & MAIN))
        {
            std::ofstream csv(out_file_prefix + ".csv");
            csv << "w,p,f,parse_lenght,dict_phrases,dict_tot_length\n";
            csv << params.w << ",";
            csv << params.p << ",";
            csv << statistics.parse_length << ",";
            csv << statistics.num_of_phrases_dictionary << ",";
            csv << statistics.total_dictionary_length << ",";
            csv << "\n";
            csv.close();
        }
    }
    
    const std::string& get_file_name() const { return this->out_file_name; }
    const Statistics& get_statistics() const { return this->statistics; }
    void register_worker(ParserVCF& parser) { this->registered_workers.push_back(std::ref(parser)); }
    void set_working_genotype(std::size_t genotype) { this->working_genotype = genotype; }
    
    void operator()(const Sample& sample);
    void close();
};

//------------------------------------------------------------------------------

class ParserFasta
{

private:
    
    std::ofstream out_file;
    std::string out_file_prefix;
    std::string out_file_name;
    std::string in_file_path;
    std::vector<std::string> sequences_processed;

    std::ofstream out_lift;
    std::string out_lift_prefix;
    std::string out_lift_name;
    std::string tmp_out_lift_name;

    std::ofstream out_len;
    std::string out_len_prefix;
    std::string out_len_name;
    std::string tmp_out_len_name;

    std::string reference; // TODO: Replace it with as many parses as references.
    std::vector<std::string> references;
    std::vector<std::string> references_name;
    std::unordered_map<std::string, std::size_t> references_id;
    
    std::vector<std::vector<Variation>> variations;
    std::vector<Contig> contigs;
    std::vector<Sample> samples;
    std::unordered_map<std::string, std::size_t> samples_id;
    
    std::size_t max_samples = 0;
    std::set<std::string> input_samples;
    std::vector<std::size_t> populated_samples;

    std::vector<std::size_t> ref_sum_lengths;
    
    Params params;
    Statistics statistics;
    
    Dictionary dictionary;
    
    hash_type w, p;
    size_type parse_size, tags;
    
    bool closed = false;
    
public:
    
    void init(const Params& params, const std::string& prefix);

    void init_ref(const std::string& ref_path, const size_t w, bool last);
    
    ParserFasta(const Params& params, const std::string& file_path, const std::string& out_prefix)
    {
        this->in_file_path = file_path;
        if (out_prefix.empty()) { this->init(params, "out"); }
        else { this->init(params, out_prefix); }
    }
    
    ParserFasta() = default;
    
    ~ParserFasta()
    {
        close();
        size_type total_length = 0;
        for (auto& entry : dictionary.hash_string_map) { total_length += entry.second.phrase.size(); }
        
        // Fill out statistics
        this->statistics.parse_length = this->parse_size;
        this->statistics.num_of_phrases_dictionary = this->dictionary.sorted_phrases.size();
        this->statistics.total_dictionary_length = total_length;
        std::string name = "Parser Fasta";
        spdlog::info("{} -\tParse size: {}\tDic Size: {} Dic Total Length: {}", name ,parse_size, dictionary.size(), total_length);
        
        if (params.print_out_statistics_csv)
        {
            std::ofstream csv(out_file_prefix + ".csv");
            csv << "w,p,f,parse_lenght,dict_phrases,dict_tot_length\n";
            csv << params.w << ",";
            csv << params.p << ",";
            csv << statistics.parse_length << ",";
            csv << statistics.num_of_phrases_dictionary << ",";
            csv << statistics.total_dictionary_length << ",";
            csv << "\n";
            csv.close();
        }
    }
    
    const std::string& get_file_name() const { return this->out_file_name; }
    const Statistics& get_statistics() const { return this->statistics; }
    
    void operator()();
    void close();
};

//------------------------------------------------------------------------------

class ParserText
{

private:
    
    std::ofstream out_file;
    std::string out_file_prefix;
    std::string out_file_name;
    std::string in_file_path;
    std::vector<std::string> sequences_processed;
    
    Params params;
    Statistics statistics;
    
    Dictionary dictionary;
    
    hash_type w, p;
    size_type parse_size, tags;
    
    bool closed = false;

public:
    
    void init(const Params& params, const std::string& prefix);
    
    ParserText(const Params& params, const std::string& file_path, const std::string& out_prefix)
    {
        this->in_file_path = file_path;
        if (out_prefix.empty()) { this->init(params, "out"); }
        else { this->init(params, out_prefix); }
    }
    
    ParserText() = default;
    
    ~ParserText()
    {
        close();
        size_type total_length = 0;
        for (auto& entry : dictionary.hash_string_map) { total_length += entry.second.phrase.size(); }
        
        // Fill out statistics
        this->statistics.parse_length = this->parse_size;
        this->statistics.num_of_phrases_dictionary = this->dictionary.sorted_phrases.size();
        this->statistics.total_dictionary_length = total_length;
        std::string name = "Parser Text";
        spdlog::info("{} -\tParse size: {}\tDic Size: {} Dic Total Length: {}", name ,parse_size, dictionary.size(), total_length);
        
        if (params.print_out_statistics_csv)
        {
            std::ofstream csv(out_file_prefix + ".csv");
            csv << "w,p,f,parse_lenght,dict_phrases,dict_tot_length\n";
            csv << params.w << ",";
            csv << params.p << ",";
            csv << statistics.parse_length << ",";
            csv << statistics.num_of_phrases_dictionary << ",";
            csv << statistics.total_dictionary_length << ",";
            csv << "\n";
            csv.close();
        }
    }
    
    const std::string& get_file_name() const { return this->out_file_name; }
    const Statistics& get_statistics() const { return this->statistics; }
    
    void operator()();
    void close();
};

//------------------------------------------------------------------------------

class ParserUtils
{
public:
    
    static void read_parse(std::string parse_file_name, std::vector<size_type>& parse);
    static void read_dictionary(std::string dic_file_name, std::vector<std::string>& dictionary_vector);
    static void read_compressed_dictionary(std::string dic_file_name, std::string len_file_name, std::vector<std::string> &dictionary_vector);
    static void merge(std::string left_prefix, std::string right_prefix, std::string out_prefix, const Params &params);
    static void parse_fasta(std::string fasta_file_name, std::string out_prefix, const Params& params);
    
    static std::vector<std::size_t> compute_occurrences(std::vector<std::string>& dictionary_vector, std::vector<size_type>& parse);
};

} // end namespace pfp
} // end namespace vcfbwt

#endif //pfp_algo_hpp
