#include "common.hpp"
#include "arg_parse.hpp"
#include "gff_reader.hpp"

#include <omp.h>

#include <unordered_map>

#include <regex>
#include <string>
#include <vector>
#include <tuple>

#include <stdio.h>
#include <sys/wait.h>

struct AnnotateWithMSACLIParams
{
    std::string output_path = "";
    std::string mmseqs2_bin = "mmseqs";
    unsigned threads = omp_get_max_threads();
//    algorithm_t strategy = MLE;
//    bool comp_bls = true;

    std::string reference_genome_name;
    std::string reference_genome_path;

    std::vector<std::tuple<std::string, std::string> > aligning_genomes;
};

// TODO: make this function prettier (e.g., no local struct def)
void mmseqs_fasta_to_maf(const std::string & src, const std::string & dest, const AnnotateWithMSACLIParams & params,
                         const std::unordered_map<std::string, uint32_t> & lookup_genome_ids, const uint8_t phase)
{
    struct maf_object
    {
        std::string chr;
        uint64_t begin;
        uint64_t end;
        char strand;

        std::string seq;

        std::vector<std::tuple<std::string, std::string> > aln;

        void print(FILE *fp, const uint8_t phase)
        {
            // create format string (depends on max sequence name length)s
            char format_str[50];
            uint16_t max_seq_name_len = chr.size();
            for (auto const & a : aln)
            {
                if (std::get<0>(a).size() > max_seq_name_len)
                {
                    max_seq_name_len = std::get<0>(a).size();
                }
            }
            sprintf(format_str, "s %%-%ds %%10ld %%10ld %%c %%ld %%s\n", max_seq_name_len);

            fprintf(fp, "a score=NAN\n");
            // seq name, begin pos (0-indexed), length, strand, total chrom len, seq
            // NOTE: an entry on the reverse strand already has the sequence reverse-complemented, i.e., it would actually now be on the forward strand
            // but we need to remember what strand the original sequence was on to easily map them to the CDS in the GFF file
            seq.erase(0, phase); // remove the first 0-2 bases

            fprintf(fp, format_str, chr.c_str(), begin - 1, end - (begin - 1), strand, 0, seq.c_str());
            for (auto & a : aln)
            {
                // we don't care about locations or strandness
                std::get<1>(a).erase(0, phase); // remove the first 0-2 bases
                fprintf(fp, format_str, std::get<0>(a).c_str(), 0ul, 0, '+', 0, std::get<1>(a).c_str());
            }
            fprintf(fp, "\n");
        }
    };

    FILE *f_in = fopen(src.c_str(), "r");
    if (f_in == NULL)
    {
        printf(OUT_ERROR "Cannot open file for reading: %s\n" OUT_RESET, src.c_str());
        exit(-1);
    }

    FILE *f_out = fopen(dest.c_str(), "w");
    if (f_out == NULL)
    {
        printf(OUT_ERROR "Cannot open file for reading: %s\n" OUT_RESET, dest.c_str());
        exit(-1);
    }

    char * line = NULL;
    size_t len = 0;
    ssize_t read;

    // TODO: get rid of regex and use sscanf
    const std::regex pieces_regex("(.+):([0-9]+)-([0-9]+)#([\\+-])");
    std::smatch pieces_match;

    maf_object m;

    std::string id;

    while ((read = getline(&line, &len, f_in)) != -1)
    {
        if (read == 0) // empty line
            continue;

        if (line[read - 1] == '\n')
            line[read - 1] = 0; // remove newline

        // mmseqs2 outputs a 0x00 char before the beginning of each alignment (i.e., before the identifier)
        char *line2 = line;
        if (line2[0] == 0)
            line2 = line + 1;

        if (line2[0] == '>')
        {
            id = line2;
            id.erase(0, 1); // remove ">"

            if (std::regex_match(id, pieces_match, pieces_regex))
            {
                if (m.aln.size() > 0)
                    m.print(f_out, phase);

//                printf("ref id\n");

                m.chr = params.reference_genome_name + "." + pieces_match[1].str();
                m.begin = std::stoi(pieces_match[2].str());
                m.end = std::stoi(pieces_match[3].str());
                m.strand = pieces_match[4].str()[0];
                m.seq = "";
                m.aln.clear();

                id = ""; // NOTE: this is impotant. last line in file might only by 0x00, hence it will try to add a new sequence to the alignment

                //printf("XXX: %s\t%ld\t%ld\t%c\n", m.chr.c_str(), m.begin, m.end, m.strand);
            }
            else
            {
//                printf("other id\n");

                const size_t space = id.find_first_of(' ');
                if (space != std::string::npos)
                    id = id.substr(0, space); // remove everything after the first space
            }
        }
        else
        {
            if (m.chr != "" && m.seq == "")
            {
//                printf("ref seq\n");
                m.seq = line2;
            }
            else
            {
                if (id != "")
                {
//                    printf("other seq: x%sx (%ld)\n", line2, strlen(line2));
                    auto genome_id = lookup_genome_ids.find(id);
                    if (genome_id != lookup_genome_ids.end())
                        m.aln.emplace_back(std::get<0>(params.aligning_genomes[genome_id->second]) + ".UNK", line2);
                    else
                        printf(OUT_ERROR "Could not match sequence id to genome (this should never happen!) in file %s: %s\n" OUT_RESET, src.c_str(), id.c_str());
                    id = ""; // at the end of the file we might encounter a newline (i.e, it does does not start with '>') and might think a sequene is following
                }
            }
        }
    }

    if (m.aln.size() > 0)
        m.print(f_out, phase);

    free(line);

    fclose(f_in);
    fclose(f_out);
}

void load_fasta_file(const std::string & reference_path, std::unordered_map<std::string, std::string> & genome)
{
    FILE *file = fopen(reference_path.c_str(), "r");
    if (file == NULL)
    {
        printf(OUT_ERROR "Cannot open genomic fasta file %s\n" OUT_RESET, reference_path.c_str());
        exit(-1);
    }

    std::string id;
    std::string seq;
    seq.reserve(250000000); // reserve 250MB of seq (length of hg38.chr1)

    char line[BUFSIZ];
    while (fgets(line, sizeof line, file) != NULL)
    {
        if (line[0] == '>')
        {
            // write previous sequence out
            if (seq != "")
            {
                genome.emplace(id, seq);
                // printf("'%s'\t%ld\t%ld\n", id.c_str(), genome.find(id)->second.size(), genome.find(id)->second.capacity());
            }

            seq = "";

            // extract new identifier
            id = line;
            id = id.substr(1); // remove first character '>'
            while (id[0] == ' ')
                id = id.substr(1); // remove space(s) after '>'

            const size_t next_space = id.find(' ');
            if (next_space != std::string::npos)
                id = id.substr(0, next_space); // remove everything after the next space

            if (id.back() == '\n')
                id.erase(id.size() - 1);
        }
        else
        {
            char *newline = strchr(line, '\n');
            *newline = 0;
            seq += line;
        }
    }

    if (seq != "")
        genome.emplace(id, seq);

    fclose(file);
}

void load_genome_file(const std::string & genome_file, std::vector<std::tuple<std::string, std::string> > & genome_file_paths,
                      std::string & reference_genome_name, std::string & reference_genome_path)
{
    FILE *file = fopen(genome_file.c_str(), "r");
    if (file == NULL)
    {
        printf(OUT_ERROR "Cannot open genome-file %s\n" OUT_RESET, genome_file.c_str());
        exit(-1);
    }

    bool first_entry = true;
    char line[BUFSIZ];
    char name[BUFSIZ];
    char path[BUFSIZ];
    while (fgets(line, sizeof line, file) != NULL)
    {
        sscanf(line, "%s\t%s", name, path); // e.g., Human \t hg38

        // remove leading digits (e.g., hg38 -> hg)
        for (size_t i = 0; i < strlen(name); ++i)
        {
            if (isdigit(name[i]))
            {
                name[i] = 0;
                break;
            }
        }

        if (first_entry)
        {
            reference_genome_name = name;
            reference_genome_path = path;
            first_entry = false;
        }
        else
        {
            genome_file_paths.emplace_back(name, path);
        }
    }
    fclose(file);
}

void run_annotate_with_msa(const std::string & gff_path, const AnnotateWithMSACLIParams & params,
                           const Model & model, const ScoreMSACLIParams & scoring_params,
                           std::unordered_set<std::string> & missing_sequences)
{
    constexpr bool debug_index_genomes                  = 0;
    constexpr bool debug_extract_cds                    = 0;
    constexpr bool debug_align_sequences                = 0;
    constexpr bool debug_transform_and_score_alignments = 1;

    const std::string dir_mmseqs_work = params.output_path;

    const std::string dir_genomesdb = dir_mmseqs_work + "/genomesDB";
    if (create_directory(dir_genomesdb))
        printf(OUT_INFO "Created the genomesDB directory.\n" OUT_RESET);

    const std::string dir_cds_files = dir_mmseqs_work + "/cds";
    if (create_directory(dir_cds_files))
        printf(OUT_INFO "Created the cds directory.\n" OUT_RESET);

    const std::string dir_tmp = dir_mmseqs_work + "/tmp";

    if (debug_extract_cds)
    {
        printf("Reading reference genome of GFF file %s ...\n", params.reference_genome_path.c_str());

        std::unordered_map<std::string, std::string> reference_genome;
        load_fasta_file(params.reference_genome_path, reference_genome);

        FILE *cds_fasta[3];
        for (uint8_t i = 0; i < 3; ++i)
        {
            const std::string cds_fasta_path = dir_cds_files + "/cds_phase" + std::to_string(i) + ".fasta";
            cds_fasta[i] = fopen(cds_fasta_path.c_str(), "w");
            if (cds_fasta[i] == NULL)
            {
                printf(OUT_ERROR "Could not create file %s.\n" OUT_RESET, cds_fasta_path.c_str());
                exit(-1);
            }
        }

        printf("Reading GFF file and extracting CDS coordinates ...\n");

        gff_reader reader(gff_path.c_str());
        reader.setup_progressbar(1, 1);

        gff_transcript t;
        uint32_t transcripts_left = 99999999; // TODO: remove

        std::unordered_set<std::string> processed_cds; // don't extract the same CDS sequence twice

        while (reader.get_next_transcript<false>(t) && transcripts_left > 0)
        {
            if (t.CDS.size() == 0)
                continue;

            const auto ref_chr = reference_genome.find(t.chr);
            if (ref_chr == reference_genome.end())
            {
                // only report a missing sequence once
                if (missing_sequences.find(t.chr) == missing_sequences.end())
                {
                    missing_sequences.insert(t.chr);
                    printf(OUT_DEL OUT_INFO "Sequence %s from the GFF file does not occur in the reference fasta file. Skipping ...\n" OUT_RESET, t.chr.c_str());
                }
                continue;
            }

            for (const auto & c : t.CDS)
            {
                // do not search a cds twice
                std::string key = t.chr + ":" + std::to_string(c.begin) + "-" + std::to_string(c.end) + "#" + t.strand;
                if (processed_cds.find(key) != processed_cds.end())
                    continue;

                processed_cds.emplace(key);

                std::string cds_seq = ref_chr->second.substr(c.begin - 1, c.end - (c.begin - 1));

                if (cds_seq.size() < 3ul + c.phase) // NOTE c.phase == 1 means that 1 base has to be thrown away because the previous CDS had 2 bases "too many"
                    continue;

                if (t.strand == '-')
                {
                    std::reverse(cds_seq.begin(), cds_seq.end());
                    for (uint64_t i = 0; i < cds_seq.size(); ++i)
                        cds_seq[i] = complement(cds_seq[i]);
                }

                // TODO: do not need end position
                fprintf(cds_fasta[c.phase], ">%s:%ld-%ld#%c\n%s\n", t.chr.c_str(), c.begin, c.end, t.strand, cds_seq.c_str());
            }

            --transcripts_left;
            reader.print_progress();
        }

        for (uint8_t i = 0; i < 3; ++i)
            fclose(cds_fasta[i]);
    }

    // do all the stuff wth mmseqs
    // TODO: keep/delete intermediary files?
    // TODO: --verbose to show or hide mmseqs output

    std::string cmd;

    // 1. index genomes
    if (debug_index_genomes)
    {
        std::string all_genomes_paths = "";
        for (uint64_t genome_id = 0; genome_id < params.aligning_genomes.size(); ++genome_id)
        {
            // the order is important, because the genomes/indices will be indexed in the same order by mmseqs
            all_genomes_paths += std::get<1>(params.aligning_genomes[genome_id]) + " ";
        }

        cmd = params.mmseqs2_bin + " createdb " + all_genomes_paths + " " + dir_genomesdb + "/genbankseqs";
        if (system_with_return(cmd.c_str()))
            exit(2);

        for (uint16_t i = 0; i < params.aligning_genomes.size(); ++i)
        {
            cmd = "bash -c $'" + params.mmseqs2_bin + " createsubdb <(awk \\'$3 == " + std::to_string(i) + "\\' " + dir_genomesdb + "/genbankseqs.lookup) " + dir_genomesdb + "/genbankseqs " + dir_genomesdb + "/genbankseqs_" + std::to_string(i) + "'";
            if (system_with_return(cmd.c_str()))
                exit(3);

            cmd = params.mmseqs2_bin + " createindex " + dir_genomesdb + "/genbankseqs_" + std::to_string(i) + " " + dir_tmp + " --search-type 2 --min-length 15" + " --threads " + std::to_string(params.threads);
            if (system_with_return(cmd.c_str()))
                exit(4);
        }
    }

    // load lookup table (seq id -> genome id)
    std::unordered_map<std::string, uint32_t> lookup_genome_ids;
    {
        const std::string lookup_path = dir_genomesdb + "/genbankseqs.lookup";
        FILE *lookup_file = fopen(lookup_path.c_str(), "r");
        if (lookup_file == NULL)
        {
            printf(OUT_ERROR "Could not open file %s.\n" OUT_RESET, lookup_path.c_str());
            exit(-1);
        }

        uint32_t lookup_id;
        char lookup_seq_name[255];
        uint32_t lookup_genome_id;
        while (!feof(lookup_file))
        {
            fscanf(lookup_file, "%u\t%s\t%u\n", &lookup_id, lookup_seq_name, &lookup_genome_id);
            lookup_genome_ids.emplace(lookup_seq_name, lookup_genome_id);
        }

        fclose(lookup_file);
    }

    for (uint8_t phase = 1; phase < 3; ++phase)
    {
        const std::string aln_dir = dir_mmseqs_work + "/aln_" + std::to_string(phase);
        if (create_directory(aln_dir))
            printf(OUT_INFO "Created the aln_dir directory.\n" OUT_RESET);

        const std::string cds_fasta_path = dir_cds_files + "/cds_phase" + std::to_string(phase) + ".fasta";
        const std::string exon_index_path = dir_cds_files + "/cds_phase" + std::to_string(phase) + ".index";
        const std::string aln_all_tophit_file = aln_dir + "/aln_all_tophit";
        const std::string msa_file = aln_dir + "/msa";
        const std::string maf_file = aln_dir + "/msa.maf";

        // index query sequences
        if (debug_align_sequences)
        {
            cmd = params.mmseqs2_bin + " createdb " + cds_fasta_path + " " + exon_index_path;
            if (system_with_return(cmd.c_str()))
                exit(5);

            std::string all_top_hit_files = "";

            for (uint16_t genome_id = 0; genome_id < params.aligning_genomes.size(); ++genome_id)
            {
                const std::string indexed_genome_path = dir_genomesdb + "/genbankseqs_" + std::to_string(genome_id);
                const std::string aln_output = aln_dir + "/aln_" + std::to_string(genome_id);
                const std::string aln_tophit_output = aln_dir + "/aln_tophit_" + std::to_string(genome_id);

                cmd = params.mmseqs2_bin + " search " + exon_index_path + " " + indexed_genome_path + " " + aln_output + " " + dir_tmp + " -a --search-type 4 --min-length 15 --remove-tmp-files --forward-frames " + std::to_string(phase + 1) + " --reverse-frames 0" + " --threads " + std::to_string(params.threads);
                if (system_with_return(cmd.c_str()))
                    exit(6);

                cmd = params.mmseqs2_bin + " filterdb " + aln_output + " " + aln_tophit_output + " --extract-lines 1" + " --threads " + std::to_string(params.threads);
                if (system_with_return(cmd.c_str()))
                    exit(7);
                all_top_hit_files = aln_tophit_output + " " + all_top_hit_files;
            }


            cmd = params.mmseqs2_bin + " mergedbs " + exon_index_path + " " + aln_all_tophit_file + " " + all_top_hit_files;
            if (system_with_return(cmd.c_str()))
                exit(8);

            cmd = params.mmseqs2_bin + " result2dnamsa " + exon_index_path + " " + dir_genomesdb + "/genbankseqs" + " " + aln_all_tophit_file + " " + msa_file + " --threads " + std::to_string(params.threads);
            if (system_with_return(cmd.c_str()))
                exit(9);
        }

        // score alignments
        if (debug_transform_and_score_alignments)
        {
            // TODO: cut off phase before even starting to score. when reading from the maf we have to consider this when computing the end-coordinate for lookup
            mmseqs_fasta_to_maf(msa_file, maf_file, params, lookup_genome_ids, phase);
            run_scoring_msa(maf_file, model, scoring_params, 1, 1);
        }
    }

    // TODO: don't store scores in output files but emplace them directly into the map
    // read in score files
    std::unordered_map<std::string, std::tuple<float, float> > computed_scores; // key (chrom:from-to#strand#phase), val: <phylo score, bls>
    for (uint8_t phase = 0; phase < 3; ++phase)
    {
        const std::string score_file_path = dir_mmseqs_work + "/aln_" + std::to_string(phase) + "/msa.maf.scores";

        FILE *score_file = fopen(score_file_path.c_str(), "r");
        if (score_file == NULL)
        {
            printf(OUT_ERROR "Could not open score file %s.\n" OUT_RESET, score_file_path.c_str());
            exit(-1);
        }

        char line[BUFSIZ];
        char chr[BUFSIZ];
        uint64_t start;
        uint64_t end;
        char strand;
        float score;

        fgets(line, sizeof line, score_file); // skip comment line
        fgets(line, sizeof line, score_file); // skip column name line

        while (fgets(line, sizeof line, score_file) != NULL)
        {
            sscanf(line, "%s\t%lu\t%lu\t%c\t%f", chr, &start, &end, &strand, &score);
            const std::string key = std::string(chr) + ':' + std::to_string(start) + '-' + std::to_string(end)
                                  + '#' + strand + '#' + std::to_string(phase);
            computed_scores.emplace(key, std::tuple<float, float>(score, NAN));
//            printf("%s\t%lu\t%lu\t%c\t%f\n", chr, start, end, strand, score);
        }
        fclose(score_file);
    }

    // read (again) and output with gff scores
    gff_reader gff_in(gff_path.c_str());
    gff_in.setup_progressbar(1, 1);

    // store output next to GFF file if no output_path was specified
    std::string output_file_path = params.output_path;
    if (output_file_path == "")
    {
        output_file_path = gff_path;
    }
    else
    {
        const size_t pos = gff_path.find_last_of('/');
        if (pos != std::string::npos) // gff_path == ".../file.gff"
            output_file_path += "/" + gff_path.substr(pos + 1);
        else                          // gff_path == "file.gff"
            output_file_path += "/" + gff_path;
    }

    // insert/append ".PhyloCSF++" (in)to filename
    const size_t output_file_path_last_dot = output_file_path.find_last_of('.');
    if (output_file_path_last_dot == std::string::npos)
        output_file_path += ".PhyloCSF++"; // NoFileEnding -> NoFileEnding.PhyloCSF++
    else
        output_file_path.insert(output_file_path_last_dot, ".PhyloCSF++"); // Human.Genes.gtf -> Human.Genes.PhyloCSF++.gtf

    FILE *gff_out = fopen(output_file_path.c_str(), "w");
    if (gff_out == NULL)
    {
        printf(OUT_ERROR "Error creating file %s!\n" OUT_RESET, output_file_path.c_str());
        exit(1);
    }

    // TODO: move header line to the end of the already existing header block
    fprintf(gff_out, "# PhyloCSF scores computed with PhyloCSF++ %s (%s, %s) and MMseqs2\n", ArgParse::version.c_str(), ArgParse::git_hash.c_str(), ArgParse::git_date.c_str());

    gff_transcript t;
    while (gff_in.get_next_transcript<true>(t))
    {
        uint64_t bases_with_scores = 0;
        float weighted_phylo_score = 0.0;
        float weighted_phylo_power = 0.0;

        if (t.CDS.size() > 0)
        {
            for (auto & c : t.CDS)
            {
                const std::string key = std::string(t.chr) + ':' + std::to_string(c.begin) + '-' + std::to_string(c.end)
                                        + '#' + t.strand + '#' + std::to_string(c.phase);

                auto scores_it = computed_scores.find(key);
                if (scores_it != computed_scores.end())
                {
                    const uint64_t cds_length = c.end - c.begin + 1;
                    bases_with_scores += cds_length;
                    c.phylo_score = std::get<0>(scores_it->second);
                    weighted_phylo_score += c.phylo_score / cds_length; // TODO: does this make sense?

                    if (scoring_params.comp_bls)
                    {
                        c.phylo_power = std::get<1>(scores_it->second);
                        weighted_phylo_power += c.phylo_power / cds_length; // TODO: does this make sense?
                    }
                }
            }

            if (bases_with_scores == 0)
            {
                t.phylo_score = NAN;
                t.phylo_power = NAN;
            }
            else
            {
                t.phylo_score = weighted_phylo_score / bases_with_scores;
                t.phylo_power = weighted_phylo_power / bases_with_scores;
            }
        }

        // output gtf file (with annotation)
        uint16_t CDS_id = 0;
        for (auto line_tuple : t.lines)
        {
            const feature_t f = std::get<0>(line_tuple);
            const std::string & line = std::get<1>(line_tuple);
            if (f == OTHER || t.CDS.size() == 0) // no annotation necessary
            {
                fprintf(gff_out, "%s\n", line.c_str());
            }
            else
            {
                float score;
                float power;
                if (f == TRANSCRIPT)
                {
                    score = t.phylo_score;
                    power = t.phylo_power;
                }
                else
                {
                    score = t.CDS[CDS_id].phylo_score;
                    power = t.CDS[CDS_id].phylo_power;
                    ++CDS_id;
                }

                // TODO: use correct format (gff/gtf)
                if (scoring_params.comp_bls)
                    fprintf(gff_out, "%s phylocsf_mean \"%.3f\"; phylocsf_power_mean \"%.3f\";\n", line.c_str(), score, power);
                else
                    fprintf(gff_out, "%s phylocsf_mean \"%.3f\";\n", line.c_str(), score);
            }
        }

        gff_in.print_progress();
    }

    fclose(gff_out);

    (void)gff_path;
    (void)missing_sequences;
    (void)model;
    (void)scoring_params;
}

int main_annotate_with_msa(int argc, char** argv)
{
    ArgParse args("phylocsf++ annotate-with-msa",
                  "Computes PhyloCSF scores for CDS features in GFF/GTF files and outputs them in \n"
                  "the same format. Requires MMseqs2 to be installed and reference genomes to compute multiple sequence"
                  "alignments from scratch.");

    AnnotateWithMSACLIParams params;
    ScoreMSACLIParams scoring_params; // parameters for the scoring part

    const std::string model_list = get_list_of_models();

    std::string default_strategy;
    switch (scoring_params.strategy)
    {
        case MLE:
            default_strategy = "MLE";
            break;
        case FIXED:
            default_strategy = "FIXED";
            break;
        case OMEGA:
            default_strategy = "OMEGA";
            break;
    }

    args.add_option("threads", ArgParse::Type::INT, "Parallelize scoring of multiple MSAs in a file using multiple threads. Default: " + std::to_string(params.threads), ArgParse::Level::GENERAL, false);
    args.add_option("output", ArgParse::Type::STRING, "Path where tracks in wig format will be written. If it does not exist, it will be created. Default: output files are stored next to the input files.", ArgParse::Level::GENERAL, true);

    args.add_option("strategy", ArgParse::Type::STRING, "PhyloCSF scoring algorithm: MLE, FIXED or OMEGA. Default: " + default_strategy, ArgParse::Level::GENERAL, false);
    args.add_option("comp-power", ArgParse::Type::BOOL, "Output confidence score (branch length score). Default: " + std::to_string(scoring_params.comp_bls), ArgParse::Level::GENERAL, false);
    args.add_option("mmseqs-bin", ArgParse::Type::STRING, "Path to MMseqs2 binary. Default: " + params.mmseqs2_bin, ArgParse::Level::GENERAL, false);

    args.add_positional_argument("genome-file", ArgParse::Type::STRING, "Two-column text file with species name and path to its genomic fasta file. First line has to be the reference genome of the GFF file.", true);
    args.add_positional_argument("model", ArgParse::Type::STRING, "Path to parameter files, or one of the following predefined models: " + model_list + ".", true);
    args.add_positional_argument("gff-files", ArgParse::Type::STRING, "One or more GFF/GTF files with coding exons to be scored.", true, true);
    args.parse_args(argc, argv);

    args.check_args(); // check here because if "--model-info" is set, we don't want to require mandatory arguments

    // check whether MMseqs2 is installed
    printf(OUT_INFO "Checking whether MMseqs2 is installed ...\n" OUT_RESET);
    if (args.is_set("mmseqs-bin"))
        params.mmseqs2_bin = args.get_string("mmseqs-bin");

    if (system_with_return(params.mmseqs2_bin + " --help", false))
    {
        printf("ERROR: MMseqs2 seems not to be installed. `%s --help` does not seem to work.\n", params.mmseqs2_bin.c_str());
        return -1;
    }

    if (system_with_return(params.mmseqs2_bin + " result2dnamsa --help", false))
    {
        printf("ERROR: A version of MMseqs2 seems to be installed that is too old. `%s result2dnamsa --help` does not seem to work.\n", params.mmseqs2_bin.c_str());
        return -1;
    }

    // retrieve flags/options that were set
    if (args.is_set("threads"))
        params.threads = args.get_int("threads");
    scoring_params.threads = params.threads;

    if (args.is_set("strategy"))
    {
        std::string strategy = args.get_string("strategy");
        std::transform(strategy.begin(), strategy.end(), strategy.begin(), toupper);

        if (strategy == "MLE")
            scoring_params.strategy = MLE;
        else if (strategy == "FIXED")
            scoring_params.strategy = FIXED;
        else if (strategy == "OMEGA")
            scoring_params.strategy = OMEGA;
        else
        {
            printf(OUT_ERROR "Please choose a valid strategy (MLE, FIXED or OMEGA)!\n" OUT_RESET);
            return -1;
        }
    }

    if (args.is_set("output"))
    {
        params.output_path = args.get_string("output");
        // create a directory if it doesn't exist yet
        if (create_directory(params.output_path))
            printf(OUT_INFO "Created the output directory.\n" OUT_RESET);
    }
    scoring_params.output_path = ""; // this will automatically save scores next to maf file

    scoring_params.comp_phylo = true;
    scoring_params.comp_anc = false;
    if (args.is_set("comp-power"))
        scoring_params.comp_bls = args.get_bool("comp-power");
    else
        scoring_params.comp_bls = false;

    const std::string genome_file = args.get_positional_argument("genome-file");
    load_genome_file(genome_file, params.aligning_genomes, params.reference_genome_name, params.reference_genome_path);

    // TODO: test whether species in genome-file can be mapped to model.

    std::string species = ""; // TODO: set species subset according to provided reference genomes in genome-file

    // load and prepare model
    Model model;
    load_model(model, args.get_positional_argument("model"), species, false, 0, "");

    // run for every gff file
    std::unordered_set<std::string> missing_sequences;
    for (uint16_t i = 2; i < args.positional_argument_size(); ++i)
    {
        if (args.positional_argument_size() > 2)
            printf("Processing GFF %s\n", args.get_positional_argument(i).c_str());

        run_annotate_with_msa(args.get_positional_argument(i), params, model, scoring_params, missing_sequences);

        if (args.positional_argument_size() > 2)
            printf("----------------------------\n");
    }

    printf(OUT_DEL "Done!\n");

    return 0;
}
