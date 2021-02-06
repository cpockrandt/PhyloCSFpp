#include "common.hpp"
#include "arg_parse.hpp"
#include "gff_reader.hpp"

#include <omp.h>

#include <unordered_map>

#include <stdio.h>
#include <sys/wait.h>

struct AnnotateWithMSACLIParams
{
    std::string output_path = "";
    std::string mmseqs2_bin = "mmseqs";
    unsigned threads = omp_get_max_threads();
    algorithm_t strategy = MLE;
    bool comp_bls = true;

    std::unordered_map<std::string, std::string> genome_file_paths;
    std::string reference_genome_name;
};

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

void load_genome_file(const std::string & genome_file, std::unordered_map<std::string, std::string> & genome_file_paths,
                      std::string & reference_genome_name)
{
    FILE *file = fopen(genome_file.c_str(), "r");
    if (file == NULL)
    {
        printf(OUT_ERROR "Cannot open genome-file %s\n" OUT_RESET, genome_file.c_str());
        exit(-1);
    }

    bool first_line = true;
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

        if (first_line)
        {
            reference_genome_name = name;
            first_line = false;
        }
        genome_file_paths.emplace(name, path);
    }
    fclose(file);
}

void run_annotate_with_msa(const std::string & gff_path, const AnnotateWithMSACLIParams & params,
                           const uint32_t file_id, const uint32_t files)
{
    printf("Reading reference genome of GFF file ...\n");

    const std::string & reference_path = params.genome_file_paths.find(params.reference_genome_name)->second;
    std::unordered_map<std::string, std::string> reference_genome;
    load_fasta_file(reference_path, reference_genome);

    printf("Reading GFF file and extracting CDS coordinates ...\n");

    gff_reader reader(gff_path.c_str());
    reader.setup_progressbar(file_id, files);

    uint32_t transcripts_left = 15;

    gff_transcript t;
    while (reader.get_next_transcript<false>(t) && transcripts_left > 0)
    {
        if (t.CDS.size() == 0 || t.strand == '+' /*|| t.CDS.size() != 1 || t.CDS[0].end - t.CDS[0].begin > 1500*/)
            continue;

        const auto ref_chr = reference_genome.find(t.chr);
        if (ref_chr == reference_genome.end())
        {
            // TODO: warnings
            continue;
        }

        printf("%s\t%ld\t%ld\t%c\n", t.chr.c_str(), t.begin, t.end, t.strand);
        for (const auto & c : t.CDS)
        {
            const std::string cds_seq = ref_chr->second.substr(c.begin - 1, c.end - (c.begin - 1));
            printf("\t%ld\t%ld\t%s\n", c.begin, c.end, cds_seq.c_str());
        }

        printf("\n\n\n");

        --transcripts_left;
//        reader.print_progress();
    }

    // do all the stuff wth mmseqs
    // score alignments
    // read, copy and annotate gff

    (void)file_id;
    (void)files;
}

int main_annotate_with_msa(int argc, char** argv)
{
    ArgParse args("phylocsf++ annotate-with-msa",
                  "Computes PhyloCSF scores for CDS features in GFF/GTF files and outputs them in \n"
                  "the same format. Requires MMseqs2 to be installed and reference genomes to compute multiple sequence"
                  "alignments from scratch.");

    AnnotateWithMSACLIParams params;

    const std::string model_list = get_list_of_models();

    std::string default_strategy;
    switch (params.strategy)
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
    args.add_option("output", ArgParse::Type::STRING, "Path where tracks in wig format will be written. If it does not exist, it will be created. Default: output files are stored next to the input files.", ArgParse::Level::GENERAL, false);

    args.add_option("strategy", ArgParse::Type::STRING, "PhyloCSF scoring algorithm: MLE, FIXED or OMEGA. Default: " + default_strategy, ArgParse::Level::GENERAL, false);
    args.add_option("comp-power", ArgParse::Type::BOOL, "Output confidence score (branch length score). Default: " + std::to_string(params.comp_bls), ArgParse::Level::GENERAL, false);
    args.add_option("mmseqs-bin", ArgParse::Type::BOOL, "Path to MMseqs2 binary. Default: " + params.mmseqs2_bin, ArgParse::Level::GENERAL, false);

    args.add_positional_argument("genome-file", ArgParse::Type::STRING, "Two-column text file with species name and path to its genomic fasta file. First line has to be the reference genome of the GFF file.", true);
    args.add_positional_argument("model", ArgParse::Type::STRING, "Path to parameter files, or one of the following predefined models: " + model_list + ".", true);
    args.add_positional_argument("gff-files", ArgParse::Type::STRING, "One or more GFF/GTF files with coding exons to be scored.", true, true);
    args.parse_args(argc, argv);

    args.check_args(); // check here because if "--model-info" is set, we don't want to require mandatory arguments

    // phylocsf++ annotate --genomes-file PATH --strategy STRING --mmseqs-bin PATH <gtf file(s)>
    // both: --output, --confidence

    // check whether MMseqs2 is installed
    if (args.is_set("mmseqs-bin"))
        params.mmseqs2_bin = args.get_int("mmseqs-bin");

    if (system_with_return(params.mmseqs2_bin + " --help > /dev/null 2>&1"))
    {
        printf("ERROR: MMseqs2 seems not to be installed. `%s --help` does not seem to work.\n", params.mmseqs2_bin.c_str());
        return -1;
    }

    if (system_with_return(params.mmseqs2_bin + " result2dnamsa --help > /dev/null 2>&1"))
    {
        printf("ERROR: A version of MMseqs2 seems to be installed that is too old. `%s result2dnamsa --help` does not seem to work.\n", params.mmseqs2_bin.c_str());
        return -1;
    }

    // retrieve flags/options that were set
    if (args.is_set("threads"))
        params.threads = args.get_int("threads");

    if (args.is_set("output"))
    {
        params.output_path = args.get_string("output");
        // create a directory if it doesn't exist yet
        if (create_directory(params.output_path))
            printf(OUT_INFO "Created the output directory.\n" OUT_RESET);
    }

    const std::string genome_file = args.get_positional_argument("genome-file");
    load_genome_file(genome_file, params.genome_file_paths, params.reference_genome_name);

    // TODO: test whether species in genome-file can be mapped to model.

    std::string species = ""; // TODO: set species subset according to provided reference genomes in genome-file

    // load and prepare model
    Model model;
    load_model(model, args.get_positional_argument("model"), species, false, 0, "");

    // run for every gff file
    for (uint16_t i = 2; i < args.positional_argument_size(); ++i)
    {
         run_annotate_with_msa(args.get_positional_argument(i), params, i - 1, args.positional_argument_size() - 2);
    }

    printf(OUT_DEL "Done!\n");

    return 0;
}
