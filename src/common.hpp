#pragma once

#include <cassert>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstring>

#include <string>
#include <unordered_map>

#define OUT_INFO    "\033[1;34m"
#define OUT_ERROR   "\033[1;31m"
#define OUT_RESET   "\033[0m"
#define OUT_DEL     "\x1b[K"
#define OUT_BOLD    "\033[1m"

void append(FILE * file_dest, const char path_src[])
{
    FILE *file_src = fopen(path_src, "rb");
    if (!file_src)
        abort();

    char buf[BUFSIZ];
    size_t n;
    while ((n = fread(buf, 1, sizeof buf, file_src)) > 0)
        if (fwrite(buf, 1, n, file_dest) != n)
            abort();

    if (ferror(file_src))
        abort();

    fclose(file_src);
}

bool create_directory(const std::string & path)
{
    struct stat st;

    if (stat(path.c_str(), &st) == -1)
    {
        mkdir(path.c_str(), 0764);
        return true;
    }
    return false;
}

void my_fprintf(FILE* f, const char *format_str, const float d)
{
    char buf[12];
    sprintf(buf, format_str, d);
    for (int8_t i = strlen(buf); i >= 0; --i)
    {
        if (buf[i] == '.')
        {
            buf[i + 1] = '0';
            break;
        }
        if (isdigit(buf[i]))
        {
            if (buf[i] != '0')
                break;
            else
                buf[i] = 0; // cut off the 0
        }
    }
    fprintf(f, "%s\n", buf);
}

int system_with_return(std::string cmd, bool verbose = true)
{
    if (verbose)
        printf(OUT_INFO "%s\n" OUT_RESET, cmd.c_str());
    else
        cmd  += " > /dev/null 2>&1";
    int ret = system(cmd.c_str());
    return WEXITSTATUS(ret);
}

void str_to_lower(char * str)
{
    for (size_t i = 0; str[i]; ++i)
        str[i] = tolower(str[i]);
}

void str_to_lower(std::string & str)
{
    for (size_t i = 0; i < str.size(); ++i)
        str[i] = std::tolower(str[i]);
}

void str_to_upper(std::string & str)
{
    for (size_t i = 0; i < str.size(); ++i)
        str[i] = std::toupper(str[i]);
}

bool is_gff_format(const std::string & line)
{
    uint8_t col = 1;
    size_t line_pos = 0;
    while (col > 0 && line_pos < line.size())
    {
        if (col == 9)
        {
            while (line_pos < line.size())
            {
                if (line[line_pos] == ' ') // key "value"
                    return false;
                else if (line[line_pos] == '=') // key=vale
                    return true;
                ++line_pos;
            }
        }

        if (line[line_pos] == '\t')
            ++col;

        ++line_pos;
    }
    return true; // assume by default it's a gff file
}

void load_fasta_file(const std::string & path_to_fasta, std::unordered_map<std::string, std::string> & genome)
{
    FILE *file = fopen(path_to_fasta.c_str(), "r");
    if (file == NULL)
    {
        printf(OUT_ERROR "Cannot open genomic fasta file %s\n" OUT_RESET, path_to_fasta.c_str());
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
                genome.emplace(id, seq);

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