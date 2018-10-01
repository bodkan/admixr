#include <iostream>
#include <sstream>
#include <fstream>
#include <regex>

#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(BH)]]

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>

std::vector<std::string> split_line(std::string& line)
{
    std::istringstream ss(line);
    std::vector<std::string> tokens{
        std::istream_iterator<std::string>{ss},
        std::istream_iterator<std::string>{}
    };
    return tokens;
}

void write_snp(std::ofstream& file, const std::vector<std::string>& tokens)
{
    file << tokens[2] << "\t" << // ID
            tokens[0] << "\t" << // chromosome
            "0.0"     << "\t" << // genetic distance
            tokens[1] << "\t" << // position
            tokens[3] << "\t" << // 
            tokens[4] << "\n";
}

void write_ind(std::ofstream& file, const std::vector<std::string> inds)
{
    for (auto const & ind : inds) file << ind << "\tU\t" << ind << "\n";
}

// test case for VCF GT -> EIGENSTRAT conversion
/* std::vector<std::string> gts{".|.", "./.", ".", "0|0", "0/0", "0", "0|1", "1|0", "0/1", "1|1", "1/1", "1", ".|.", "./.", ".", "0|0", "0/0", "0", "0|1", "1|0", "0/1", "1|1", "1/1", "1"}; */
/* std::string str_gts = std::accumulate( */
/*    gts.begin(), */
/*    gts.end(), */
/*    std::string("GT"), */
/*    [](const std::string & a, const std::string & b) { return a + "\t" + b; } */
/* ); */
/* std::cout << str_gts << "\n"; */
/* std::cout << convert_genotypes(str_gts) << "\n"; */
std::string convert_genotypes(const std::string& s)
{
    std::string converted(s);

    converted = std::regex_replace(converted, std::regex("\t0(\t|$)"),         "\t2$1");
    converted = std::regex_replace(converted, std::regex("\t1(\t|$)"),         "\t0$1");
    converted = std::regex_replace(converted, std::regex("\t\\./\\.(\t|$)"),   "\t9$1");

    converted = std::regex_replace(converted, std::regex("\\.\\|\\."), "9");
    converted = std::regex_replace(converted, std::regex("\\."),       "9");

    converted = std::regex_replace(converted, std::regex("0\\|0"),     "2");
    converted = std::regex_replace(converted, std::regex("0/0"),       "2");

    converted = std::regex_replace(converted, std::regex("0\\|1"),     "1");
    converted = std::regex_replace(converted, std::regex("1\\|0"),     "1");
    converted = std::regex_replace(converted, std::regex("0/1"),       "1");

    converted = std::regex_replace(converted, std::regex("1\\|1"),     "0");
    converted = std::regex_replace(converted, std::regex("1/1"),       "0");
    return converted;
}

void write_geno(std::ofstream& file, const std::string& line)
{
    // extract part of the line after the GT format field
    std::string l = line.substr(line.rfind("GT"));
    // remove the GT field and keep just the individual genotypes
    l = l.substr(l.find("\t") + 1);
    // convert the genotypes on the line and concatenate them
    std::string converted = convert_genotypes(l);
    file << std::regex_replace(converted, std::regex("\t"), "") << "\n";
}

//' Convert VCF file into an EIGENSTRAT format
//'
//' @param vcf A path to a VCF file.
//' @param eigenstrat A path to an output EIGENSTRAT triplet.
//'
//' @export
//'
// [[Rcpp::export]]
void vcf_to_eigenstrat(const char* vcf, const char* eigenstrat) {
    std::ofstream ind_file(std::string(eigenstrat) + ".ind"),
                  snp_file(std::string(eigenstrat) + ".snp"),
                  geno_file(std::string(eigenstrat) + ".geno");

    std::unique_ptr<std::istream> vcf_file;

    // setup boost machinery that will be used for reading gzipped input
    boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
    std::ifstream file(vcf, std::ios_base::in | std::ios_base::binary);
    in.push(boost::iostreams::gzip_decompressor());
    in.push(file);

    // check whether the VCF is gzipped and setup the input stream accordingly
    if (std::regex_search(vcf, std::regex(".vcf.gz"))) {
        vcf_file = std::make_unique<std::istream>(&in);
    } else {
        vcf_file = std::make_unique<std::ifstream>(vcf);
    }

    std::string line;

    // parse the header part of the VCF
    while (std::getline(*vcf_file, line)) {
        // skip the header
        if (line.find("##") == 0) continue;

        // write the ind file
        if (line.find("#") == 0) {
            auto elems = split_line(line);
            auto inds = std::vector<std::string>(elems.begin() + 9, elems.end());
            write_ind(ind_file, inds);
            break;
        }
    }

    // parse the genotype part of the VCF
    while (std::getline(*vcf_file, line)) {
        auto elems = split_line(line);

        // write snp and geno records
        if (line.find("#") != 0) {
            write_snp(snp_file, elems);
            write_geno(geno_file, line);
        }
    }
}

int main(int argc, char* argv[])
{
    if (argc != 3) {
        std::cout << "Usage:\n\t./vcf_to_eigenstrat <path to VCF> <output EIGENSTRAT prefix>\n";
        return 0;
    }
    vcf_to_eigenstrat(argv[1], argv[2]);
    return 0;
}

