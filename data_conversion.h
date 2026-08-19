#ifndef data_conversion_h
#define data_conversion_h

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <filesystem>
#include <system_error>
#include <cstdlib>
#include <cctype>
using namespace std;

// Opens filename for writing, creating its containing folder first. ofstream
// does not create missing directories - it just fails to open - and writes to a
// failed stream are silent no-ops, so without this a bad path produces no file
// and no error while the caller happily reports success.
inline ofstream openForWrite(const string &filename) {
    filesystem::path outPath(filename);
    if (outPath.has_parent_path()) {
        error_code ec;
        filesystem::create_directories(outPath.parent_path(), ec);
        if (ec) {
            cerr << "ERROR: could not create " << outPath.parent_path().string()
                 << ": " << ec.message() << "\n";
            exit(EXIT_FAILURE);
        }
    }

    ofstream outputFile(filename);
    if (!outputFile) {
        cerr << "ERROR: could not open " << filename << " for writing\n";
        exit(EXIT_FAILURE);
    }
    return outputFile;
}

inline void state2text(const vector<vector<double>> &U, const string &filename){
    ofstream outputFile = openForWrite(filename);
    for(int iElem = 0; iElem < U.size(); iElem++){
        for(int iU = 0; iU < 4; iU++){
            outputFile << U[iElem][iU] << " ";
        }
        outputFile << "\n";
    }
}

// The second-order solver carries the state as one contiguous [nElem*4] block
// rather than a vector of 4-element vectors, so it gets its own overload here
// instead of being unflattened just to be written out. Same file format.
inline void state2text(const vector<double> &U, const string &filename){
    ofstream outputFile = openForWrite(filename);
    for(size_t iElem = 0; iElem * 4 < U.size(); iElem++){
        for(int iU = 0; iU < 4; iU++){
            outputFile << U[iElem * 4 + iU] << " ";
        }
        outputFile << "\n";
    }
}

// Conversions between the two layouts, used at the seam between the first-order
// solver (nested) and the second-order solver (flat).
inline vector<double> flattenState(const vector<vector<double>> &U){
    vector<double> flat(U.size() * 4);
    for(size_t iElem = 0; iElem < U.size(); iElem++){
        for(int iU = 0; iU < 4; iU++){
            flat[iElem * 4 + iU] = U[iElem][iU];
        }
    }
    return flat;
}

inline vector<vector<double>> unflattenState(const vector<double> &U){
    vector<vector<double>> nested(U.size() / 4, vector<double>(4));
    for(size_t iElem = 0; iElem < nested.size(); iElem++){
        for(int iU = 0; iU < 4; iU++){
            nested[iElem][iU] = U[iElem * 4 + iU];
        }
    }
    return nested;
}

// Inverse of state2text: reads a [nElem x 4] state file into U. The element
// count is not stored in the file, so rows are appended until EOF.
inline vector<vector<double>> text2state(const string &filename){
    ifstream inputFile(filename);
    if (!inputFile) {
        cerr << "ERROR: could not open " << filename << " for reading\n";
        exit(EXIT_FAILURE);
    }

    vector<vector<double>> U;
    string line;
    while(getline(inputFile, line)){
        istringstream lineStream(line);
        vector<double> state(4);
        for(int iU = 0; iU < 4; iU++){
            if(!(lineStream >> state[iU])){
                // Trailing blank line at EOF is fine; a partial row is not.
                if(iU == 0) break;
                cerr << "ERROR: " << filename << " line " << U.size() + 1
                     << " has fewer than 4 values\n";
                exit(EXIT_FAILURE);
            }
        }
        if(lineStream.fail()) continue;
        U.push_back(state);
    }
    return U;
}

// Recovers the flux option and CFL from a filename written by the naming scheme
// in main ("dat/<flux>_CFL<value>_firstOrder.dat"). The state file itself stores
// only the four conserved variables, so the filename is the only record of which
// settings produced it. Returns false - leaving opt and CFL untouched - if the
// name does not follow the scheme, so the caller can ask the user instead of
// silently marching the second-order solve with the wrong flux.
inline bool filename2settings(const string &filename, int &opt, double &CFL){
    static const char* fluxName[] = { "", "roe", "rusanov", "hlle" };

    // stem() drops both the directory and the ".dat", leaving "roe_CFL0.7_firstOrder"
    string stem = filesystem::path(filename).stem().string();
    for(size_t i = 0; i < stem.size(); i++){
        stem[i] = static_cast<char>(tolower(static_cast<unsigned char>(stem[i])));
    }

    int optParsed = 0;
    string flux = stem.substr(0, stem.find('_'));
    for(int i = 1; i <= 3; i++){
        if(flux == fluxName[i]) optParsed = i;
    }

    double cflParsed = 0.0;
    size_t cflPos = stem.find("cfl");
    if(cflPos != string::npos){
        // extraction stops at the '_' before "firstorder", so no trimming is needed
        istringstream cflStream(stem.substr(cflPos + 3));
        if(!(cflStream >> cflParsed)) cflParsed = 0.0;
    }

    if(optParsed == 0 || cflParsed <= 0) return false;
    opt = optParsed;
    CFL = cflParsed;
    return true;
}

inline void L1Residual2text(vector<double> &L1Residual, const string &filename){
    ofstream outputFile = openForWrite(filename);
    for(int iR = 0; iR < L1Residual.size(); iR++){
        outputFile << L1Residual[iR] << "\n";
    }
}

#endif /* data2text_h */
