#ifndef data2text_h
#define data2text_h

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <filesystem>
#include <system_error>
#include <cstdlib>
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

inline void state2text(vector<vector<double>> &U, const string &filename){
    ofstream outputFile = openForWrite(filename);
    for(int iElem = 0; iElem < U.size(); iElem++){
        for(int iU = 0; iU < 4; iU++){
            outputFile << U[iElem][iU] << " ";
        }
        outputFile << "\n";
    }
}

inline void L1Residual2text(vector<double> &L1Residual, const string &filename){
    ofstream outputFile = openForWrite(filename);
    for(int iR = 0; iR < L1Residual.size(); iR++){
        outputFile << L1Residual[iR] << "\n";
    }
}

#endif /* data2text_h */
