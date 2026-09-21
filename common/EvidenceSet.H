#ifndef _EVIDENCE_SET
#define _EVIDENCE_SET

#include <vector>
#include <unordered_map>

using namespace std;

class EvidenceSet
{
    public:
        EvidenceSet(vector<double>& dataset, size_t sampleCount);

        int getSampleCount();

        int getVariableCount();

        double getEvidenceAt(int varIndex, int sampleIndex);

    private:
        vector<double>* dataset;
        size_t sampleCount;
        size_t varCount;

};

#endif