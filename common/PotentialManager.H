#ifndef _POTENTIAL_MANAGER
#define _POTENTIAL_MANAGER

using namespace std;

class Potential;
class EvidenceManager;
class Matrix;

class PotentialManager
{
	public:
		PotentialManager();
		~PotentialManager();

		// When init is called PotentialManager clears its internal data and loads
		// new means and covariances from the provided EvidenceManager.
		int init(EvidenceManager* evMgr, bool random, vector<int>& varIDs);

		// Creates a Potential of a single variable.
		Potential* createPotential(int factorID);

		// Computes the conditional Gaussian log likelihood of the factor's assignments
		// conditioned upon its parents. It also creates a Potential representing this conditional
		// Gaussian PDF.
		double computeLL(int factorID, vector<int>& parentIDs, int sampleSize, Potential** newPot);

		void computeLLs(int factorID, int sampleSize, vector<int>& existingParents, vector<int>&candidateParents, unordered_map<int, double>&scores);

		void computeLLs(int factorID, int sampleSize, vector<int>&candidateParents, unordered_map<int, double>&scores);

	private:
		gsl_matrix* ludecomp;
		gsl_permutation* perm;
		vector<double> globalMeans;
		Matrix *globalCovariances;
};
#endif
