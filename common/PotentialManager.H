/* This class calculates and caches the means and covariances for each variable
present in the evidences of the provided EvidenceManager.*/

#ifndef _POTENTIAL_MANAGER
#define _POTENTIAL_MANAGER

#include <vector>
#include <map>

using namespace std;

class Potential;
class EvidenceManager;
class SlimFactor;

class PotentialManager
{
	public:
		PotentialManager();
		~PotentialManager();

		int init(int);
		int setEvidenceManager(EvidenceManager*);
		int setOutputDir(const char*);
		int setRandom(bool);
		int reset();

		Error::ErrorCode populatePotentialsSlimFactors(map<int,SlimFactor*>&,VSET&);

		int populatePotential(Potential*);

	private:
		int readAllMeanCov(const char*, const char*);
		int estimateAllMeanCov(INTDBLMAP&,map<int,INTDBLMAP*>&,INTINTMAP&);
		int estimateCovariance(INTDBLMAP*,int,int);

		EvidenceManager* evMgr;
		INTDBLMAP globalMean;
		map<int,INTDBLMAP*> globalCovar;
		char outputDir[1024];
		gsl_matrix* ludecomp;
		gsl_permutation* perm;
		bool randomData;
};
#endif
