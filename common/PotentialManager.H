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
class Matrix;

class PotentialManager
{
	public:
		PotentialManager();
		~PotentialManager();

		int init();
		int setEvidenceManager(EvidenceManager*);
		int setRandom(bool);
		int reset();

		Error::ErrorCode populatePotentialsSlimFactors(map<int,SlimFactor*>&,VSET&);

		int populatePotential(Potential*);

	private:
		int initData();
		int estimateCovariance(int,int);

		EvidenceManager* evMgr;
		gsl_matrix* ludecomp;
		gsl_permutation* perm;
		bool randomData;
		Matrix* data;
		Matrix* meanMat;
		Matrix* covMat;
};
#endif
