#ifndef _FACTOR_MANAGER
#define _FACTOR_MANAGER

#include <vector>
#include <fstream>
using namespace std;
#include "gsl/gsl_randist.h"
#include "CommonTypes.H"
#define SAMPLE_CNT 1000

class Error;
class VariableManager;
class EvidenceManager;
class PotentialManager;
class SlimFactor;
class LatticeStructure;
class Graph;
class FactorGraph;

class FactorManager
{
	public:
		FactorManager();
		~FactorManager();

		int setVariableManager(VariableManager*);
		int setEvidenceManager(EvidenceManager*);
		int setPotentialManager(PotentialManager*);
		int setMaxFactorSize(int);
		int setMaxFactorSize_Approx(int);
		int setOutputDir(const char*);
		
		int allocateFactorSpace();

		FactorGraph* createInitialFactorGraph();

		int learnStructure();
		int showStructure(int);
		int readStructure(int);

		SlimFactor* getFactorAt(int);

		int readRestrictedVarlist(const char*);

	private:	
		int initFactorSet();
		int populateFactorSet();
		int combCnt(int,int);
		int combCnt(int,int,int);
		int getFactorIndex(int*, int);
		int getFactorKey(int*, int,string&);
		int addToLattice(SlimFactor*,int**);
		int readClusterProperties(const char*);
		Error::ErrorCode estimateClusterProperties();
		SlimFactor* getFactorFromVars(int*, int);

		int globalFactorID;	
		int maxFactorSize;
		int maxFactorSize_Approx;
		VariableManager* vMgr;
		PotentialManager* potMgr;
		EvidenceManager* evMgr;

		//This is used to store the ids of factors which have only one variable
		//This will be used to compute the multi-information
		map<int,SlimFactor*> slimFactorSet;

		map<string,int> factorNameToIDMap;
		map<int,string> factorIDToNameMap;

		LatticeStructure lattice;

		char outputDir[1024];

		map<int,Variable*> restrictedNeighborList;
};
#endif
