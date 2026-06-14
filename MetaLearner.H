#ifndef _META_LEARNER_
#define _META_LEARNER_

#include <map>
#include <unordered_map>
#include "CommonTypes.H"
#include "HierarchicalCluster.H"

using namespace std;

//This is gamma map for each variable within each datapoint
typedef map<int,INTDBLMAP*> DATA_GAMMA_MAP;
typedef map<int,INTDBLMAP*>::iterator DATA_GAMMA_MAP_ITER;

class VariableManager;
class EvidenceManager;
class PotentialManager;
class MetaMove;
class HierarchicalCluster;
class Matrix;

class MetaLearner
{
	public:
		MetaLearner();
		~MetaLearner();
		int setMaxFactorSize_Approx(int);
		int setBeta1(double aval);
		int setBeta_Motif(double bval);
		int setConvergenceThreshold(double);
		int setRestrictedList(const char*);
		int setPreRandomizeSplit();
		int initEdgePriorMeta(const string& priorName, map<string, map<string, double>*>& graph, unordered_map<int, unordered_map<int, double>*>& edgePrior);
		void doCrossValidation(int);
		int setGlobalEvidenceManager(EvidenceManager*);
		int setVariableManager(VariableManager*);
		int setOutputDirName(const char*);
		int setClusteringThreshold(double);
		int setSpecificFold(int);
		int getPredictionError_CrossValid(int);
		int setRandom(bool);
		int readModuleMembership(const char*);
		int setDefaultModuleMembership();
		int setPriorGraph_All(const char*);

	private:
		void start(int);
		int initEdgePriorMeta_All();
		double getInitPLLScore();
		double getInitPLLScore(int vId);
		double getPLLScore();
		double getEdgePrior(int tfID, int targetID);
		int clearFoldSpecData();
		int initEdgeSet();
		bool getNextMove(int currK, int vID, MetaMove& outMove);
		void makeMove(MetaMove& move, int currIteration);
		int dumpAllGraphs(int,int,int);
		int setPriorGraph(const char* aFName, map<string, map<string, double>*>& priorGraph);
		void redefineModules(int currFold);
		double getModuleContribLogistic(string& tgtName, string& tfName);
		void initPhysicalDegree();
		void getEnrichedTFs(map<string, int>& tfSet, map<string, int>* genes, map<string, map<string, double>*>& edgeSet);
		void initDistances();
		void updateSharedParentDistances();

		FactorGraph* factorGraph;
		PotentialManager* potManager;
		EvidenceManager* evidenceManager;
		VariableManager* varManager;
		map<string, double> betaMap;
		map<string, unordered_map<int, unordered_map<int, double>*>*> edgePriorMap;
		map<string, map<string, map<string, double>*>*> priorGraphMap;
		int maxFactorSizeApprox;
		double convThreshold;
		char restrictedFName[1024];
		unordered_map<int, unordered_map<int, int>> edgeMap;
		map<string,int> restrictedVarList;
		unordered_map<int, double>* currPLL;
		map<int,double> finalScores;
		bool preRandomizeSplit;
		char outputDirName[1024];
		map<string, double> edgePresenceProb;
		map<int, double> varNeighborhoodPrior;
		double beta1;
		double beta_motif;
		bool random;
		map<int, map<string, int>*> moduleGeneSet;
		unordered_map<string, int> geneModuleID;
		unordered_map<int, unordered_map<string, int>*> moduleIndegree;
		unordered_map<string, int> regulatorModuleOutdegree;
		double clusterThreshold;
		unordered_map<string, int> variableStatus;
		char foldoutDirName[1024];
		int specificFold;
		HierarchicalCluster hc;
		Matrix* correlationDistances;
		Matrix* sharedParentDistances;
		unordered_map<int, unordered_map<int, int>> sharedParents;
		vector<int> updatedThisIteration;
};
#endif
