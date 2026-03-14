#include <iostream>
#include <cstring>
#include <math.h>
#include <sys/timeb.h>
#include <sys/time.h>
#include <time.h>

#include "Error.H"
#include "Variable.H"
#include "SlimFactor.H"
#include "Potential.H"
#include "Evidence.H"
#include "EvidenceManager.H"
#include "VariableManager.H"
#include "LatticeStructure.H"
#include "PotentialManager.H"
#include "FactorGraph.H"
#include "FactorManager.H"

FactorManager::FactorManager()
{
	globalFactorID=0;
	maxFactorSize_Approx=-1;
}

FactorManager::~FactorManager()
{
	for(map<int,SlimFactor*>::iterator fIter=slimFactorSet.begin();fIter!=slimFactorSet.end();fIter++)
	{
		delete fIter->second;
	}
	slimFactorSet.clear();
	lattice.clear();
	factorNameToIDMap.clear();
	factorIDToNameMap.clear();
}

int 
FactorManager::setVariableManager(VariableManager* aPtr)
{
	vMgr=aPtr;
	return 0;
}

int 
FactorManager::setEvidenceManager(EvidenceManager* aPtr)
{
	evMgr=aPtr;
	return 0;
}

int
FactorManager::setPotentialManager(PotentialManager* aPtr)
{
	potMgr=aPtr;
	return 0;
}

int 
FactorManager::setMaxFactorSize(int size)
{
	maxFactorSize=size;
	if(maxFactorSize_Approx<maxFactorSize)
	{
		maxFactorSize_Approx=size;
	}
	return 0;
}

int 
FactorManager::setMaxFactorSize_Approx(int size)
{
	if(size>maxFactorSize)
	{
		maxFactorSize_Approx=size;
	}
	return 0;
}

int
FactorManager::setOutputDir(const char* aPtr)
{
	strcpy(outputDir,aPtr);
	return 0;
}

int
FactorManager::allocateFactorSpace()
{
	initFactorSet();
	populateFactorSet();
	return 0;
}

FactorGraph* 
FactorManager::createInitialFactorGraph()
{
	FactorGraph* fg=new FactorGraph;
	for(map<int,SlimFactor*>::iterator fIter=slimFactorSet.begin();fIter!=slimFactorSet.end();fIter++)
	{
		SlimFactor* factor=fIter->second;
		if(factor->vCnt>1)
		{
			break;
		}
		SlimFactor* newFactor=new SlimFactor;
		newFactor->vCnt=1;
		newFactor->vIds=new int[1];
		newFactor->vIds[0]=factor->vIds[0];
		newFactor->fId=factor->fId;
		newFactor->mutualInfo=factor->mutualInfo;
		newFactor->jointEntropy=factor->jointEntropy;
		newFactor->marginalEntropy=factor->jointEntropy;
		newFactor->mbScore=factor->mbScore;
		newFactor->moveScore=factor->moveScore;
		fg->setFactor(newFactor);
	}
	return fg;
}

int 
FactorManager::learnStructure()
{
	Error::ErrorCode err=estimateClusterProperties();
	if(err!=Error::SUCCESS)
	{
		cout <<Error::getErrorString(err) << endl;
		return 0;
	}
	return 0;
}

int 
FactorManager::showStructure(int f)
{
	VSET& variableSet=vMgr->getVariableSet();
	char aFName[1024];
	map<int,ofstream*> filePtrs;
	for(int k=2;k<=maxFactorSize;k++)
	{
		sprintf(aFName,"%s/strct_k%d_%d.txt",outputDir,k,f);
		ofstream* oFile=new ofstream(aFName);
		filePtrs[k]=oFile;
	}
	for(map<int,SlimFactor*>::iterator fIter=slimFactorSet.begin();fIter!=slimFactorSet.end();fIter++)
	{
		SlimFactor* aFactor=fIter->second;
		for(int k=2;k<=maxFactorSize;k++)
		{
			if(aFactor->vCnt<=k)
			{
				ofstream* oFile=filePtrs[k];
				(*oFile) <<fIter->first<< " " << factorIDToNameMap[fIter->first] 
				<<" " << aFactor->mutualInfo << " " 
				<< aFactor->jointEntropy<< endl;
			}
		}
	}
	for(int k=2;k<=maxFactorSize;k++)
	{
		ofstream* oFile=filePtrs[k];
		oFile->close();
		delete oFile;
	}
	filePtrs.clear();
	return 0;
}

int
FactorManager::readStructure(int f)
{
	char aFName[1024];
	sprintf(aFName,"%s/strct_k%d_%d.txt",outputDir,maxFactorSize,f);
	if(readClusterProperties(aFName)==-1)
	{
		return -1;
	}
	return 0;
}

//This reads a file of the format written by the above function, showStructure
int
FactorManager::readClusterProperties(const char* aFName)
{
	ifstream inFile(aFName);
	if(!inFile.good())
	{
		inFile.close();
		return -1;
	}
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		int tokCnt=0;
		int fID=0;
		string fKey;
		double jointEntropy=0;
		double mutualInfo=0;
		
		char* tok=strtok(buffer," ");
		while(tok!=NULL)
		{
			switch(tokCnt)
			{
				case 0:
				{
					fID=atoi(tok);
					break;
				}
				case 1:
				{
					fKey.append(tok);
					break;
				}
				case 2:
				{
					mutualInfo=atof(tok);
					break;
				}
				case 3:
				{
					jointEntropy=atof(tok);
					break;
				}
			}
			tok=strtok(NULL," ");
			tokCnt++;
		}
		SlimFactor* factor=slimFactorSet[fID];
		factor->mutualInfo=mutualInfo;
		factor->jointEntropy=jointEntropy;
		if(factor->vCnt==1)
		{
			factor->marginalEntropy=jointEntropy;
			factor->mbScore=jointEntropy;
			factor->moveScore=jointEntropy;
		}
	}
	inFile.close();
	return 0;
}

int 
FactorManager::readRestrictedVarlist(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	VSET& varSet=vMgr->getVariableSet();
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		string varName(buffer);
		int varID=vMgr->getVarID(varName.c_str());
		if(varID==-1)
		{
			continue;
		}
		restrictedNeighborList[varID]=varSet[varID];;
	}
	inFile.close();
	return 0;
}

SlimFactor*
FactorManager::getFactorAt(int fId)
{
	if(slimFactorSet.find(fId)==slimFactorSet.end())
	{
		return NULL;
	}
	return slimFactorSet[fId];
}

int 
FactorManager::getFactorIndex(int* vIds, int vCnt)
{
	string aKey;
	getFactorKey(vIds,vCnt,aKey);
	int fId=-1;
	if(factorNameToIDMap.find(aKey)!=factorNameToIDMap.end())
	{
		fId=factorNameToIDMap[aKey];
	}
	return fId;
}

int
FactorManager::getFactorKey(int* vIds, int vCnt, string& key)
{
	for(int v=0;v<vCnt;v++)
	{
		char aBuff[56];
		sprintf(aBuff,"-%d",vIds[v]);
		key.append(aBuff);
	}
	return 0;
}

SlimFactor*
FactorManager::getFactorFromVars(int* vIds, int vCnt)
{
	int fid=getFactorIndex(vIds,vCnt);
	if(slimFactorSet.find(fid)==slimFactorSet.end())
	{
		cout <<"Warning! Accessing null factor " << endl;
		return NULL;
	}
	return slimFactorSet[fid];
}

int
FactorManager::initFactorSet()
{
	struct timeval begintime;
	struct timeval endtime;
	gettimeofday(&begintime,NULL);
	
	int factorCnt=0;
	int fSize=1;
	int vCnt=vMgr->getVariableSet().size();
	int rCnt=restrictedNeighborList.size();
	while(fSize<=maxFactorSize)
	{
		int fCnt=0;
		if(rCnt==0)
		{
			fCnt=combCnt(vCnt,fSize);
		}
		else
		{
			fCnt=combCnt(rCnt,vCnt-rCnt,fSize);
		}
		factorCnt=factorCnt+fCnt;
		fSize++;
	}
	//Calculate the number of factors
	cout <<" Number of factors "  << factorCnt << endl;
	fSize=1;
	int fCnt=combCnt(vCnt,fSize);
	for(int i=0;i<factorCnt;i++)
	{
		SlimFactor* sFactor=new SlimFactor;
		slimFactorSet[i]=sFactor;

		if(i==fCnt)
		{
			fSize++;
			if(rCnt==0)
			{
				fCnt=fCnt+combCnt(vCnt,fSize);
			}
			else
			{
				fCnt=fCnt+combCnt(rCnt,vCnt-rCnt,fSize);
			}
		}
		sFactor->vIds=new int[fSize];
		sFactor->vCnt=fSize;
		//sFactor->secondPId=-1;
		sFactor->mutualInfo=0;
		sFactor->jointEntropy=0;
		sFactor->fId=globalFactorID;
		globalFactorID++;
	}
	cout << "Global factor id " << globalFactorID << endl;
	gettimeofday(&endtime,NULL);
	cout << "Time elapsed " << endtime.tv_sec-begintime.tv_sec<< " seconds and " << endtime.tv_usec-begintime.tv_usec << " micro secs" << endl;

	return 0;
}

//Here we simply write to the memory that we have allocated
int
FactorManager::populateFactorSet()
{
	struct timeval begintime;
	struct timeval endtime;
	gettimeofday(&begintime,NULL);
	int fSize=1;
	//These indices are simply used to create factors of size k+1
	//from factors of size k
	int startFid=0;
	int endFid=0;
	int currFid=0;
	map<int,Variable*>& variableSet=vMgr->getVariableSet();
	for(map<int,Variable*>::iterator vIter=variableSet.begin();vIter!=variableSet.end();vIter++)
	{
		SlimFactor* sFactor=slimFactorSet[currFid];
		currFid++;
		sFactor->vIds[sFactor->vCnt-1]=vIter->first;
		string fKey;
		getFactorKey(sFactor->vIds,sFactor->vCnt,fKey);
		factorNameToIDMap[fKey]=sFactor->fId;
		factorIDToNameMap[sFactor->fId]=fKey;
	}
	endFid=currFid;
	fSize++;
	
	//This int matrix is allocated before hand and reused
	//This is done to avoid the allocation and deallocation repeatedly
	//The maximum of subsets is equal to the maxFactorSize
	//Each subset can be at most maxFactorSize-1
	int** subsetSpace=new int*[maxFactorSize];
	for(int i=0;i<maxFactorSize;i++)
	{
		subsetSpace[i]=new int[maxFactorSize-1];
	}
	VSET* neighborSet=&vMgr->getVariableSet();
	if(restrictedNeighborList.size()>0)
	{
		neighborSet=&restrictedNeighborList;
	}
	while(fSize<=maxFactorSize)
	{
		int pFidCnt=endFid-startFid;
		//The number of newFids added
		//For each parent factor, iterate over the set of variables
		//and add the variable in the parent to the new factor 
		for(int i=0;i<pFidCnt;i++)
		{
			SlimFactor* pFactor=slimFactorSet[startFid+i];
			for(map<int,Variable*>::iterator vIter=neighborSet->begin();vIter!=neighborSet->end();vIter++)
			{
				int newVId=vIter->first;
				if(restrictedNeighborList.size()==0)
				{
					if(newVId<=pFactor->vIds[pFactor->vCnt-1])
					{
						continue;
					}
				}
				else
				{
					int pvar=pFactor->vIds[pFactor->vCnt-1];
					if(restrictedNeighborList.find(pvar)!=restrictedNeighborList.end())
					{
						if(newVId<=pvar)
						{
							continue;
						}
					}
					else
					{
						//Dont think I need this check
						if(pFactor->isMemberVariable(newVId))
						{
							continue;
						}
					}
				}
				SlimFactor* sFactor=slimFactorSet[currFid];
				int fIter=0;
				int dIter=0;
				while((fIter!=pFactor->vCnt) && (pFactor->vIds[fIter]<newVId))
				{
					sFactor->vIds[dIter]=pFactor->vIds[fIter];
					dIter++;
					fIter++;
				}
				sFactor->vIds[dIter]=newVId;
				dIter++;
				while(fIter<pFactor->vCnt)
				{
					sFactor->vIds[dIter]=pFactor->vIds[fIter];
					fIter++;
					dIter++;
				}
				//Here we need to add the super-set and sub-set relationships
				//Specifically, sFactor is a super-set of pFactor
				//pFactor is a sub-set of sFactor
				/*for(int j=0;j<pFactor->vCnt;j++)
				{
					sFactor->vIds[j]=pFactor->vIds[j];
				}
				sFactor->vIds[sFactor->vCnt-1]=newVId;*/
				string fKey;
				getFactorKey(sFactor->vIds,sFactor->vCnt,fKey);
				factorNameToIDMap[fKey]=sFactor->fId;
				factorIDToNameMap[sFactor->fId]=fKey;
				//sFactor->secondPId=startFid+i;
				currFid++;
				addToLattice(sFactor,subsetSpace);
			}
		}
		fSize++;
		startFid=endFid;
		endFid=currFid;
		cout <<"Populated " << endFid-startFid << " new factors "<< endl; 
	}
	for(int i=0;i<maxFactorSize;i++)
	{
		delete [] subsetSpace[i];
	}
	delete[] subsetSpace;
	cout << "Global factor id " << globalFactorID << endl;
	gettimeofday(&endtime,NULL);
	cout << "Time elapsed " << endtime.tv_sec-begintime.tv_sec<< " seconds and " << endtime.tv_usec-begintime.tv_usec << " micro secs" << endl;
	return 0;
}

//Computes n choose k, where k is very small
int
FactorManager::combCnt(int n, int k)
{
	int fCnt=1;
	int start=n;
	int i=0;
	while(i<k)
	{
		fCnt=fCnt*start;
		start--;
		i++;
	}
	start=k;
	i=1;
	while(i<=k)
	{
		fCnt=fCnt/i;
		i++;
	}
	return fCnt;
}

int
FactorManager::combCnt(int n1, int n2, int k)
{
	int fCnt=0;
	if(k==1)
	{
		fCnt=n1+n2;
	}
	else
	{
		fCnt=combCnt(n1,k)+(combCnt(n1,k-1)*n2);
	}
	return fCnt;
}

Error::ErrorCode
FactorManager::estimateClusterProperties()
{
	struct timeval begintime;
	struct timeval endtime;
	gettimeofday(&begintime,NULL);

	Error::ErrorCode err=potMgr->populatePotentialsSlimFactors(slimFactorSet,vMgr->getVariableSet());

	gettimeofday(&endtime,NULL);
	cout << "Time elapsed " << endtime.tv_sec-begintime.tv_sec<< " seconds and " << endtime.tv_usec-begintime.tv_usec << " micro secs" << endl;

	if(err!=Error::SUCCESS)
	{
		return err;
	}
	for(map<int,SlimFactor*>::iterator fIter=slimFactorSet.begin();fIter!=slimFactorSet.end();fIter++)
	{
		SlimFactor* sFactor=fIter->second;
		if(sFactor->vCnt>1)
		{
			//The multiinfo is just the difference of the joint entropy and the sum of the marginal
			//entropies of the variables
			double mi=(-1)*sFactor->jointEntropy;
			for(int j=0;j<sFactor->vCnt;j++)
			{
				SlimFactor* subFactor=getFactorFromVars(sFactor->vIds+j,1);
				mi=mi+subFactor->jointEntropy;
			}
			sFactor->mutualInfo=mi;
			if(mi<0)
			{
				cout <<"Negative mi for "<< fIter->first << "setting to 0" << endl;
				mi=0;
			}
		}
		else
		{
			sFactor->mbScore=sFactor->jointEntropy;
		}
	}
	return Error::SUCCESS;
}

//This function generates all maximal subsets of this factor
//Then it adds the subset and superset relationships
int
FactorManager::addToLattice(SlimFactor* aFactor,int** subsetSpace)
{
	aFactor->generateMaximalSubsets(subsetSpace);
	for(int i=0;i<aFactor->vCnt;i++)
	{
		int subsetId=getFactorIndex(subsetSpace[i],aFactor->vCnt-1);
		lattice.addSubset(subsetId,aFactor->fId);
		//lattice.addSuperset(aFactor->fId,subsetId);
	}
	return 0;
}
