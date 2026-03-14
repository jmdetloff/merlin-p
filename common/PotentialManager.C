#include <iostream>
#include <cstring>
#include <math.h>

#include "CommonTypes.H"
#include "Error.H"
#include "Variable.H"
#include "Potential.H"
#include "Evidence.H"
#include "EvidenceManager.H"
#include "SlimFactor.H"
#include "PotentialManager.H"

PotentialManager::PotentialManager()
{
	ludecomp=NULL;
	perm=NULL;
	randomData=false;
}

PotentialManager::~PotentialManager()
{
	if(ludecomp!=NULL)
	{
		gsl_matrix_free(ludecomp);
	}
	if(perm!=NULL)
	{
		gsl_permutation_free(perm);
	}
}

int 
PotentialManager::setEvidenceManager(EvidenceManager* aPtr)
{
	evMgr=aPtr;
	return 0;
}

int
PotentialManager::setOutputDir(const char* aDirName)
{
	strcpy(outputDir,aDirName);
	return 0;
}

int
PotentialManager::setRandom(bool flag)
{
	randomData=flag;
	return 0;
}

int
PotentialManager::init(int f)
{
	char mFName[1024];
	char sdFName[1024];
	sprintf(mFName,"%s/gauss_mean_%d.txt",outputDir,f);
	sprintf(sdFName,"%s/gauss_std_%d.txt",outputDir,f);
	ifstream inFile(mFName);
	if(inFile.good())
	{
		readAllMeanCov(mFName,sdFName);
	}
	else
	{
		INTINTMAP& trainEvidSet=evMgr->getTrainingSet();
		estimateAllMeanCov(randomData,globalMean,globalCovar,trainEvidSet);
	}
	ludecomp=gsl_matrix_alloc(MAXFACTORSIZE_ALLOC,MAXFACTORSIZE_ALLOC);
	perm=gsl_permutation_alloc(MAXFACTORSIZE_ALLOC);
	return 0;
}

int
PotentialManager::reset()
{
	globalMean.clear();
	for(map<int,INTDBLMAP*>::iterator cIter=globalCovar.begin();cIter!=globalCovar.end();cIter++)
	{
		cIter->second->clear();
		delete cIter->second;
	}
	globalCovar.clear();
	if(ludecomp!=NULL)
	{
		gsl_matrix_free(ludecomp);
		ludecomp=NULL;
	}
	if(perm!=NULL)
	{
		gsl_permutation_free(perm);
		perm=NULL;
	}
	return 0;
}

int
PotentialManager::estimateAllMeanCov(bool random, INTDBLMAP& gMean, map<int,INTDBLMAP*>& gCovar,INTINTMAP& trainEvidSet)
{
	int evidCnt=trainEvidSet.size();
	//First get the mean and then the variance
	int dId=0;
	for(INTINTMAP_ITER eIter=trainEvidSet.begin();eIter!=trainEvidSet.end();eIter++)
	{
		EMAP* evidMap=NULL;
		if(random)
		{
			evidMap=evMgr->getRandomEvidenceAt(eIter->first);
		}
		else
		{
			evidMap=evMgr->getEvidenceAt(eIter->first);
		}
		for(EMAP_ITER vIter=evidMap->begin();vIter!=evidMap->end(); vIter++)
		{
			int vId=vIter->first;
			Evidence* evid=vIter->second;
			double val=evid->getEvidVal();
			if(gMean.find(vId)==gMean.end())
			{
				gMean[vId]=val;
			}
			else
			{
				gMean[vId]=gMean[vId]+val;
			}
		}
		dId++;	
	}
	//Now estimate the mean
	for(INTDBLMAP_ITER idIter=gMean.begin();idIter!=gMean.end();idIter++)
	{
		idIter->second=idIter->second/(double) evidCnt;
		INTDBLMAP* vcov=new INTDBLMAP;
		gCovar[idIter->first]=vcov;
	}
	return 0;
	int covPair=0;
	//Now the variance
	for(INTINTMAP_ITER eIter=trainEvidSet.begin();eIter!=trainEvidSet.end();eIter++)
	{
		EMAP* evidMap=NULL;
		if(random)
		{
			evidMap=evMgr->getRandomEvidenceAt(eIter->first);
		}
		else
		{
			evidMap=evMgr->getEvidenceAt(eIter->first);
		}
		for(EMAP_ITER vIter=evidMap->begin();vIter!=evidMap->end(); vIter++)
		{
			int vId=vIter->first;
			Evidence* evid=vIter->second;
			double vval=evid->getEvidVal();
			double vmean=gMean[vId];
			INTDBLMAP* vcov=NULL;
			if(gCovar.find(vId)==gCovar.end())
			{
				vcov=new INTDBLMAP;
				gCovar[vId]=vcov;
			}
			else
			{
				vcov=gCovar[vId];
			}
			for(EMAP_ITER uIter=vIter;uIter!=evidMap->end();uIter++)
			{
				int uId=uIter->first;
				Evidence* evid1=uIter->second;
				double uval=evid1->getEvidVal();
				double umean=gMean[uId];
				double diffprod=(vval-vmean)*(uval-umean);
				INTDBLMAP* ucov=NULL;
				if(gCovar.find(uId)==gCovar.end())
				{
					ucov=new INTDBLMAP;
					gCovar[uId]=ucov;
				}
				else
				{
					ucov=gCovar[uId];
				}
				if(vcov->find(uId)==vcov->end())
				{
					covPair++;
					(*vcov)[uId]=diffprod;
				}
				else
				{
					(*vcov)[uId]=(*vcov)[uId]+diffprod;
				}
				if(uId!=vId)
				{
					if(ucov->find(vId)==ucov->end())
					{
						(*ucov)[vId]=diffprod;
					}
					else
					{
						(*ucov)[vId]=(*ucov)[vId]+diffprod;
					}
				}
			}
		}

	}
	//Now estimate the variance
	for(map<int,INTDBLMAP*>::iterator idIter=gCovar.begin();idIter!=gCovar.end();idIter++)
	{
		INTDBLMAP* var=idIter->second;
		for(INTDBLMAP_ITER vIter=var->begin();vIter!=var->end();vIter++)
		{
			if(vIter->first==idIter->first)
			{
				vIter->second=(0.001+vIter->second)/((double)(evidCnt-1));
			}
			else
			{
				vIter->second=vIter->second/((double)(evidCnt-1));
			}
		}
	}
	return 0;
}

int
PotentialManager::estimateCovariance(bool random,INTDBLMAP* vcov, int uId, int vId)
{
	INTINTMAP& trainEvidSet=evMgr->getTrainingSet();
	int evidCnt=trainEvidSet.size();
	INTDBLMAP* ucov=globalCovar[uId];
	for(INTINTMAP_ITER eIter=trainEvidSet.begin();eIter!=trainEvidSet.end();eIter++)
	{
		EMAP* evidMap=NULL;
		if(random)
		{
			evidMap=evMgr->getRandomEvidenceAt(eIter->first);
		}
		else
		{
			evidMap=evMgr->getEvidenceAt(eIter->first);
		}
		Evidence* evid=(*evidMap)[vId];
		double vval=evid->getEvidVal();
		double vmean=globalMean[vId];
		Evidence* evid1=(*evidMap)[uId];
		double uval=evid1->getEvidVal();
		double umean=globalMean[uId];
		double diffprod=(vval-vmean)*(uval-umean);
		if(vcov->find(uId)==vcov->end())
		{
			(*vcov)[uId]=diffprod;
		}
		else
		{
			(*vcov)[uId]=(*vcov)[uId]+diffprod;
		}
		if(uId!=vId)
		{
			if(ucov->find(vId)==ucov->end())
			{
				(*ucov)[vId]=diffprod;
			}
			else
			{
				(*ucov)[vId]=(*ucov)[vId]+diffprod;
			}
		}
	}
	//Now estimate the variance
	if(uId==vId)
	{
		double ssd=(*ucov)[uId];
		(*ucov)[uId]=(0.001+ssd)/((double)(evidCnt-1));
	}
	else
	{
		double ssduv=(*ucov)[vId];
		(*ucov)[vId]=ssduv/((double)(evidCnt-1));
		(*vcov)[uId]=ssduv/((double)(evidCnt-1));
	}
	return 0;
}

int
PotentialManager::readAllMeanCov(const char* mFName, const char* sdFName)
{
	ifstream mFile(mFName);
	ifstream sdFile(sdFName);
	char buffer[1024];
	while(mFile.good())
	{
		mFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		int vId;
		double mean=0;
		while(tok!=NULL)
		{	
			if(tokCnt==0)
			{
				vId=atoi(tok);	
			}
			else if(tokCnt==1)
			{
				mean=atof(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		globalMean[vId]=mean;
	}
	mFile.close();
	int lineNo=0;
	while(sdFile.good())
	{
		sdFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		int vId=0;
		int uId=0;
		double covariance=0;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				uId=atoi(tok);
			}
			else if(tokCnt==1)
			{
				vId=atoi(tok);
			}
			else if(tokCnt==2)
			{
				covariance=atof(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		INTDBLMAP* ucov=NULL;
		if(globalCovar.find(uId)==globalCovar.end())
		{
			ucov=new INTDBLMAP;
			globalCovar[uId]=ucov;
		}
		else
		{
			ucov=globalCovar[uId];
		}
		(*ucov)[vId]=covariance;
		if(uId==0 && vId==0)
		{
			cout << "Found uId=0 vId=0 covar="<< covariance << " at lineno " << lineNo  << endl;
		}
		lineNo++;
	}
	sdFile.close();
	return 0;
}

Error::ErrorCode
PotentialManager::populatePotentialsSlimFactors(map<int,SlimFactor*>& factorSet,VSET& varSet)
{
	//The set of flags to keep status of the potentials that have been calculated
	map<int,bool> doneFlag;
	for(map<int,SlimFactor*>::iterator fIter=factorSet.begin();fIter!=factorSet.end();fIter++)
	{
		doneFlag[fIter->first]=false;
	}
	int popFId=0;
	for(map<int,SlimFactor*>::reverse_iterator rIter=factorSet.rbegin();rIter!=factorSet.rend();rIter++)
	{
		//If we have computed the potential for this flag move one
		if(doneFlag[rIter->first])
		{
			popFId++;
			continue;
		}
		SlimFactor* sFactor=rIter->second;
		if(sFactor->fId==176)
		{
			cout <<"Stop here " << endl;
		}
		//Otherwise create the potential
		Potential* aPotFunc=new Potential;
		for(int j=0;j<sFactor->vCnt;j++)
		{
			Variable* aVar=varSet[sFactor->vIds[j]];
			if(j==sFactor->vCnt-1)
			{
				aPotFunc->setAssocVariable(aVar,Potential::FACTOR);
			}
			else
			{
				aPotFunc->setAssocVariable(aVar,Potential::MARKOV_BNKT);
			}
		}
		aPotFunc->potZeroInit();
		populatePotential(aPotFunc,false);
		aPotFunc->calculateJointEntropy();
		sFactor->jointEntropy=aPotFunc->getJointEntropy();
		if(sFactor->jointEntropy<0)
		{
		//	sFactor->jointEntropy=0;
		//	cout <<"Negative entropy for " << sFactor->fId << endl;
		}
		doneFlag[rIter->first]=true;
		delete aPotFunc;
		if(popFId%100000==0)
		{
			cout <<"Done with " << factorSet.size()-popFId << " factors " << endl;
		}
		popFId++;
	}
	return Error::SUCCESS;
}

int
PotentialManager::populatePotential(Potential* aPot, bool random)
{
	VSET& potVars=aPot->getAssocVariables();
	for(VSET_ITER vIter=potVars.begin();vIter!=potVars.end(); vIter++)
	{
		if(globalMean.find(vIter->first)==globalMean.end())
		{
			cerr <<"No var with id " << vIter->first << endl;
			exit(-1);
		}

		double mean=globalMean[vIter->first];
		INTDBLMAP* covar=globalCovar[vIter->first];
		aPot->updateMean(vIter->first,mean);
		
		for(VSET_ITER uIter=vIter;uIter!=potVars.end();uIter++)
		{
			if(covar->find(uIter->first)==covar->end())
			{
				estimateCovariance(random,covar,uIter->first,vIter->first);
			}
			double cval=(*covar)[uIter->first];
			aPot->updateCovariance(vIter->first,uIter->first,cval);
			aPot->updateCovariance(uIter->first,vIter->first,cval);
		}
	}

	aPot->makeValidJPD(ludecomp,perm);
	return 0;
}
