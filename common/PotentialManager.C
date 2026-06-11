#include <iostream>
#include <cstring>
#include <math.h>

#include "CommonTypes.H"
#include "Error.H"
#include "Potential.H"
#include "Evidence.H"
#include "EvidenceManager.H"
#include "PotentialManager.H"

PotentialManager::PotentialManager()
{
	ludecomp=NULL;
	perm=NULL;
	globalCovariances=nullptr;
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
	if(globalCovariances!=nullptr)
	{
		delete globalCovariances;
	}
}

int
PotentialManager::init(EvidenceManager* evMgr, bool randomData, vector<int>& varIDs)
{
	if (globalCovariances != nullptr)
	{
		delete globalCovariances;
	}

	INTINTMAP& trainEvidSet = evMgr->getTrainingSet();
	EMAP* evidMap = evMgr->getEvidenceAt(trainEvidSet.begin()->first);
	int varCount = evidMap->size();
	int sampleCount = trainEvidSet.size();

	globalMeans.clear();

	globalCovariances = new Matrix(varCount, varCount);
	globalCovariances->setAllValues(-1);

	// Stores the deviations from the mean for each variable and sample.
	vector<double> deviations(varCount * sampleCount, 0);

	// Copy all the samples into the data matrix
	int sampleIndex = 0;
	for (INTINTMAP_ITER eIter = trainEvidSet.begin(); eIter != trainEvidSet.end(); eIter++)
	{
		EMAP* evidMap=NULL;
		if(randomData)
		{
			evidMap=evMgr->getRandomEvidenceAt(eIter->first);
		}
		else
		{
			evidMap=evMgr->getEvidenceAt(eIter->first);
		}
		for (EMAP_ITER vIter = evidMap->begin(); vIter != evidMap->end(); vIter++)
		{
			int vId = vIter->first;
			Evidence* evid = vIter->second;
			double val = evid->getEvidVal();
			deviations[vId * sampleCount + sampleIndex] = val;
		}
		sampleIndex++;
	}

	// Done copying. Now we can go over data and get the means
	for (int i = 0; i < varCount; i++)
	{
		double sampleSum = 0;
		for(int j = 0; j < sampleCount; j++)
		{
			sampleSum += deviations[i * sampleCount + j];
		}
		globalMeans.push_back(sampleSum / sampleCount);
	}

	// Finally, use the means to pre-center the data
	for (int i = 0; i < trainEvidSet.size(); i++)
	{
		for (int j = 0; j < varCount; j++)
		{
			deviations[j * sampleCount + i] -= globalMeans[j];
		}
	}

	int norm = sampleCount - 1;

	// Set covariances along the diagonal.
	for (int i = 0; i < varCount; i++)
	{
		double ssd = 0.001;
		for (int j = 0; j < sampleCount; j++)
		{
			double dev = deviations[i * sampleCount + j];
			ssd += dev * dev;
		}
		globalCovariances->setValue(ssd / norm, i, i);
	}

	// Set covariances between regulators and all other variables.
	for (int i = 0; i < varIDs.size(); i++)
	{
		int regID = varIDs[i];
		for (int j = 0; j < varCount; j++)
		{
			if (regID == j)
			{
				continue;
			}

			double ssd = 0;
			for (int k = 0; k < sampleCount; k++)
			{
				double devI = deviations[regID * sampleCount + k];
				double devJ = deviations[j * sampleCount + k];
				ssd += devI * devJ;
			}

			double covariance = ssd / norm;
			globalCovariances->setValue(covariance, regID, j);
			globalCovariances->setValue(covariance, j, regID);
		}
	}

	ludecomp = gsl_matrix_alloc(MAXFACTORSIZE_ALLOC, MAXFACTORSIZE_ALLOC);
	perm = gsl_permutation_alloc(MAXFACTORSIZE_ALLOC);

	return 0;
}

Potential*
PotentialManager::createPotential(int factorID)
{
	int varCount = globalMeans.size();
	double variance = globalCovariances->getValue(factorID, factorID);
	double bias = globalMeans[factorID];
	INTDBLMAP weights;
	return new Potential(factorID, variance, bias, weights);
}

void
PotentialManager::computeLLs(int factorID, int sampleSize, vector<int>&candidateParents, unordered_map<int, double>&scores) {

	double variance = globalCovariances->getValue(factorID, factorID);

	for (int i = 0; i < candidateParents.size(); i++) {
		int candidateID = candidateParents[i];

		// Var(C)
		double candidateVariance = globalCovariances->getValue(candidateID, candidateID);

		double factorCandidateCov = globalCovariances->getValue(factorID, candidateID);

		double finalVariance = variance - factorCandidateCov * factorCandidateCov / candidateVariance;

		if (finalVariance < 1e-5) {
			finalVariance = 1e-5;
		}

		if(isnan(finalVariance) || isinf(finalVariance)) {
			continue;
		}
		
		scores[candidateID] = -0.5 * ((sampleSize - 1) + sampleSize * log(2 * PI) + sampleSize * log(finalVariance));
	}
}

void
PotentialManager::computeLLs(int factorID, int sampleSize, vector<int>& existingParents, vector<int>&candidateParents, unordered_map<int, double>&scores) {

	if (existingParents.size() == 0) {
		computeLLs(factorID, sampleSize, candidateParents, scores);
		return;
	}

	double variance = globalCovariances->getValue(factorID, factorID);

	int parentCount = existingParents.size();

	// A : factor
	// B : existing parents
	// C : candidate parent
	// Create Cov(AB) and Cov(BB)

	gsl_matrix* existingParentCovariances = gsl_matrix_alloc(parentCount, parentCount);
	gsl_vector* existingParentMarginalVariances = gsl_vector_alloc(parentCount);

	for (int i = 0; i < parentCount; i++) {
		int varAID = existingParents[i];
		double marginalCovariance = globalCovariances->getValue(factorID, varAID);
		gsl_vector_set(existingParentMarginalVariances, i, marginalCovariance);

		for (int j = i; j < parentCount; j++) {
			int varBID = existingParents[j];
			double covariance = globalCovariances->getValue(varAID, varBID);
			gsl_matrix_set(existingParentCovariances, i, j, covariance);
			gsl_matrix_set(existingParentCovariances, j, i, covariance);
		}
	}

	// Create Cov(BB)^-1

	gsl_permutation* permutation = gsl_permutation_alloc(parentCount);

	int signum=0;
	gsl_linalg_LU_decomp(existingParentCovariances, permutation, &signum);

	gsl_matrix* parentCovInverse = gsl_matrix_alloc(parentCount, parentCount);

	gsl_linalg_LU_invert(existingParentCovariances, permutation, parentCovInverse);

	gsl_matrix_free(existingParentCovariances);
	gsl_permutation_free(permutation);

	// Create Cov(AB) * Cov(BB)^-1

	gsl_vector* prod = gsl_vector_alloc(parentCount);
	gsl_blas_dgemv(CblasTrans, 1, parentCovInverse, existingParentMarginalVariances, 0, prod);

	// Make variance hold Var(A|B) = Var(A) - Cov(AB)Var(BB)^-1 Cov(AB)

	for (int i = 0; i < parentCount; i++) {
		int vID = existingParents[i];
		double aVal = gsl_vector_get(prod, i);
		double bVal = gsl_vector_get(existingParentMarginalVariances, i);
		variance -= aVal * bVal;
	}

	gsl_vector_free(existingParentMarginalVariances);

	for (int i = 0; i < candidateParents.size(); i++) {
		int candidateID = candidateParents[i];

		// Var(C)
		double candidateVariance = globalCovariances->getValue(candidateID, candidateID);

		// Cov(BC)
		gsl_vector* candidateMarginalVariances = gsl_vector_alloc(parentCount);

		for (int i = 0; i < parentCount; i++) {
			int varAID = existingParents[i];
			double marginalCovariance = globalCovariances->getValue(candidateID, varAID);
			gsl_vector_set(candidateMarginalVariances, i, marginalCovariance);
		}

		// Create Cov(BC) * Cov(BB)^-1
		gsl_vector* candidateProd = gsl_vector_alloc(parentCount);
		gsl_blas_dgemv(CblasTrans, 1, parentCovInverse, candidateMarginalVariances, 0, candidateProd);

		// Calc Var(C|B) = Var(BC) - Cov(CB)Var(BB)^-1 Cov(BC)

		for (int i = 0; i < parentCount; i++) {
			int vID = existingParents[i];
			double aVal = gsl_vector_get(candidateProd, i);
			double bVal = gsl_vector_get(candidateMarginalVariances, i);
			candidateVariance -= aVal * bVal;
		}

		// Cov(AC)
		double candidateFactorCovariance = globalCovariances->getValue(candidateID, factorID);

		// Cov(AB) * Var(BB)^-1 * Cov(BC)
		double dot = 0.0;
		gsl_blas_ddot(prod, candidateMarginalVariances, &dot);

		gsl_vector_free(candidateProd);
		gsl_vector_free(candidateMarginalVariances);

		double factorAndCandidateVarianceConditionedOnParents = candidateFactorCovariance - dot;

		double finalVariance = variance - factorAndCandidateVarianceConditionedOnParents * factorAndCandidateVarianceConditionedOnParents / candidateVariance;
		
		if(finalVariance < 1e-5) {
			finalVariance = 1e-5;
		}

		if(isnan(finalVariance) || isinf(finalVariance)) {
			continue;
		}
		
		scores[candidateID] = -0.5 * ((sampleSize - 1) + sampleSize * log(2 * PI) + sampleSize * log(finalVariance));
	}

	gsl_matrix_free(parentCovInverse);
	gsl_vector_free(prod);
}

double
PotentialManager::computeLL(int factorID, vector<int>& parentIDs, int sampleSize, Potential** newPot)
{
	double variance = globalCovariances->getValue(factorID, factorID);
	double bias = globalMeans[factorID];
	INTDBLMAP weights;

	int parentCount = parentIDs.size();

	// Start by collecting a matrix of all the covariances of the conditioning variables,
	// and the marginal variances of the conditioning variables.

	gsl_matrix* parentCovariances = gsl_matrix_alloc(parentCount, parentCount);
	gsl_vector* parentMarginalVariances = gsl_vector_alloc(parentCount);

	for (int i = 0; i < parentCount; i++) {
		int varAID = parentIDs[i];
		double marginalCovariance = globalCovariances->getValue(factorID, varAID);
		gsl_vector_set(parentMarginalVariances, i, marginalCovariance);

		for (int j = i; j < parentCount; j++) {
			int varBID = parentIDs[j];
			double covariance = globalCovariances->getValue(varAID, varBID);
			gsl_matrix_set(parentCovariances, i, j, covariance);
			gsl_matrix_set(parentCovariances, j, i, covariance);
		}
	}

	// Compute the final values for the variance of the conditional gaussian,
	// plus the regression parameters for the mean of the conditional guassian.

	gsl_vector* x = gsl_vector_alloc(parentCount);

	gsl_permutation* permutation = gsl_permutation_alloc(parentCount);

	int signum = 0;
	gsl_linalg_LU_decomp(parentCovariances, permutation, &signum);

	gsl_linalg_LU_solve(parentCovariances, permutation, parentMarginalVariances, x);

	for (int i = 0; i < parentCount; i++) {
		int vID = parentIDs[i];
		double aVal = gsl_vector_get(x, i);
		double bVal = gsl_vector_get(parentMarginalVariances, i);
		double cVal = globalMeans[vID];
		weights[vID] = aVal;
		variance -= aVal * bVal;
		bias -= cVal * aVal;
	}

	gsl_vector_free(x);
	gsl_vector_free(parentMarginalVariances);
	gsl_permutation_free(permutation);
	gsl_matrix_free(parentCovariances);

	if(variance < 1e-5) {
		variance = 1e-5;
	}

	// If the variance is invalid, then we don't want to attempt adding this edge,
	// so we should just bail out before computing the LL
	if(isnan(variance) || isinf(variance)) {
		return -1;
	}

	// Now that the conditional Gaussian params are computed, we can create the potential.
	*newPot = new Potential(factorID, variance, bias, weights);

	// Finally, compute the conditional log likelihood.
	return -0.5 * ((sampleSize - 1) + sampleSize * log(2 * PI) + sampleSize * log(variance));
}
