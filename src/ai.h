#ifndef __AI_H__
#define __AI_H__


#include <vector>
#include <numeric>
#include <tuple>
#include <unordered_map>
#include <iostream>
#include "mapping.h"
#include "evaluation.h"



class Particle {
public:
	Particle(std::vector<double> inpPos,
	         std::vector<std::tuple<double, double>> inpBounds,
	         double w = 0.7,     // w:   constant inertia weight (how much to weigh the previous velocity)
	         double c1 = 2,      // c1:  cognative constant
	         double c2 = 1);     // c2:  social constant

	void VelocityUpdate(Particle& globalBestParticle);
	void PositionUpdate();
	double Score();
	virtual void Evaluate() = 0;
	int dimensions;

protected:
	std::vector<double> vecPosition, vecVelocity, vecBestPos;
	std::vector<std::tuple<double, double>> vecBounds;
	double w_vel, c_cog, c_social,
	       currScore, bestScore;
};


class ProteinParticle : public Particle {
public:
	ProteinParticle(std::string protein_sequence);
	std::string Sequence();
	void Evaluate();
private:
	std::vector<double> Vectorize(std::string protein_sequence);
};



class AI {
public:
	std::string Initial_State(int length = PROT_LEN);
	std::string Stochastic_Hill_Climbing(std::string curr_state, int max_try = 500);
	std::string Simulated_Annealing(std::string curr_state, double T, double coef, int max_try = 500);
	ProteinParticle ParticleSwarmOptimization(int num_particles, int max_iterations);
	std::string AntColonyOptimization(int NUMBEROFANTS, int ITERATIONS,
                                  double ALPHA, double BETA, double Q, double RHO);

private:
	std::string RandomReplace(std::string protein_sequence, int maxChanges);
	std::string FragmentReplace(std::string protein_sequence);
	std::string FillGaps(std::string protein_sequence);
	void MinimumSegmentDensity(int N, int L, std::vector<double> &scores, int &ansL, int &ansR);
	void TetraScore(std::string target_sequence, std::vector<double> &vecScore);
	bool FragmentExists(int i, int j);
	std::string FragmentFetch(int i, int j);
};


#endif // __AI_H__
