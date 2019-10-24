/*

THIS FILE IS NOT COMPATIBLE WITH THE NEWER VERSION

*/




#include "ai.h"



std::string AI::Stochastic_Hill_Climbing(std::string curr_state, int max_try) {
    double score = O_Fitna(curr_state);

    for (int i = 0; i < max_try; ++i) {
        std::string new_state;

        if (UniformRand(0, 1) < 0.6)
            new_state = FragmentReplace(curr_state);
        else
            new_state = RandomReplace(curr_state, 7);

        double new_score = O_Fitna(new_state);
        if (new_score > score) {
            score = new_score;
            curr_state = new_state;
            i = 0;
        }
    }
    return curr_state;
}


std::string AI::Simulated_Annealing(std::string curr_state, double T, double coef, int max_try) {
    double score = O_Fitna(curr_state);
    for (int i = 0; i < max_try; ++i) {
        std::string new_state;

        if (UniformRand(0, 1) < 0.6)
            new_state = FragmentReplace(curr_state);
        else
            new_state = RandomReplace(curr_state, 7);

        double new_score = O_Fitna(new_state);
        if (new_score >= score ||
                exp((new_score - score) / T) >= UniformRand(0, 1)) {

            score = new_score;
            curr_state = new_state;
            i = 0;
        }
        T *= coef;
    }
    return curr_state;
}


void AI::TetraScore(std::string protein_sequence, std::vector<double> &vecScore) {
    int n = protein_sequence.size();
    vecScore.clear(); vecScore.resize(n);
    std::vector<int> atom_id(n);
    double sc[4];

    for (int i = 0; i < n; ++i) {
        std::string key = std::string() + protein_sequence[i] + " CA";
        atom_id[i] = atom_map[key];
    }

    for (int i = 0; i < n; ++i) { assert(vecScore[i] == 0); assert(symbol_map.count(protein_sequence[i]) > 0);}

    for (int i = 0; i < n; ++i) {
        char sym = protein_sequence[i];
        int idx = symbol_map[sym];
        sc[0] = burial_potential[idx][target_solvent_accessibility[i]];
        sc[1] = structure_potential[idx][target_secondary_structure[i]];
        sc[2] = alignment_potential[i][AA_classification[sym]];
        sc[3] = 0;
        for (int j = 0; j < n; ++j) {
            int b = target_atom_distance[i][j] * 2;
            if (b >= num_bins || b == 0)
                continue;
            sc[3] += edDFIRE[atom_id[i]][atom_id[j]][b];
        }
        sc[3] /= dDFIRE_COEF;

        for (int j = 0; j < 4; ++j) {
            vecScore[i] += eNormalaize(sc[j], sc_weight[j][0], sc_weight[j][1]) *
                           sc_weight[j][2];
        }
    }
}


std::string AI::Initial_State(int length) {
    std::string seq('-', length);
    int pos = 0, lb = 10, ub = 30;
    while (pos != length) {
        int i = pos - rand() % std::min(pos + 1, lb),
            j = pos + rand() % std::min(length - pos, ub);

        if (FragmentExists(i, j) == false)
            continue;
        seq.replace(i, j - i + 1, FragmentFetch(i, j));
        pos = j + 1;
    }

    return FillGaps(seq);
}


std::string AI::FillGaps(std::string protein_sequence) {
    for (auto &c : protein_sequence)
        if (c == '-' || c == 'X')
            c = index_to_symbol[rand() % 20];
    return protein_sequence;
}


std::string AI::RandomReplace(std::string protein_sequence, int maxChanges) {
    int cnt = rand() % (maxChanges - 1) + 1;
    for (int i = 0; i < cnt; ++i) {
        int j = rand() % protein_sequence.length();
        protein_sequence[j] = index_to_symbol[rand() % 20];
    }
    return protein_sequence;
}


std::string AI::FragmentReplace(std::string protein_sequence) {
    std::vector<double> vecScore;
    int n = protein_sequence.length(),
        l, r;

    TetraScore(protein_sequence, vecScore);
    MinimumSegmentDensity(n, FRAG_MIN_LEN * UniformRand(1, 4), vecScore, l, r);
    double sss = 0;
    for (int i = l; i <= r; ++i)
        sss += vecScore[i];
    // std::cerr << r - l + 1 << "    " << l << " " << r << "\t" << sss / (r-l+1)<< "\n";

    while (!FragmentExists(l, r))
        l++;
    protein_sequence.replace(l, r - l + 1, FragmentFetch(l, r));
    return FillGaps(protein_sequence);
}

bool AI::FragmentExists(int i, int j) {
    return !matAMR[i][j].empty();
}


std::string AI::FragmentFetch(int i, int j) {
    return matAMR[i][j][rand() % matAMR[i][j].size()];
}


void AI::MinimumSegmentDensity(int N, int L, std::vector<double> &scores, int &ansL, int &ansR) {
    std::vector<double> P(N + 1), deque(N + 1);
    double ans = DBL_MAX,
           *qf, *qb;

    partial_sum(scores.begin(), scores.end(), P.begin() + 1);
    auto slope = [&P](int x, int y) { return (P[y + 1] - P[x]) / (y + 1 - x); };

    ansL = 0;
    ansR = L - 1;
    qf = qb = deque.data();

    for (int j = L - 1; j < N; ++j)
    {
        while (qb - qf >= 2 && slope(qb[-2], qb[-1] - 1) <= slope(qb[-2], j - L)) {
            qb--;
        }
        *qb++ = j - L + 1;

        while (qb - qf >= 2 && slope(qf[0], qf[1] - 1) >= slope(qf[0], j)) {
            qf++;
        }
        if (slope(qf[0], j) < ans) {
            ans = slope(qf[0], j);
            ansL = qf[0];
            ansR = j;
        }
    }
}


ProteinParticle AI::ParticleSwarmOptimization(int num_particles, int max_iterations) {
    std::vector<ProteinParticle> Swarm;
    for (int i = 0; i < num_particles; ++i)
        Swarm.emplace_back(ProteinParticle(Initial_State()));
    // std::cerr << "1\n";
    ProteinParticle globalBestParticle("");
    // std::cerr << "2\n";

    for (int k = 0; k < max_iterations; ++k) {
        for (int i = 0; i < num_particles; ++i) {
            Swarm[i].Evaluate();
            // std::cerr << "3\n";
            // std::cerr << Swarm[i].Score() << "\n";
            if (Swarm[i].Score() > globalBestParticle.Score()) {
                globalBestParticle = Swarm[i];
            }
        }
        std::cerr << k << " " << globalBestParticle.Score() << "\t" << globalBestParticle.Sequence() << " " << globalBestParticle.dimensions << "\n";

        for (int i = 0; i < num_particles; ++i) {
            Swarm[i].VelocityUpdate(globalBestParticle);
            // std::cerr << "5\n";

            Swarm[i].PositionUpdate();
        }
    }
    return globalBestParticle;
}



Particle::Particle(std::vector<double> inpPos,
                   std::vector<std::tuple<double, double>> inpBounds,
                   double w,       // w:   constant inertia weight (how much to weigh the previous velocity)
                   double c1,      // c1:  cognative constant
                   double c2)      // c2:  social constant
    : vecPosition( inpPos ), vecBounds( inpBounds ), dimensions( inpPos.size() ),
      w_vel ( w ), c_cog( c1 ), c_social ( c2 ), currScore( 0 ), bestScore( 0 )
{
    assert(inpPos.size() == inpBounds.size());
    vecBestPos.resize(dimensions);
    vecVelocity.resize(dimensions);
}


void Particle::PositionUpdate() {
    double lB, uB;
    for (int i = 0; i < dimensions; ++i) {
        std::tie(lB, uB) = vecBounds[i];
        vecPosition[i] = std::max(lB, std::min(uB, vecPosition[i] + vecVelocity[i]));
    }
}


void Particle::VelocityUpdate(Particle& globalBestParticle) {
    // std::cerr << "H1\n";
    for (int i = 0; i < dimensions; ++i) {
        double r1 = UniformRand(0, 1),
               r2 = UniformRand(0, 1);
        // std::cerr << "H2\n";

        double vel_cognitive = c_cog * r1 * (vecBestPos[i] - vecPosition[i]);
        // std::cerr << "H3\n";

        double vel_social = c_social * r2 * (globalBestParticle.vecPosition[i] - vecPosition[i]);
        // std::cerr << "H4\n";

        vecVelocity[i] = w_vel * vecVelocity[i] + vel_cognitive + vel_social;
    }
}


double Particle::Score() { return currScore; }



ProteinParticle::ProteinParticle(std::string protein_sequence)
    : Particle(Vectorize(protein_sequence),
               std::vector<std::tuple<double, double>> (protein_sequence.length(), std::make_tuple(0, AA_TYPE - 1))
              )
{
    for (auto& v : vecVelocity)
        v = UniformRand(0, 20);
}


std::string ProteinParticle::Sequence() {
    std::string sequence;
    for (auto p : vecPosition)
        sequence += index_to_symbol[std::lrint(p)];
    return sequence;
}


void ProteinParticle::Evaluate() {
    currScore = O_Fitna(Sequence());
    if (currScore > bestScore) {
        bestScore = currScore;
        vecBestPos = vecPosition;
    }
}


std::vector<double> ProteinParticle::Vectorize(std::string sequence) {
    std::vector<double> vPosition;
    for (auto c : sequence)
        vPosition.push_back(symbol_map[c]);
    return vPosition;
}


std::string AI::AntColonyOptimization(int NUMBEROFANTS, int ITERATIONS,
                                  double ALPHA, double BETA, double Q, double RHO) {
    
    std::unordered_map<std::string, double> PHEROMONES; // phermones of states
    std::string bestSeq;
    double bestFit = 0;
    double MAXSC = 0;
    for (int i = 0; i < 6; ++i)
        MAXSC += sc_weight[i][2];


    auto heuristic = [&](std::string seq) {  // heuristic
        std::vector<double> vc;
        TetraScore(seq, vc);
        return (double)std::accumulate(vc.begin(), vc.end(), 0) / seq.length();
    };

    for (int iterations = 1; iterations <= ITERATIONS; ++iterations) {
        std::cout << "ITERATION " << iterations << ":\n";

        for (int k = 0; k < NUMBEROFANTS; ++k) {
            std::string currentSeq;

            // find a path
            for (int idx = 0; idx < PROT_LEN; ++idx) {
                // calcualte probabilities
                std::array<double, 20> prob;
                for (int aa = 0; aa < 20; ++aa) {
                    std::string state = currentSeq + index_to_symbol[aa]; // next state
                    double ETA = std::pow(heuristic(state), BETA);
                    double TAU = std::pow(PHEROMONES[state] + 1e-2, ALPHA);
                    // std::cerr << ETA << " " << TAU << "\n";
                    prob[aa] = ETA * TAU;
                    prob[aa] += (aa != 0 ? prob[aa - 1] : 0); // sum with previous (dp)
                }

                // pick a node
                double randProb = UniformRand(0, 1) * prob[19]; // rand between 0 and SUM (=prob[19])
                int choice = 0;
                while (randProb >= prob[choice] && choice < 19)
                    choice++;

                currentSeq += index_to_symbol[choice];
            }
            double fitness = O_Fitna(currentSeq);

            // update pheromones
            for (int idx = 0; idx < PROT_LEN; ++idx) {
                std::string state = currentSeq.substr(0, idx + 1);
                PHEROMONES[state] = (1 - RHO) * PHEROMONES[state] + (Q / (MAXSC - fitness));
            }

            // keep best
            if (bestFit < fitness) {
                bestFit = fitness;
                bestSeq = currentSeq;
            }

            // log
            std::cout << "\tANT " << k << "\n\t\tFit: " << fitness << "\n\t\tSeq: " << currentSeq << "\n";

        }
        std::cout << std::endl << "ITERATION " << iterations << " HAS ENDED!\n\n";
    }
    std::cerr << "BEST: " << bestFit << "\t" << bestSeq << "\n";
    std::cerr << "SIZE: " << PHEROMONES.size() << "\n";
    return bestSeq;
}