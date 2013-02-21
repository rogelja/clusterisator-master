/*
* main.cpp
*
*  Created on: 15 déc. 2012
*      Author: manuel
*/

#include "src/common.h"
#include "src/Partition.hpp"

#include "src/KMAlgo.hpp"
#include "src/KMInstance.hpp"
#include "src/KMConstraints.hpp"
#include "src/KMPartition.hpp"

#include "src/Number.hpp"
#include "src/RegisteredInstance.hpp"

#include "projet_zz/MultiLevelKMInstance.h"
#include <random>
#include <iostream>
#include <fstream>

void usage() {
	std::cout << "Available instances : \n";
	for (size_t i(0); i < AvailableInstances::SIZE; ++i) {
		AvailableInstances id(static_cast<AvailableInstances>(i));
		RegisteredInstance instance(id);
		std::cout << std::setw(3) << i;
		std::cout << " : ";
		std::cout << std::setw(30) << std::left << instance.name << std::right;

		std::cout << "\n";
	}
	std::cout << "<exe> <id of selected instance> <number of labels>\n";
	std::cout << "The program launch the multi level algorithm\n";

}
template<class T> void WriteCsv(std::ostream & stream, T const & t) {
	stream << t << ";";
}
template<class T> void WriteCsv(std::ostream & stream, T const & t, size_t n) {
	stream << std::setprecision(n) << t << ";";
}

int main(int argc, char ** argv) {

	//liste des instances a tester
	std::list<size_t> list_instance;
	//	list_instance.push_back(0);
	//list_instance.push_back(1);
	//list_instance.push_back(4);
	//list_instance.push_back(7);
	list_instance.push_back(8);
	//list_instance.push_back(9);
	//list_instance.push_back(12);
	//list_instance.push_back(15);
	//list_instance.push_back(16);
	//list_instance.push_back(19);
	//list_instance.push_back(22);
	//list_instance.push_back(25);
	//list_instance.push_back(28);
	//list_instance.push_back(29);

	// préparation de l'output qu iszera effacer à chaque run !!

	std::ofstream file("output.csv");
	std::ofstream stats("Stats.csv");

	file << "id;nom;n;k;amax;Run;nbLevel;graine;";
	file << "ite0;score0;CPU0;";
	file << "start;";
	file << "ite;score;CPU;";
	//file << "Rapport Scores (%);Rapport iterations total(%); rapport iterations moyen(%);Rapport temps(%)";
	file << std::endl;
	// ce programme lance les tests sur toutes les instances.
	std::ofstream debug("debug.log");

	size_t const nbLancer = 50;
	size_t const kmax = 100;
	size_t const pasLevel =10 ;




	for (auto const & i : list_instance) {
		AvailableInstances id(static_cast<AvailableInstances>(i));
		RegisteredInstance instance(id);
		instance.out();
		// on va tester notre algo pour un nombre de classe compris [2 , 15% Nbpoint] en incrementant de 1


		std::map<size_t, double> sumsIte_moyenne[kmax/pasLevel];
		std::map<size_t, size_t> nbIte_moyenne[kmax/pasLevel];
		std::map<size_t, double> score_moyenne[kmax/pasLevel];
		std::map<size_t, double> sumsTime_moyenne[kmax/pasLevel];
		std::map<size_t, double> sumsIte_variance[kmax/pasLevel];
		std::map<size_t, size_t> nbIte_variance[kmax/pasLevel];
		std::map<size_t, double> score_variance[kmax/pasLevel];
		std::map<size_t, double> sumsTime_variance[kmax/pasLevel];



		std::map<size_t, MultiLevelAlgoStats> allStats[kmax/pasLevel][nbLancer];
		std::map<size_t, double> sumsIte[kmax/pasLevel][nbLancer];
		std::map<size_t, size_t> nbIte[kmax/pasLevel][nbLancer];
		std::map<size_t, double> score[kmax/pasLevel][nbLancer];
		std::map<size_t, double> sumsTime[kmax/pasLevel][nbLancer];
		size_t NBLevel[kmax/pasLevel];
		size_t Graine[nbLancer];

		//size_t const kmax((size_t) std::ceil(instance.nbObs() * 0.15));


		size_t const amax((size_t) std::ceil(instance.nbObs() * 0.10));
		std::cout << "amax = " << amax << std::endl;
		// mais on arrete si notre algo n'est plus en multi-level.

		for( size_t Lancer=0; Lancer < nbLancer ;Lancer++){

			Graine[Lancer]=rand()+1;
			Number::SetSeed(Graine[Lancer]);
			for (size_t k(10); k <= kmax; k+= pasLevel) {
				MultiLevelAlgo algo(instance, k);
				algo.setOut(debug);
				// on agrege 20k des noeuds par palier de 5% des noeuds totaux 
				algo.buildMultiLevelData(20 * k, amax);
				Partition start(instance.nbObs(), k);
				//generation du point de depart
				algo.getStartPoint(start);


				NBLevel[k/pasLevel-1] = algo.nbLevels();




				//on ne s'occupe aps encore du pas
				//size_t const stepMax((size_t) std::ceil(algo.nbLevels() * 0.10));


				size_t stepMax=1;

				for (size_t step(1) ; step < stepMax+1 ; ++step)
				{
					algo.setStep(step);

					for (size_t level(0); level <= algo.nbLevels(); ++level) {
						algo.setStartPoint(start);
						algo.setStartLevel(level);
						algo.launch();
						allStats[k/pasLevel-1][Lancer][level] = algo.stats();
					}

					for (auto const & stat : allStats[k/pasLevel-1][Lancer]) {
						size_t const level(stat.first);
						nbIte[k/pasLevel-1][Lancer][level] = allStats[k/pasLevel-1][Lancer][level].size();
						score[k/pasLevel-1][Lancer][level] = allStats[k/pasLevel-1][Lancer][level].rbegin()->second._cost;
						for (auto & stat : allStats[k/pasLevel-1][Lancer][level]) {
							sumsIte[k/pasLevel-1][Lancer][level] += stat.second._ite;
							sumsTime[k/pasLevel-1][Lancer][level] += stat.second._time;




						}
					}

				}
			}
		}

		for (size_t k(10); k <= kmax; k+= pasLevel) {
			size_t const levelMax(NBLevel[k/pasLevel-1]);
			for( size_t level=0; level < levelMax ;level++){

				for( size_t Lancer=0; Lancer < nbLancer ;Lancer++){


					//double RapCost	( (scoreTot[level]-sumsCost[algo.nbLevels()])/sumsCost[algo.nbLevels()]);
					//double RapIte	( (sumsIteTot[level]-sumsIteTot[algo.nbLevels()])/sumsIteTot[algo.nbLevels()]);
					//double RapTime	( (sumsTimeTot[level]-sumsTimeTot[algo.nbLevels()])/sumsTimeTot[algo.nbLevels()]);
					//double RapIteMoy( (sumsIte[level]-sumsIte[algo.nbLevels()])/sumsIte[algo.nbLevels()]);


					WriteCsv(file, i, 6);
					WriteCsv(file, instance.name);
					WriteCsv(file, instance.nbObs());
					WriteCsv(file, k);
					WriteCsv(file, amax);
					WriteCsv(file, Lancer);
					WriteCsv(file, levelMax);
					WriteCsv(file, Graine[Lancer]);
					WriteCsv(file, sumsIte[k/pasLevel-1][Lancer][0]);
					WriteCsv(file, score[k/pasLevel-1][Lancer][0], 15);
					WriteCsv(file, sumsTime[k/pasLevel-1][Lancer][0], 6);
					WriteCsv(file, level);
					WriteCsv(file, sumsIte[k/pasLevel-1][Lancer][level]);
					WriteCsv(file, score[k/pasLevel-1][Lancer][level], 15);
					WriteCsv(file, sumsTime[k/pasLevel-1][Lancer][level], 6);
					file << std::endl;


					sumsIte_moyenne[k/pasLevel-1][level] += sumsIte[k/pasLevel-1][Lancer][level]/nbLancer;
					nbIte_moyenne[k/pasLevel-1][level] += nbIte[k/pasLevel-1][Lancer][level]/nbLancer;
					score_moyenne[k/pasLevel-1][level] += score[k/pasLevel-1][Lancer][level]/nbLancer;
					sumsTime_moyenne[k/pasLevel-1][level] += sumsTime[k/pasLevel-1][Lancer][level]/nbLancer;


					// calcul de la variance par la formule E[X²] - E[X]²

					sumsIte_variance[k/pasLevel-1][level] += (sumsIte[k/pasLevel-1][Lancer][level]*sumsIte[k/pasLevel-1][Lancer][level])/nbLancer-
						(sumsIte_moyenne[k/pasLevel-1][level]*sumsIte_moyenne[k/pasLevel-1][level]);
					nbIte_variance[k/pasLevel-1][level] += (nbIte[k/pasLevel-1][Lancer][level]*nbIte[k/pasLevel-1][Lancer][level])/nbLancer-
						(nbIte_moyenne[k/pasLevel-1][level]*nbIte_moyenne[k/pasLevel-1][level]);
					score_variance[k/pasLevel-1][level] += (score[k/pasLevel-1][Lancer][level]*score[k/pasLevel-1][Lancer][level])/nbLancer-
						(score_moyenne[k/pasLevel-1][level]*score_moyenne[k/pasLevel-1][level]);
					sumsTime_variance[k/pasLevel-1][level] += (sumsTime[k/pasLevel-1][Lancer][level]*sumsTime[k/pasLevel-1][Lancer][level])/nbLancer-
						(sumsTime_moyenne[k/pasLevel-1][level]*sumsTime_moyenne[k/pasLevel-1][level]);




				}
			}
		}
		file << std::endl;






		stats << "id;nom;n;k;amax;NbRun;nbLevel;";
		stats << "ite0_moyenne;ite0_variance;score0_moyenne;score0_variance;CPU0_moyenne;CPU0_variance;";
		stats << "start;";
		stats << "ite_moyenne;ite_variance;score_moyenne;score_variance;CPU_moyenne;CPU_variance;";
		stats << std::endl;

		for (size_t k(10); k <= kmax; k+= pasLevel) {
			size_t const levelMax(NBLevel[k/pasLevel-1]);
			for( size_t level=0; level < levelMax ;level++){

				WriteCsv(stats, i, 6);
				WriteCsv(stats, instance.name);
				WriteCsv(stats, instance.nbObs());
				WriteCsv(stats, k);
				WriteCsv(stats, amax);
				WriteCsv(stats, nbLancer);
				WriteCsv(stats, levelMax);

				WriteCsv(stats, sumsIte_moyenne[k/pasLevel-1][0]);
				WriteCsv(stats, sumsIte_variance[k/pasLevel-1][0]);
				WriteCsv(stats, score_moyenne[k/pasLevel-1][0], 15);
				WriteCsv(stats, score_variance[k/pasLevel-1][0], 15);
				WriteCsv(stats, sumsTime_moyenne[k/pasLevel-1][0], 6);
				WriteCsv(stats, sumsTime_variance[k/pasLevel-1][0], 6);
				WriteCsv(stats, level);
				WriteCsv(stats, sumsIte_moyenne[k/pasLevel-1][level]);
				WriteCsv(stats, sumsIte_variance[k/pasLevel-1][level]);
				WriteCsv(stats, score_moyenne[k/pasLevel-1][level], 15);
				WriteCsv(stats, score_variance[k/pasLevel-1][level], 15);
				WriteCsv(stats, sumsTime_moyenne[k/pasLevel-1][level], 6);
				WriteCsv(stats, sumsTime_variance[k/pasLevel-1][level], 6);

				stats << std::endl;
			}
		}


	}


	file.close();
	stats.close();
	debug.close();
	//	std::cout << "Press ENTER to continue...";
	//	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	return 0;
}
