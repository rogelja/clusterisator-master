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


double CalculSeuil(double ,RegisteredInstance );


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
	//usage();
	//liste des instances a tester
	std::list<AvailableInstances> list_instance;

	//on travaille sur l'instance segmentation
	list_instance.push_back(AvailableInstances::segmentation);




	//fichiers relevant les données de sortie et les statistiques.
	std::ofstream file("output.csv");
	std::ofstream stats("Stats.csv");



	// les informations qui nous intéressent
	file << "id;nom;n;k;Run;RunMax;Seed;LevelMax;";
	file << "ite0;score0;CPU0;";
	file << "start;";
	file << "ite;score;CPU;";
	//file << "Rapport Scores (%);Rapport iterations total(%); rapport iterations moyen(%);Rapport temps(%)";
	file << std::endl;



	//fichier de sortie pour les erreurs
	std::ofstream debug("debug.log");

	//Nb de Run lancer pour les stats
	size_t const nbLancer ( 10);


	// définit le nombre de classe de départ, celle de fin et le pas d'incrément	
	size_t const kmin (3 );
	size_t const kstep( 1 );
	size_t const kmax ( 10);



	//On lance les runs sur toutes les instances de la liste.
	for (auto const & id : list_instance) {
		RegisteredInstance instance(id);
		instance.out();

		std::map<size_t, double> sumsIte_moyenne[kmax/kstep];
		std::map<size_t, size_t> nbIte_moyenne[kmax/kstep];
		std::map<size_t, double> score_moyenne[kmax/kstep];
		std::map<size_t, double> sumsTime_moyenne[kmax/kstep];
		std::map<size_t, double> sumsIte_variance[kmax/kstep];
		std::map<size_t, size_t> nbIte_variance[kmax/kstep];
		std::map<size_t, double> score_variance[kmax/kstep];
		std::map<size_t, double> sumsTime_variance[kmax/kstep];

		std::map<size_t, MultiLevelAlgoStats> allStats[kmax/kstep][nbLancer];
		std::map<size_t, double> sumsIte[kmax/kstep][nbLancer];
		std::map<size_t, size_t> nbIte[kmax/kstep][nbLancer];
		std::map<size_t, double> score[kmax/kstep][nbLancer];
		std::map<size_t, double> sumsTime[kmax/kstep][nbLancer];
		size_t NBLevel[kmax/kstep];
		size_t seed[nbLancer];


		for( size_t Lancer=0; Lancer < nbLancer ;Lancer++){



			for (size_t k(kmin); k <= kmax; k+= kstep) {	




				seed[Lancer]=(k-kmin)/kstep+1;

				Number::SetSeed(seed[Lancer]);

				std::cout << " ------------------------"<<k<<std::endl;
				MultiLevelAlgo algo(instance, k);
				//algo.setSeuil(seuil);				
				algo.setOut(debug);
				algo.buildMultiLevelData_tab(20 * k);
				NBLevel[k/kstep-1]=algo.nbLevels();


				//algo.buildMultiLevelData_tab(20);
				// on agrege 20k des noeuds par palier de 5% des noeuds totaux 
				//algo.buildMultiLevelData(20 * k, amax);
				Partition start(instance.nbObs(), k);
				//generation du point de depart
				algo.getStartPoint(start);		
				//on ne s'occupe aps encore du pas
				//size_t const stepMax((size_t) std::ceil(algo.nbLevels() * 0.10));
				for (size_t level(0); level <= algo.nbLevels(); ++level) {
					algo.setStartPoint(start);
					algo.setStartLevel(level);
					algo.launch();
					allStats[k/kstep-1][Lancer][level] = algo.stats();
				}


				//on calcule les différentes valeurs intéréssantes pour nos analyses
				for (auto const & stat : allStats[k/kstep-1][Lancer]) {
					size_t const level(stat.first);
					nbIte[k/kstep-1][Lancer][level] = allStats[k/kstep-1][Lancer][level].size();
					score[k/kstep-1][Lancer][level] = allStats[k/kstep-1][Lancer][level].begin()->second._cost;
					for (auto & stat : allStats[k/kstep-1][Lancer][level]) {
						sumsIte[k/kstep-1][Lancer][level]  += stat.second._ite;
						sumsTime[k/kstep-1][Lancer][level] += stat.second._time;
					}
				}
			}
		}


		//Cette partie calcule les moyennes et variances des données relevées pendant les tests et les
		//rédige dans les fichiers output.csv et stats.csv .

		for (size_t k(kmin); k <= kmax; k+= kstep) {
			size_t const levelMax(NBLevel[k/kstep-1]);
			for( size_t level=0; level <= levelMax ;level++){
				for( size_t Lancer=0; Lancer < nbLancer ;Lancer++){
					WriteCsv(file, id, 6);
					WriteCsv(file, instance.name);
					WriteCsv(file, instance.nbObs());
					WriteCsv(file, k);
					WriteCsv(file, Lancer);
					WriteCsv(file, nbLancer);
					WriteCsv(file, levelMax);
					WriteCsv(file, seed[Lancer]);

					WriteCsv(file, sumsIte[k/kstep-1][Lancer][0]);
					WriteCsv(file, score[k/kstep-1][Lancer][0], 15);
					WriteCsv(file, sumsTime[k/kstep-1][Lancer][0], 6);
					WriteCsv(file, level);
					WriteCsv(file, sumsIte[k/kstep-1][Lancer][level]);
					WriteCsv(file, score[k/kstep-1][Lancer][level], 15);
					WriteCsv(file, sumsTime[k/kstep-1][Lancer][level], 6);
					file << std::endl;


					sumsIte_moyenne[k/kstep-1][level] += sumsIte[k/kstep-1][Lancer][level]/((double)nbLancer);
					nbIte_moyenne[k/kstep-1][level] += nbIte[k/kstep-1][Lancer][level]/nbLancer;
					score_moyenne[k/kstep-1][level] += score[k/kstep-1][Lancer][level]/((double)nbLancer);
					sumsTime_moyenne[k/kstep-1][level] += sumsTime[k/kstep-1][Lancer][level]/((double)nbLancer);
				}
			}
		}
		file << std::endl;

		// calcul de la variance par la formule E[(X- E[X])²]


		for (size_t k(kmin); k <= kmax; k+= kstep) {
			size_t const levelMax(NBLevel[k/kstep-1]);
			for( size_t level=0; level <= levelMax ;level++){

				for( size_t Lancer=0; Lancer < nbLancer ;Lancer++){

					sumsIte_variance[k/kstep-1][level] += (sumsIte[k/kstep-1][Lancer][level]-sumsIte_moyenne[k/kstep-1][level])*
						(sumsIte[k/kstep-1][Lancer][level]-sumsIte_moyenne[k/kstep-1][level])/
						((double)nbLancer);
					nbIte_variance[k/kstep-1][level] += (nbIte[k/kstep-1][Lancer][level]-nbIte_moyenne[k/kstep-1][level])*
						(nbIte[k/kstep-1][Lancer][level]-nbIte_moyenne[k/kstep-1][level])/
						nbLancer;
					score_variance[k/kstep-1][level] += (score[k/kstep-1][Lancer][level]-score_moyenne[k/kstep-1][level])*
						(score[k/kstep-1][Lancer][level]-score_moyenne[k/kstep-1][level])/
						((double)nbLancer);
					sumsTime_variance[k/kstep-1][level] += (sumsTime[k/kstep-1][Lancer][level]-sumsTime_moyenne[k/kstep-1][level])*
						(sumsTime[k/kstep-1][Lancer][level]-sumsTime_moyenne[k/kstep-1][level])/
						((double)nbLancer);


				}
			}
		}

		stats << "id;nom;n;k;NbRun;nbLevel;";
		stats << "ite0_moyenne;ite0_variance;score0_moyenne;score0_variance;CPU0_moyenne;CPU0_variance;";
		stats << "start;";
		stats << "ite_moyenne;ite_variance;score_moyenne;score_variance;CPU_moyenne;CPU_variance;";
		stats << std::endl;

		for (size_t k(kmin); k <= kmax; k+= kstep) {
			size_t const levelMax(NBLevel[k/kstep-1]);
			for( size_t level=0; level <= levelMax ;level++){

				WriteCsv(stats, id, 6);
				WriteCsv(stats, instance.name);
				WriteCsv(stats, instance.nbObs());
				WriteCsv(stats, k);
				WriteCsv(stats, nbLancer);
				WriteCsv(stats, levelMax);

				WriteCsv(stats, sumsIte_moyenne[k/kstep-1][0]);
				WriteCsv(stats, sumsIte_variance[k/kstep-1][0]);
				WriteCsv(stats, score_moyenne[k/kstep-1][0], 15);
				WriteCsv(stats, score_variance[k/kstep-1][0], 15);
				WriteCsv(stats, sumsTime_moyenne[k/kstep-1][0], 6);
				WriteCsv(stats, sumsTime_variance[k/kstep-1][0], 6);
				WriteCsv(stats, level);
				WriteCsv(stats, sumsIte_moyenne[k/kstep-1][level]);
				WriteCsv(stats, sumsIte_variance[k/kstep-1][level]);
				WriteCsv(stats, score_moyenne[k/kstep-1][level], 15);
				WriteCsv(stats, score_variance[k/kstep-1][level], 15);
				WriteCsv(stats, sumsTime_moyenne[k/kstep-1][level], 6);
				WriteCsv(stats, sumsTime_variance[k/kstep-1][level], 6);

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
