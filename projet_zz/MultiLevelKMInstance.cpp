/*
* MultiLevelKMInstance.cpp
*
*  Created on: 15 dÃƒÂ©c. 2012
*      Author: manuel
*/

#include "../projet_zz/MultiLevelKMInstance.h"
#include "../src/HMeans.hpp"
#include "../src/KMAlgo.hpp"
#include "../src/Timer.hpp"
#include <iostream>
#include <fstream>

MultiLevelAlgo::MultiLevelAlgo(KMInstance const & instance, size_t k) :
	_instance(instance), _input(_instance, k), _startPoint(
	_instance.nbObs(), k) {
		// par defaut sur la sortie standard.
		setOut();
		_step = 1;
}

MultiLevelAlgo::~MultiLevelAlgo() {
	for (auto & ptr : _multiLevelConstraints)
		delete ptr;
}
void MultiLevelAlgo::buildInstance(size_t level, KMInstance & instance,
								   Aggregations & aggregations) {
									   _instance.mustLinks().clear();
									   _instance.cannotLinks().clear();
									   for (size_t i(0); i < level; ++i) {
										   for (auto const & ctr : *_multiLevelConstraints[i]) {
											   _instance.addMustLink(ctr.first, ctr.second);
										   }
									   }
									   _instance.buildMustLink(aggregations);
									   instance = KMInstance(_instance, aggregations);
}





void MultiLevelAlgo::buildMultiLevelData_tab(double nbNodes) {
#if 1
	// multimap pour avoir des doublons, on en a une par noeud qui contient le classement
	std::vector< std::multimap<Double,size_t > > TabDistance(_instance.nbObs());	
	Double const p(10);
	// On calcule toutes les distances entres les points 
	for (size_t i(0); i < _instance.nbObs(); ++i){
		for (size_t j(i+1); j < _instance.nbObs(); ++j){
			// calcul de la distance
			Double const d(_instance.distance(i,j));
			// insertion dans la map
			TabDistance[i].insert(std::make_pair(d, j));
			TabDistance[j].insert(std::make_pair(d, i));
			

			size_t n(0);
			std::multimap<Double,size_t >::iterator it(TabDistance[i].begin());
			// on ne garde que les 10 plus proches voisins 
			while( n<p && it != TabDistance[i].end()){
				++it;
				++n;
			}
			TabDistance[i].erase(it, TabDistance[i].end());
		}
	}




	KMPartition partition(_instance, _instance.nbObs());

	// on crée les singletons
	for (size_t i(0); i < _instance.nbObs(); ++i)
		partition.shift(i, i);
	
	
	//on choisit de construire au maximum 3 niveaux ou un nb de points restants supérieur supérieur a nbNodes
	
	while (partition.nbLabels() > nbNodes &&  nbLevels() < 2) {
		bool create(false);
		IndexedList points(_instance.nbObs(), true);
		while(! points.empty() && partition.nbLabels() > nbNodes &&  nbLevels() < 2 ){
			size_t const point(points.pop_random());
			for(auto const & voisin : TabDistance[point]){
				if(points.contains(voisin.second)){
					if(!create){
						create = true;
						_multiLevelConstraints.push_back(new KMConstraints(_input.nbObs()));
					}
					points.erase(voisin.second);					
					partition.fusion(partition.label(point), partition.label(voisin.second));
					_multiLevelConstraints.back()->newCtr(point, voisin.second);
					break;
				}
				points.erase(point);
			}
		}
	};
	std::cout << "nbNodes     " << partition.nbLabels() << std::endl;
	std::cout << "Built       " << nbLevels() << " aggregation levels"	<< std::endl;
#else
	_multiLevelConstraints.push_back(new KMConstraints(_input.nbObs()));
	_multiLevelConstraints.back()->newCtr(0,2);
	_multiLevelConstraints.back()->newCtr(1,2);
#endif
}








//
// _step: 
// _startLevel : niveau de dÃƒÂ©part pour le raffinement
// _startPoint : (attention doit ÃƒÂªtre compatible avec le niveau de plus agrÃƒÂ©gÃƒÂ©)
void MultiLevelAlgo::refine() {


	KMInstance instance;
	Aggregations aggregations;
	Timer timer;

	//Pour chaque level construit 

	for (int level(_startLevel); level >= 0; level -= _step) {
		
		//on reconstruit le graphe avec les points agrégés
		buildInstance(level, instance, aggregations);


		KMInput input(instance, _input.maxNbLabels());


		// On applique le dernier partitionnement trouvé
		for (size_t i(0); i < aggregations.v.size(); ++i) {
			input.shiftForced(i, _input.label(*aggregations.v[i].begin()));
		}


		//on initialise le compteur
		timer.restart();

		// On applique les K-means pour rafiner le partitionnement
		HMeans<true>()(input);

		// On sauvegarde les mesures dans un tableau recapitulatif
		_stats[level] = MultiLevelAlgoStat(level, input.ite(), timer.elapsed(),
			input.cost());

		//On replace les points dans la nouvelle classe qu'on lui a trouvé
		for (size_t i(0); i < aggregations.v.size(); ++i) {
			for (auto const & j : aggregations.v[i])
				_input.shiftForced(j, input.label(i));
		}

	}


}
// @brief lance l'algorithme multi-niveau Ã  partir d'une suite d'agrÃ©gation, d'une partition de dÃ©part et d'un niveau de dÃ©part et avec un pas donnÃ©
// _step: 
// _startLevel : niveau de dÃƒÂ©part pour le raffinement
// _startPoint : (attention doit ÃƒÂªtre compatible avec le niveau de plus agrÃƒÂ©gÃƒÂ©)

void MultiLevelAlgo::setOut(std::ostream & stream) {
	_out = &stream;
}
void MultiLevelAlgo::launch() {
	_stats.clear();
	
	std::cout << "Step  : " << _step << std::endl;
	std::cout << "Level : " << _startLevel << std::endl;
	
	
	// initialisation au point de depart
	for (size_t i(0); i < _input.nbObs(); ++i) {
		_input.shiftForced(i, _startPoint.label(i));
	}


	std::cout << "Starting point value : " << _input.computeCost() << std::endl;

	// lancement du run 
	refine();

}


size_t MultiLevelAlgo::nbLevels() const {
	return _multiLevelConstraints.size();
}



void MultiLevelAlgo::getStartPoint(Partition & point) {
	KMInstance instance;
	Aggregations aggregations;
	buildInstance(nbLevels(), instance, aggregations);
	KMInput input(instance, _input.maxNbLabels());
	input.random(0);
	std::cout << "Starting point value : " << input.computeCost() << std::endl;
	for (size_t i(0); i < _input.nbObs(); ++i) {
		point.shift(i, input.label(aggregations.newIds[i]));
	}

}

void MultiLevelAlgo::setStartLevel(size_t level) {
	_startLevel = level;
}

void MultiLevelAlgo::setStep(size_t step) {
	_step = step;
}

void MultiLevelAlgo::setStartPoint(Partition & point) {
	_startPoint = point;
}

std::ostream & MultiLevelAlgo::out() {
	return *_out;
}
MultiLevelAlgoStats const & MultiLevelAlgo::stats() const {
	return _stats;
}


