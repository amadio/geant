#ifndef GEANT_RUN_MANAGER_H
#define GEANT_RUN_MANAGER_H

#include "GeantVApplication.h"
#include "GeantVApplication.h"

class GeantRunManager
{
private:

	//general config
	int fNpropagator=0;
	int fNthreads=0;
	int fNevents=0;
	int fNtotal=0;


	GeantVApplication *fApplication;
	GeantVApplication *fStdApplication;

	GeantRunManager();

public:
	/**
	 * @brief      Conctructor of the RunManager to manage many GeantPropagator 
	 *
	 * @param      config  An array with the config {nbPropagator, nbThreadForEachPropagator, nbTrackForEachThread,NumberOfBufferedTarck}
	 */
	GeantRunManager(unsigned int config[4]=NULL);
	GeantRunManager(unsigned int nbPropagator,unsigned int  nbThreadForEachPropagator,unsigned int  nbTrackForEachThread, unsigned int NumberOfBufferedTarck);

	void RunSimulation(PhysicsProcess *process,PrimaryGenerator *generator,GeantVApplication application, const char *geomfile = "geometry.root", bool graphics = false, bool single = false);

	~GeantRunManager();
};


#endif // GEANT_RUN_MANAGER_H
