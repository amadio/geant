#ifndef GEANT_RUN_MANAGER_H
#define GEANT_RUN_MANAGER_H

class GeantRunManager
{
private:
  int fNpropagator = 0; /** Number of propagators */
  int fNthreads    = 0; /** Number of threads per propagator */

public:
	GeantRunManager() {}
	GeantRunManager(unsigned int npropagators, unsigned int nthreads);

  bool Initialize();
	void RunSimulation();

	~GeantRunManager();
};

#endif // GEANT_RUN_MANAGER_H
