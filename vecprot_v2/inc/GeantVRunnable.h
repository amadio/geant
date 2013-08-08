#ifndef GEANT_VRUNNABLE
#define GEANT_VRUNNABLE

//______________________________________________________________________________
//   GeantVRunnable - A general interface to concurrent runnables. A runnable
//                    is an object embedding data which is processed by the Run()
//                    method. Running the object may produce a new runnable.
//______________________________________________________________________________

class GeantVRunnable {
public:
   GeantVRunnable();
   ~GeantVRunnable();
   virtual GeantVRunnable *Run() = 0;

   ClassDef(GeantVRunnable, 1)      // ABC for concurrent runnables
};
