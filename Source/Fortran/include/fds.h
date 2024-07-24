#ifndef FDS_H
#define FDS_H

extern "C" {

void start(const char *);

void mainLoopBegin();

void runPredictor();
void processExternalVariables();
void predictorFinitDifference();

// change time step loop functions
void computeDensity();
void exchangeValues();
void computeVelocity();
void hvacSolver();
void computeWallBC();
void exchangeMeshDivergence();
void updateGlobalPressure();
void computeDivergencePhase2();
void solvePressure();
void predictVelocity();
int verifyInstability();
void updateTimeStepIndex();
int timeStepReduced();

void runCorrector();
void dumpOutputFiles();
int stopMainLoop();

void startPredictor();
void finishPredictor();

void finish();

}

#endif
