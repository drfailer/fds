#ifndef FDS_H
#define FDS_H

extern "C" {

void start(const char *);
void mainLoopBegin();
void runPredictor();
void runCorrector();
void dumpOutputFiles();
int stopMainLoop();
void finish();

}

#endif
