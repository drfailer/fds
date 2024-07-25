#ifndef CONSTANTS_H
#define CONSTANTS_H
#include <cstddef>

constexpr size_t LABEL_LENGTH = 60;
constexpr size_t MAX_MATERIALS = 20;
constexpr size_t MAX_DIM = 3;
constexpr size_t CC_MAXCROSS_EDGE = 10;
constexpr size_t LOW_IND = 1;
constexpr size_t HIGH_IND = 2;
constexpr size_t MESH_STRING_LENGTH = 2;

// note should be moved somewhere else
constexpr int dim(int a, int b) { return b - a + 1; }

#endif
