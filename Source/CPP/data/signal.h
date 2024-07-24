#ifndef SIGNAL_H
#define SIGNAL_H

enum class Sigs {
  Stop,
};

template <Sigs Sig>
struct Signal {};

#endif
