#ifndef LOG_H
#define LOG_H
#include <iostream>

// log levels (comment to deactivate)
#define LOG
#define LOG_INFO
#define LOG_WARN
#define LOG_ERR
#define LOG_TODO

// info groups active
#define INFO_GRPS 0

// colors and metafunctions
#if defined(LOG)
#define BRED "\033[1;31m"
#define BYEL "\033[1;33m"
#define BGRN "\033[1;32m"
#define BBLU "\033[1;34m"
#define CRESET "\033[0m"

template <size_t ...Ids>
struct LOG_GRP_LST {};

template <size_t Id, typename Lst>
struct LOG_INFO_ID_ACTIVE;

template <size_t Id> struct LOG_INFO_ID_ACTIVE<Id, LOG_GRP_LST<>> {
  static constexpr bool value = false;
};

template <size_t Id, size_t Head, size_t... Ids>
struct LOG_INFO_ID_ACTIVE<Id, LOG_GRP_LST<Head, Ids...>> {
  static constexpr bool value =
      Id == Head || LOG_INFO_ID_ACTIVE<Id, LOG_GRP_LST<Ids...>>::value;
};
#endif

// INFO is displayed all the time when LOG_INFO is defined
// INFO_GRP takes the id of a group, only the groups that are in INFO_GRPS are
// displayed
#if defined(LOG) && defined(LOG_INFO)
#define INFO_GRP(msg, id)                                                      \
  {                                                                            \
    if constexpr (LOG_INFO_ID_ACTIVE<id, LOG_GRP_LST<INFO_GRPS>>::value) {     \
      std::cout << BBLU "INFO: " CRESET msg << std::endl;                      \
    }                                                                          \
  }
#define INFO(msg) INFO_GRP(msg, 0)
#else
#define INFO(msg)
#define INFO_GRP(msg, id)
#endif

// WARN
#if defined(LOG) && defined(LOG_WARN)
#define WARN(msg) std::cout << BYEL "WARN: " CRESET msg << std::endl;
#else
#define WARN(msg)
#endif

// ERROR
#if defined(LOG) && defined(LOG_ERR)
#define ERROR(msg) std::cout << BRED "ERROR: " CRESET msg << std::endl;
#else
#define ERROR(msg)
#endif

// TODO
#if defined(LOG) && defined(LOG_TODO)
#define TODO(msg) std::cout << BGRN "TODO: " CRESET msg << std::endl;
#else
#define TODO(msg)
#endif

#endif
