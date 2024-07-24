#ifndef LOG_H
#define LOG_H
#include <iostream>
#include <type_traits>

// log levels (comment to deactivate)
#define LOG
#define LOG_INFO
#define LOG_WARN
#define LOG_ERR
#define LOG_TODO
#define LOG_DBG

// groups
#define MAIN_LOOP_GRP 1
#define TIME_STEP_LOOP_INFO_GRP 2
#define PREDICTOR_GRP 3
#define CORRECTOR_GRP 4

// info groups active
#define INFO_GRPS 0

// colors and metafunctions
#if defined(LOG)
#define BRED "\033[1;31m"
#define BYEL "\033[1;33m"
#define BGRN "\033[1;32m"
#define BBLU "\033[1;34m"
#define MAG "\033[0;35m"
#define CRESET "\033[0m"

namespace logh {

template <size_t... Ids> struct GroupList {};

template <size_t Id, typename Lst> struct isInfoIdActive;

template <size_t Id> struct isInfoIdActive<Id, GroupList<>> {
  static constexpr bool value = false;
};

template <size_t Id, size_t Head, size_t... Ids>
struct isInfoIdActive<Id, GroupList<Head, Ids...>> {
  static constexpr bool value =
      Id == Head || isInfoIdActive<Id, GroupList<Ids...>>::value;
};

} // namespace logh
#endif

// INFO is displayed all the time when LOG_INFO is defined
// INFO_GRP takes the id of a group, only the groups that are in INFO_GRPS are
// displayed
#if defined(LOG) && defined(LOG_INFO)
#define INFO_GRP(msg, id)                                                      \
  if constexpr (logh::isInfoIdActive<id, logh::GroupList<INFO_GRPS>>::value) { \
    std::cout << BBLU "INFO: " CRESET << msg << std::endl;                     \
  }
#define INFO(msg) INFO_GRP(msg, 0)
#else
#define INFO(msg)
#define INFO_GRP(msg, id)
#endif

// WARN
#if defined(LOG) && defined(LOG_WARN)
#define WARN(msg) std::cout << BYEL "WARN: " CRESET << msg << std::endl;
#else
#define WARN(msg)
#endif

// ERROR
#if defined(LOG) && defined(LOG_ERR)
#define ERROR(msg) std::cout << BRED "ERROR: " CRESET << msg << std::endl;
#else
#define ERROR(msg)
#endif

// TODO
#if defined(LOG) && defined(LOG_TODO)
#define TODO(msg) std::cout << BGRN "TODO: " CRESET << msg << std::endl;
#else
#define TODO(msg)
#endif

// DBG
#if defined(LOG) && defined(LOG_DBG)
#define DBG(var)                                                               \
  if constexpr (!std::is_same_v<const char *, decltype(var)>) {                \
    std::cout << MAG "DBG: " CRESET #var " = " << var << std::endl;            \
  } else { /* if not variable */                                               \
    std::cout << MAG "DBG: " CRESET << var << std::endl;                       \
  }
#else
#define DBG(var)
#endif

#endif
