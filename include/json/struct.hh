#ifndef NLOHMANN_JSON_STRUCT_SERIALIZER_HH
#define NLOHMANN_JSON_STRUCT_SERIALIZER_HH

#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/seq/for_each.hpp>

#define JSON_SERIALIZER_GET(r, data, val) \
  try { \
    j.at( BOOST_PP_STRINGIZE(val) ).get_to(x.val); \
  } catch(const json::out_of_range&) { }

#define JSON_SERIALIZER_SET(r, data, val) \
  j[ BOOST_PP_STRINGIZE(val) ] = x.val;

#define JSON_SERIALIZER(NAME, VALUES) \
namespace nlohmann { \
template <> struct adl_serializer<NAME> { \
  static void from_json(const json& j, NAME& x) { \
    BOOST_PP_SEQ_FOR_EACH(JSON_SERIALIZER_GET, , VALUES) \
  } \
  static void to_json(json& j, const NAME& x) { \
    BOOST_PP_SEQ_FOR_EACH(JSON_SERIALIZER_SET, , VALUES) \
  } \
}; \
}

#endif
