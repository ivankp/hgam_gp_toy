#ifndef NLOHMANN_JSON_BOOST_OPTIONAL_HH
#define NLOHMANN_JSON_BOOST_OPTIONAL_HH

namespace nlohmann {
template <typename T> struct adl_serializer<boost::optional<T>> {
  static void from_json(const json& j, boost::optional<T>& opt) {
    if (j.is_null()) opt = boost::none;
    else opt = j.get<T>();
  }
  static void to_json(json& j, const boost::optional<T>& opt) {
    if (opt == boost::none) j = nullptr;
    else j = *opt;
  }
};
}

#endif
