#ifndef GENERATOR_HH
#define GENERATOR_HH

template <typename It, typename Fcn>
class generator_ {
  const It _begin, _end;
  mutable Fcn f;

  class iterator {
    friend class generator_;
    Fcn& f;
    It it;

    iterator(Fcn& f, It it): f(f), it(it) { }

    template <typename T, typename = void>
    struct impl_ {
      using diff_type = typename std::iterator_traits<T>::difference_type;
      static decltype(auto) get(const iterator& it) { return it.f(*it.it); }
      static diff_type distance(const It& a, const It& b) {
        using std::distance;
        return distance(a,b);
      }
      static T next(const T& it, diff_type n = 1) {
        using std::next;
        return next(it,n);
      }
    };
    template <typename T>
    struct impl_<T, std::enable_if_t<std::is_integral<T>::value>> {
      using diff_type = T;
      static decltype(auto) get(const iterator& it) { return it.f(it.it); }
      static diff_type distance(T a, T b) { return b - a; }
      static T next(T it, diff_type n = 1) { return it + n; }
    };
    using impl = impl_<It>;

  public:
    using difference_type = typename impl::diff_type;

    decltype(auto) operator*() const { return impl::get(*this); }

    bool operator==(const iterator& o) const noexcept { return it == o.it; }
    bool operator!=(const iterator& o) const noexcept { return it != o.it; }
    iterator& operator++() { ++it; return *this; }

    friend difference_type distance(const iterator& a, const iterator& b) {
      return impl::distance(a.it,b.it);
    }
    friend iterator next(const iterator& it, difference_type n = 1) {
      return { it.f, impl::next(it.it,n) };
    }
  };
  friend class iterator;

public:
  template <typename B, typename E, typename F>
  generator_(B&& b, E&& e, F&& f)
  : _begin(std::forward<B>(b)),
    _end(std::forward<E>(e)),
    f(std::forward<F>(f)) { }

  iterator begin() const noexcept { return { f, _begin }; }
  iterator end  () const noexcept { return { f, _end   }; }

  friend iterator begin(const generator_& g) noexcept { return g.begin(); }
  friend iterator end  (const generator_& g) noexcept { return g.end(); }
};

template <typename B, typename E, typename F>
generator_<B,F> generator(B&& b, E&& e, F&& f) {
  return { std::forward<B>(b), std::forward<E>(e), std::forward<F>(f) };
}

#endif
