#include <map>
#include <string>
#include <time.h>
#include <sys/time.h>
#include <stdint.h>
//--------------------------------------------------------------------------------------------------
#ifndef _ZEITGEIST_HPP
#define _ZEITGEIST_HPP
//--------------------------------------------------------------------------------------------------
#ifdef ZEITGEIST
//#define ZeitGeist_define(NAME) static ZeitGeist& NAME ## _zg = ZeitGeist::zg(#NAME)
#ifdef ZEITGEIST_DEV
#define ZeitGeist_define(NAME) ZeitGeist& NAME ## _zg = ZeitGeist::zg(#NAME)
#define ZeitGeist_tick(NAME) NAME ## _zg.Tick(#NAME)
#define ZeitGeist_tock(NAME) NAME ## _zg.Tock(#NAME)
#else
#define ZeitGeist_define(NAME) static ZeitGeist& NAME ## _zg = ZeitGeist::zg(#NAME)
#define ZeitGeist_tick(NAME) NAME ## _zg.Tick()
#define ZeitGeist_tock(NAME) NAME ## _zg.Tock()
#endif
#define ZeitGeist_inc(NAME) NAME ## _zg.Inc()
#define ZeitGeist_add(NAME, I) NAME ## _zg.Inc(I)
#define ZeitGeist_name(NAME) NAME ## _zg
#else
#define ZeitGeist_define(NAME)
#define ZeitGeist_tick(NAME)
#define ZeitGeist_tock(NAME)
#define ZeitGeist_inc(NAME)
#define ZeitGeist_add(NAME, I)
#define ZeitGeist_reset()
#define ZeitGeist_name(NAME) NAME ## _zg
#endif
//------------------------------------------------------------------------------
class ZeitGeist {
private:
#ifdef ZEITGEIST_DEV
  static std::string prefix;
#endif
  int64_t tick;
  bool ticktock;
  int64_t cnt;
  int64_t tot;
  int64_t min;
  int64_t max;
  timeval tmp_tv;
  inline int64_t usec () {
    return tmp_tv.tv_sec * 1000000 + tmp_tv.tv_usec;
  }
public:
  inline static ZeitGeist& zg (std::string str) {
#ifdef ZEITGEIST_DEV
    return ZeitGeist::zgMap()[prefix+str];
#else
    return ZeitGeist::zgMap()[str];
#endif
  }
  static std::map<std::string, ZeitGeist>& zgMap () {
    static std::map<std::string, ZeitGeist> _zg;
    return _zg;
  }

  ZeitGeist () {
    ticktock = false;
    cnt = 0;
    min = 0;
    max = 0;
    tot = 0;
  }
#ifdef ZEITGEIST_DEV
  inline void Tick (std::string name) {
    prefix += name + ":";
#else
  inline void Tick () {
#endif
    gettimeofday(&tmp_tv, NULL);
    tick = usec();
    ticktock = true;
  }
#ifdef ZEITGEIST_DEV
  inline void Tock(std::string name) {
#else
  inline void Tock() {
#endif
    gettimeofday(&tmp_tv, NULL);
    ++cnt;
    tot += usec() - tick;
#ifdef ZEITGEIST_DEV
    prefix.erase(prefix.size()-name.size()-1);
#endif
  }
  inline void Inc(int inc=1) {
    cnt+=inc;
  }
  int64_t Tock_stat () {
    gettimeofday(&tmp_tv, NULL);
    if (!ticktock) return 0;
    int64_t diff = usec() - tick;
    ++cnt;
    if (cnt == 1) {
      min = max = diff;
    } else {
      if (diff < min) min = diff;
      if (diff > max) max = diff;
    }
    tot += diff;
    ticktock = false;
    return diff;
  }
  void Reset () {
    ticktock = false;
    cnt = 0;
    min = 0;
    max = 0;
    tot = 0;
  }
  const int64_t& Min () const { return min; }
  const int64_t& Max () const { return max; }
  const int64_t& Total () const { return tot; }
  double Total_s () const { return double(tot)*0.000001; }
  const int64_t& Count () const { return cnt; }
  double Avg () const {
    if (cnt == 0) return 0;
    return double(tot)/double(cnt);
  }
};
//--------------------------------------------------------------------------------------------------
#endif // _ZEITGEIST_HPP
