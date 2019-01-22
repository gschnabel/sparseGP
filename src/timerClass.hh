
class Timer {
public:
  Timer() : sys_time("Sys.time") { Reset(); }
  void Start()  { start_t = getFractionalSeconds() ; }
  void step(std::string msg) {
    end_t = getFractionalSeconds();
    elapsed = end_t - start_t;              // Calculate elapsed time in seconds
    if (!silent)
      Rcpp::Rcout << msg << " " << elapsed << std::endl;
  }
  void Silent() {
    silent = true;
  }
  void Stop() {
    end_t = getFractionalSeconds();
    elapsed = end_t - start_t;              // Calculate elapsed time in seconds
    cumul += elapsed;
  }
  void Reset() { end_t = start_t = elapsed = cumul = 0.0; silent=false; }
  double ElapsedTime() { return elapsed; }
  double CumulativeTime() { return cumul; }
private:
  Rcpp::Function sys_time ;
  bool silent;
  double start_t, end_t, elapsed, cumul;
  double getFractionalSeconds(void) {
    return Rcpp::as<double>( sys_time() ) ;
  }
};
