
// A Scalar that asserts for uninitialized access.
template<typename T>
class SafeScalar {
 public:
  SafeScalar() : initialized_(false) {}
  SafeScalar(const SafeScalar& other) {
    *this = other;
  }
  SafeScalar& operator=(const SafeScalar& other) {
    val_ = T(other);
    initialized_ = true;
    return *this;
  }
  
  SafeScalar(T val) : val_(val), initialized_(true) {}
  SafeScalar& operator=(T val) {
    val_ = val;
    initialized_ = true;
  }
  
  operator T() const {
    VERIFY(initialized_ && "Uninitialized access.");
    return val_;
  }
 
 private:
  T val_;
  bool initialized_;
};
