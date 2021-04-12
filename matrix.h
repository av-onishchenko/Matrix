#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

int gcdex(int a, int b, int &x, int &y) {
  if (a == 0) {
    x = 0;
    y = 1;
    return b;
  }
  int x_1 = 0;
  int y_1 = 0;
  int d = gcdex(b % a, a, x_1, y_1);
  x = y_1 - (b / a) * x_1;
  y = x_1;
  return d;
}

template<int N>
class Finite {
public:
    Finite() {
      value_ = 0;
    }

    Finite(int x) {
      value_ = x % N;
      if (value_ < 0) {
        value_ += N;
      }
    }

    int GetValue() const {
      return value_;
    }

    unsigned GetBack() const {
      int x_1 = 0;
      int y_1 = 0;
      gcdex(value_, N, x_1, y_1);
      x_1 %= N;
      if (x_1 < 0) {
        x_1 += N;
      }
      return x_1;
    }

    Finite<N> &operator-() {
      value_ = N - value_;
      return *this;
    }

    Finite<N> &operator+=(const Finite<N> &x) {
      value_ += x.GetValue();
      normalaze();
      return *this;
    }

    Finite<N> &operator-=(const Finite<N> &x) {
      value_ -= x.GetValue();
      normalaze();
      return *this;
    }

    Finite<N> &operator*=(const Finite<N> &x) {
      value_ *= x.GetValue();
      normalaze();
      return *this;
    }

    Finite<N> &operator/=(const Finite<N> &x) {
      value_ *= x.GetBack();
      normalaze();
      return *this;
    }

    Finite<N> &operator--() {
      value_ -= 1;
      normalaze();
      return *this;
    }

    Finite<N> &operator++() {
      value_ += 1;
      normalaze();
      return *this;
    }

private:
    long long value_;

    void normalaze() {
      value_ = value_ % N;
      if (value_ < 0) {
        value_ += N;
      }
    }
};

template<int N>
Finite<N> operator+(const Finite<N> &a, const Finite<N> &b) {
  Finite<N> c(a.GetValue());
  c += b;
  return c;
}

template<int N>
Finite<N> operator-(const Finite<N> &a, const Finite<N> &b) {
  Finite<N> c(a.GetValue());
  c -= b;
  return c;
}

template<int N>
Finite<N> operator*(const Finite<N> &a, const Finite<N> &b) {
  Finite<N> c(a.GetValue());
  c *= b;
  return c;
}

template<int N>
Finite<N> operator/(const Finite<N> &a, const Finite<N> &b) {
  Finite<N> c(a.GetValue());
  c /= b;
  return c;
}

template<int N>
bool operator==(const Finite<N> &a, const Finite<N> &b) {
  return (a.GetValue() == b.GetValue());
}

template<int N>
bool operator!=(const Finite<N> &a, const Finite<N> &b) {
  return !(a == b);
}

template<int N>
bool operator<=(const Finite<N> &a, const Finite<N> &b) {
  return (a.GetValue() <= b.GetValue());
}

template<int N>
bool operator>=(const Finite<N> &a, const Finite<N> &b) {
  return (a.GetValue() >= b.GetValue());
}


const long long base = 1000000000;


class BigInteger {
public:
    BigInteger() : nums_({0}), is_negative_(false) {}

    BigInteger(int x) {
      if (x < 0) {
        is_negative_ = true;
        x = -x;
      } else {
        is_negative_ = false;
      }
      while (x > 0) {
        nums_.push_back(x % base);
        x /= base;
      }
      if (nums_.size() == 0) {
        nums_.push_back(0);
      }
      Balance(*this);
    }


    void Balance(BigInteger &x) {
      while (x.nums_.size() > 1 && x.nums_.back() == 0) {
        x.nums_.pop_back();
      }
      if (x.nums_.size() == 1 && x.nums_.back() == 0) {
        x.is_negative_ = false;
      }
      if (x.nums_.size() == 0) {
        nums_.push_back(0);
      }
    }

    BigInteger(std::string str) {
      if (str[0] == '-') {
        str = str.substr(1);
        is_negative_ = true;
      } else {
        is_negative_ = false;
      }
      for (long long i = str.length(); i > 0; i -= 9) {
        if (i < 9) {
          nums_.push_back(atoi(str.substr(0, i).c_str()));
        } else {
          nums_.push_back(atoi(str.substr(i - 9, 9).c_str()));
        }
      }
      Balance(*this);
    }

    BigInteger &operator=(const BigInteger &a) {
      nums_ = a.nums_;
      is_negative_ = a.is_negative_;
      return *this;
    }

    BigInteger operator-() const {
      BigInteger a = *this;
      a.is_negative_ = !a.is_negative_;
      if (a.nums_.size() == 1 && a.nums_.back() == 0) {
        a.is_negative_ = false;
      }
      return a;
    }

    bool smaller(const BigInteger &x) const {
      if (is_negative_ != x.is_negative_) {
        return is_negative_;
      }
      if (nums_.size() != x.nums_.size()) {
        if (is_negative_) {
          return !(nums_.size() < x.nums_.size());
        } else {
          return nums_.size() < x.nums_.size();
        }
      }
      for (int i = nums_.size() - 1; i > -1; --i) {
        if (nums_[i] != x.nums_[i]) {
          if (is_negative_) {
            return !(nums_[i] < x.nums_[i]);
          } else {
            return nums_[i] < x.nums_[i];
          }
        }
      }
      return false;
    }

    explicit operator bool() const;

    BigInteger &operator+=(const BigInteger &x) {
      BigInteger a;
      if (is_negative_ && x.is_negative_) {
        a = -*this;
        a += -x;
        *this = -a;
        Balance(*this);
        return *this;
      }
      if (is_negative_) {
        a = x;
        a -= -*this;
        Balance(a);
        return *this = a;
      }
      if (x.is_negative_) {
        *this -= -x;
        Balance(*this);
        return *this;
      }
      int d = 0;
      for (int i = 0; i < int(std::max(nums_.size(), x.nums_.size())) || d > 0; ++i) {
        if (int(nums_.size()) == i) {
          nums_.push_back(0);
        }
        nums_[i] += d;
        if (int(x.nums_.size()) > i) {
          nums_[i] += x.nums_[i];
        }
        d = nums_[i] / base;
        nums_[i] %= base;
      }
      Balance(*this);
      return *this;
    }

    BigInteger &operator-=(const BigInteger &x) {
      BigInteger a;
      if (is_negative_ && x.is_negative_) {
        a = -*this;
        a -= -x;
        a = -a;
        *this = a;
        Balance(*this);
        return *this;
      }
      if (is_negative_) {
        a = -*this;
        a += x;
        *this = -a;
        Balance(*this);
        return *this;
      }
      if (x.is_negative_) {
        *this += -x;
        Balance(*this);
        return *this;
      }
      if (smaller(x)) {
        a = x;
        a -= *this;
        *this = -a;
        Balance(*this);
        return *this;
      }
      for (int i = 0; i < int(nums_.size()); ++i) {
        if (int(x.nums_.size()) > i) {
          nums_[i] -= x.nums_[i];
        }
        if (nums_[i] < 0) {
          nums_[i] += base;
          --nums_[i + 1];
        }
      }
      Balance(*this);
      return *this;
    }

    BigInteger &operator++() {
      *this += 1;
      return *this;
    }

    BigInteger &operator--() {
      *this -= 1;
      return *this;
    }

    BigInteger operator++(int) {
      BigInteger a = *this;
      ++*this;
      return a;
    }

    BigInteger operator--(int) {
      BigInteger a = *this;
      --*this;
      return a;
    }

    BigInteger &operator*=(const BigInteger &x) {
      BigInteger a;
      if (is_negative_ != x.is_negative_) {
        a.is_negative_ = true;
      }
      std::vector<long long> result(nums_.size() + x.nums_.size());
      for (int i = 0; i < int(nums_.size()); ++i) {
        int d = 0;
        for (int j = 0; j < int(x.nums_.size()) || d > 0; ++j) {
          result[i + j] += d;
          if (int(x.nums_.size()) > j) {
            result[i + j] += nums_[i] * x.nums_[j];
          }
          d = result[i + j] / base;
          result[i + j] %= base;
        }
      }
      a.nums_ = result;
      Balance(a);
      return *this = a;
    }

    void reverse(BigInteger &x) {
      for (int j = 0; j < int(x.nums_.size()) / 2; ++j) {
        std::swap(x.nums_[j], x.nums_[int(x.nums_.size()) - j - 1]);
      }
    }


    BigInteger &operator/=(const BigInteger &x) {
      BigInteger a;
      BigInteger b;
      BigInteger c = 0;
      if (c.nums_.size() != 0) {
        c.nums_.clear();
      }
      BigInteger result = 0;
      if (result.nums_.size() != 0) {
        result.nums_.clear();
      }
      if (is_negative_) {
        a = -*this;
      } else {
        a = *this;
      }
      if (x.is_negative_) {
        b = -x;
      } else {
        b = x;
      }
      if (is_negative_ != x.is_negative_) {
        result.is_negative_ = true;
      } else {
        result.is_negative_ = false;
      }
      for (int i = int(a.nums_.size()) - 1; i > -1; --i) {
        reverse(c);
        c.nums_.push_back(a.nums_[i]);
        reverse(c);
        Balance(c);
        BigInteger d;
        if (!(c.smaller(b))) {
          int left = 1;
          int right = base;
          while (right - left > 1) {
            int middle = (right + left) / 2;
            d = b;
            d *= middle;
            if (!(c.smaller(d))) {
              left = middle;
            } else {
              right = middle;
            }
          }
          result.nums_.push_back(left);
          d = b;
          d *= left;
          c -= d;
        } else {
          result.nums_.push_back(0);
        }
      }
      reverse(result);
      Balance(result);
      *this = result;
      return *this;
    }

    BigInteger &operator%=(const BigInteger &x) {
      BigInteger a = *this;
      a /= x;
      a *= x;
      *this -= a;
      Balance(*this);
      return *this;
    }

    std::string toString() const {
      std::string str;
      if (is_negative_) {
        str = "-";
      }
      if (nums_.size() != 0) {
        str += std::to_string(nums_[nums_.size() - 1]);
      }
      for (int i = nums_.size() - 2; i >= 0; --i) {
        int a = base / 10;
        int b = nums_[i];
        for (int j = 0; j < 9; ++j) {
          str += '0' + b / a;
          b %= a;
          a /= 10;
        }
      }
      return str;
    }

    double toDouble() const {
      double a = 0;
      double b = 1;
      for (size_t i = 0; i < size_t(nums_.size()); ++i) {
        a += nums_[i] * b;
        b *= base;
      }
      if (is_negative_) {
        return -a;
      }
      return a;
    }

    BigInteger &div_2() {
      int carry = 0;
      for (int i = int(nums_.size()) - 1; i > -1; --i) {
        long long cur = nums_[i] + carry * base;
        nums_[i] = cur / 2;
        carry = cur % 2;
      }
      Balance(*this);
      return *this;
    }

    int mod_2() const {
      return nums_[0] % 2;
    }

    BigInteger abs() {
      if (!is_negative_) {
        return *this;
      }
      BigInteger ans = *this;
      ans.is_negative_ = false;
      return ans;
    }

private:
    std::vector<long long> nums_;
    bool is_negative_;

};

bool operator<(const BigInteger &a, const BigInteger &b) {
  return a.smaller(b);
}

bool operator>(const BigInteger &a, const BigInteger &b) {
  return b < a;
}

bool operator==(const BigInteger &a, const BigInteger &b) {
  return !(a < b || b < a);
}

bool operator!=(const BigInteger &a, const BigInteger &b) {
  return !(a == b);
}

bool operator<=(const BigInteger &a, const BigInteger &b) {
  return (a < b || a == b);
}

bool operator>=(const BigInteger &a, const BigInteger &b) {
  return (a > b || a == b);
}

BigInteger::operator bool() const {
  return *this != 0;
}

BigInteger operator+(const BigInteger &a, const BigInteger &b) {
  BigInteger c = a;
  c += b;
  return c;
}

BigInteger operator-(const BigInteger &a, const BigInteger &b) {
  BigInteger c = a;
  c -= b;
  return c;
}

BigInteger operator*(const BigInteger &a, const BigInteger &b) {
  BigInteger c = a;
  c *= b;
  return c;
}

BigInteger operator/(const BigInteger &a, const BigInteger &b) {
  BigInteger c = a;
  c /= b;
  return c;
}

BigInteger operator%(const BigInteger &a, const BigInteger &b) {
  BigInteger c = a;
  c %= b;
  return c;
}

BigInteger operator ""_bi(unsigned long long a) {
  return BigInteger(a);
}

std::istream &operator>>(std::istream &in, BigInteger &a) {
  std::string str;
  in >> str;
  a = BigInteger(str);
  return in;
}

std::ostream &operator<<(std::ostream &out, const BigInteger &a) {
  out << a.toString();
  return out;
}

BigInteger gcd(BigInteger a, BigInteger b) {
  a = a.abs();
  b = b.abs();
  BigInteger res = 1;
  while (a != 0 && b != 0) {
    if (a.mod_2() == 0 && b.mod_2() == 0) {
      a.div_2();
      b.div_2();
      res *= 2;
    }
    if (a.mod_2() == 1 && b.mod_2() == 0) {
      b.div_2();
    }
    if (a.mod_2() == 0 && b.mod_2() == 1) {
      a.div_2();
    }
    if (a.mod_2() == 1 && b.mod_2() == 1) {
      (a > b ? a = (a - b).div_2() : b = (b - a).div_2());
    }
  }
  return (a == 0 ? b * res : a * res);
}

class Rational {
public:

    Rational(int x) : numerator_(x), denominator_(1) {}

    Rational() : numerator_(0), denominator_(1) {}

    Rational(const BigInteger x, const BigInteger y = 1) : numerator_(x), denominator_(y) {}


    Rational &operator=(const Rational &x) {
      numerator_ = x.numerator_;
      denominator_ = x.denominator_;
      Balance();
      return *this;
    }

    Rational operator-() const {
      Rational x = *this;
      x.numerator_ = -x.numerator_;
      return x;
    }

    Rational &operator+=(const Rational &x) {
      numerator_ = numerator_ * x.denominator_ + x.numerator_ * denominator_;
      denominator_ = denominator_ * x.denominator_;
      Balance();
      return *this;
    }


    Rational &operator-=(const Rational &x) {
      *this += (-x);
      Balance();
      return *this;
    }

    Rational &operator*=(const Rational &x) {
      numerator_ *= x.numerator_;
      denominator_ *= x.denominator_;
      Balance();
      return *this;
    }

    Rational &operator/=(const Rational &x) {
      numerator_ *= x.denominator_;
      denominator_ *= x.numerator_;
      Balance();
      return *this;
    }

    std::string toString() const {
      std::string str;
      str = numerator_.toString();
      if (denominator_ == 1) {
        return str;
      }
      str += '/' + denominator_.toString();;
      return str;
    }

    explicit operator double() const {
      double a = numerator_.toDouble();
      double b = denominator_.toDouble();
      return a / b;
    }

    std::string asDecimal(size_t precision = 0) {
      BigInteger n = numerator_;
      for (size_t i = 0; i < precision + 1; ++i) {
        n *= 10;
      }
      n /= denominator_;
      bool neg = false;
      if (n < 0) {
        neg = true;
        n = -n;
      }
      if (n % 10 >= 5) {
        n += 10;
      }
      n /= 10;
      std::string str = n.toString();

      for (int i = 0; i < int(str.size()) / 2; i++) {
        std::swap(str[i], str[str.size() - i - 1]);
      }

      while (precision > str.size()) {
        str += '0';
      }

      for (int i = 0; i < int(str.size()) / 2; i++) {
        std::swap(str[i], str[str.size() - i - 1]);
      }

      if (precision != 0 && precision <= str.size())
        str.insert(str.begin() + str.size() - precision, '.');

      if (str[0] == '.') {
        str = '0' + str;
      }
      if (neg) {
        str = '-' + str;
      }
      return str;
    }

    bool smaller(const Rational &x) const {
      return numerator_ * x.denominator_ < x.numerator_ * denominator_;
    }

private:
    BigInteger numerator_;
    BigInteger denominator_;

    void Balance() {
      BigInteger g = gcd(numerator_, denominator_);
      numerator_ /= g;
      denominator_ /= g;
      if (denominator_ < 0) {
        denominator_ = -denominator_;
        numerator_ = -numerator_;
      }
    }
};

Rational operator+(const Rational &a, const Rational &b) {
  Rational c(a);
  c += b;
  return c;
}

Rational operator-(const Rational &a, const Rational &b) {
  Rational c(a);
  c -= b;
  return c;
}

Rational operator*(const Rational &a, const Rational &b) {
  Rational c(a);
  c *= b;
  return c;
}

Rational operator/(const Rational &a, const Rational &b) {
  Rational c(a);
  c /= b;
  return c;
}

bool operator<(const Rational &a, const Rational &b) {
  return a.smaller(b);
}

bool operator>(const Rational &a, const Rational &b) {
  return b < a;
}

bool operator==(const Rational &a, const Rational &b) {
  return !(a < b || b < a);
}

bool operator!=(const Rational &a, const Rational &b) {
  return !(a == b);
}

bool operator<=(const Rational &a, const Rational &b) {
  return (a < b || a == b);
}

bool operator>=(const Rational &a, const Rational &b) {
  return (a > b || a == b);
}

std::istream &operator>>(std::istream &in, Rational &a) {
  int x;
  in >> x;
  a = static_cast<Rational>(x);
  return in;
}

std::ostream &operator<<(std::ostream &out, const Rational &a) {
  out << a.toString();
  return out;
}

template<unsigned M, unsigned N, typename Field = Rational>
class Matrix {
public:
    Matrix() {
      for (unsigned i = 0; i < M; ++i) {
        std::vector<Field> a;
        for (unsigned j = 0; j < N; ++j) {
          if (i == j) {
            a.push_back(Field(1));
          } else {
            a.push_back(Field(0));
          }
        }
        nums_.push_back(a);
      }
    }


    Matrix(const Matrix<N, M, Field> &matrix_1) {
      for (unsigned i = 0; i < M; ++i) {
        std::vector<Field> String;
        for (unsigned j = 0; j < N; ++j) {
          String.push_back(matrix_1[i][j]);
        }
        nums_.push_back(String);
      }
    }

    template<typename T>
    Matrix(std::vector<std::vector<T>> x) {
      for (unsigned i = 0; i < M; ++i) {
        std::vector<Field> String;
        for (unsigned j = 0; j < N; ++j) {
          Field b(x[i][j]);
          String.push_back(b);
        }
        nums_.push_back(String);
      }
    }

    std::vector<Field> operator[](unsigned i) const {
      return nums_[i];
    }

    std::vector<Field> &operator[](unsigned i) {
      return nums_[i];
    }


    std::vector<Field> getRow(unsigned i) {
      return nums_[i];
    }

    std::vector<Field> getColumn(unsigned i) {
      std::vector<Field> ans;
      for (unsigned j = 0; j < M; ++j) {
        ans.push_back(nums_[i][j]);
      }
      return ans;
    }

    Matrix<N, M, Field> transposed() const {
      std::vector<std::vector<Field>> ans;
      for (unsigned i = 0; i < N; ++i) {
        std::vector<Field> String;
        for (unsigned j = 0; j < M; ++j) {
          String.push_back(nums_[j][i]);
        }
        ans.push_back(String);
      }
      return Matrix<N, M, Field>(ans);
    }

    Field trace() {
      static_assert(N == M);
      Field ans = 0;
      for (unsigned i = 0; i < M; ++i) {
        ans += nums_[i][i];
      }
      return ans;
    }

    Matrix<M, N, Field> &operator+=(const Matrix<M, N, Field> &matrix_1) {
      for (unsigned i = 0; i < M; ++i) {
        for (unsigned j = 0; j < N; ++j) {
          nums_[i][j] += matrix_1[i][j];
        }
      }
      return *this;
    }

    Matrix<M, N, Field> &operator-=(const Matrix<M, N, Field> &matrix_1) {
      for (unsigned i = 0; i < M; ++i) {
        for (unsigned j = 0; j < N; ++j) {
          nums_[i][j] -= matrix_1[i][j];
        }
      }
      return *this;
    }

    Matrix<N, N, Field> &operator*=(const Matrix<N, N, Field> &matrix_1) {
      static_assert(N == M);
      std::vector<std::vector<Field>> nums;
      for (unsigned i = 0; i < M; ++i) {
        std::vector<Field> String;
        for (unsigned j = 0; j < N; ++j) {
          String.push_back(nums_[i][j]);
          nums_[i][j] = 0;
        }
        nums.push_back(String);
      }
      for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < N; ++j) {
          for (unsigned t = 0; t < N; ++t) {
            nums_[i][t] += nums[i][j] * matrix_1[j][t];
          }
        }
      }
      return *this;
    }

    Field det() const {
      static_assert(N == M);
      std::vector<std::vector<Field>> nums;
      for (unsigned i = 0; i < M; ++i) {
        std::vector<Field> String;
        for (unsigned j = 0; j < N; ++j) {
          String.push_back(nums_[i][j]);
        }
        nums.push_back(String);
      }
      Field ans = 1;
      Field n = 0;
      int count = 0;
      for (unsigned i = 0; i < N; ++i) {
        unsigned j = i;
        while (j < N && nums[j][i] == n) {
          ++j;
        }
        if (j == N) {
          return n;
        }
        if (j != i) {
          ++count;
          std::swap(nums[i], nums[j]);
        }
        ans *= nums[i][i];
        Field d = nums[i][i];
        for (unsigned k = i; k < N; ++k) {
          nums[i][k] /= d;
        }
        for (unsigned r = i + 1; r < N; ++r) {
          Field w = nums[r][i];
          for (unsigned k = i; k < N; ++k) {
            nums[r][k] -= w * nums[i][k];
          }
        }
      }
      if (count % 2 == 1) {
        ans = -ans;
      }
      return ans;
    }

    int rank() const {
      std::vector<std::vector<Field>> nums;
      for (unsigned i = 0; i < M; ++i) {
        std::vector<Field> String;
        for (unsigned j = 0; j < N; ++j) {
          String.push_back(nums_[i][j]);
        }
        nums.push_back(String);
      }
      unsigned j = 0;
      unsigned i = 0;
      int counter = N;
      Field n = 0;
      while (i < M || j < N) {
        unsigned k = i;
        while (k < M && nums[k][j] == n) {
          ++k;
        }
        if (k == M) {
          ++j;
          --counter;
          break;
        }
        if (k != i) {
          std::swap(nums[i], nums[k]);
        }
        Field d = nums[i][j];
        for (unsigned t = j; t < N; ++t) {
          nums[i][t] /= d;
        }
        for (unsigned t = i + 1; t < M; ++t) {
          Field c = nums[t][j];
          for (unsigned q = j; q < N; ++q) {
            nums[t][q] -= c * nums[i][q];
          }
        }
        ++i;
        ++j;
      }
      return counter;
    }

    Matrix<N, N, Field> inverted() {
      static_assert(N == M);
      Field n = 0;
      std::vector<std::vector<Field>> nums;
      for (unsigned i = 0; i < M; ++i) {
        std::vector<Field> String;
        for (unsigned j = 0; j < N; ++j) {
          String.push_back(nums_[i][j]);
        }
        nums.push_back(String);
      }
      Matrix<N, N, Field> one;
      for (unsigned i = 0; i < N; ++i) {
        unsigned j = i;
        while (j < N && nums[j][i] == n) {
          ++j;
        }
        if (j != i) {
          std::swap(nums[i], nums[j]);
          std::swap(one[i], one[j]);
        }
        Field d = nums[i][i];
        for (unsigned k = 0; k < N; ++k) {
          nums[i][k] /= d;
          one[i][k] /= d;
        }
        for (unsigned r = i + 1; r < N; ++r) {
          Field w = nums[r][i];
          for (unsigned k = 0; k < N; ++k) {
            nums[r][k] -= w * nums[i][k];
            one[r][k] -= w * one[i][k];
          }
        }
      }
      unsigned j = N - 1;
      while (j >= 0) {
        for (unsigned i = 0; i < j; ++i) {
          Field t_2 = nums[i][j];
          for (unsigned k = 0; k < N; ++k) {
            one[i][k] -= one[j][k] * t_2;
            nums[i][k] -= nums[j][k] * t_2;
          }
        }
        if (j == 0) {
          return one;
        }
        --j;
      }
    }

    Matrix<M, N, Field> &invert() {
      *this = this->inverted();
      return *this;
    }

private:
    std::vector<std::vector<Field>> nums_;
};

template<unsigned N, typename Field = Rational>
using SquareMatrix = Matrix<N, N, Field>;

template<unsigned M, unsigned N, typename Field = Rational>
Matrix<M, N, Field> operator+(const Matrix<M, N, Field> &matrix_1, const Matrix<M, N, Field> &matrix_2) {
  Matrix<M, N, Field> m_3(matrix_1);
  m_3 += matrix_2;
  return m_3;
}

template<unsigned M, unsigned N, typename Field = Rational>
Matrix<M, N, Field> operator-(const Matrix<M, N, Field> &matrix_1, const Matrix<M, N, Field> &matrix_2) {
  Matrix<M, N, Field> m_3(matrix_1);
  m_3 -= matrix_2;
  return m_3;
}

template<unsigned M, unsigned N, unsigned K, typename Field = Rational>
Matrix<M, K, Field> operator*(const Matrix<M, N, Field> &matrix_1, const Matrix<N, K, Field> &matrix_2) {
  std::vector<std::vector<Field>> matrix(M);
  for (unsigned i = 0; i < M; ++i) {
    matrix[i].resize(K);
  }
  for (unsigned i = 0; i < M; ++i) {
    for (unsigned j = 0; j < N; ++j) {
      for (unsigned t = 0; t < K; ++t) {
        matrix[i][t] += matrix_1[i][j] * matrix_2[j][t];
      }
    }
  }
  return Matrix<M, K, Field>(matrix);
}


template<unsigned M, unsigned N, typename Field = Rational, typename T>
Matrix<M, N, Field> operator*(const Matrix<M, N, Field> &matrix_1, T x) {
  std::vector<std::vector<Field>> ans;
  for (unsigned i = 0; i < M; ++i) {
    std::vector<Field> a;
    for (unsigned j = 0; j < N; ++j) {
      a.push_back(matrix_1[i][j] * x);
    }
    ans.push_back(a);
  }
  return Matrix<M, N, Field>(ans);
}

template<unsigned M, unsigned N, typename Field = Rational, typename T>
Matrix<M, N, Field> operator*(T x, const Matrix<M, N, Field> &matrix_1) {
  std::vector<std::vector<Field>> ans;
  for (unsigned i = 0; i < M; ++i) {
    std::vector<Field> a;
    for (unsigned j = 0; j < N; ++j) {
      a.push_back(matrix_1[i][j] * x);
    }
    ans.push_back(a);
  }
  return Matrix<M, N, Field>(ans);
}


template<unsigned M_1, unsigned N_1, typename Field_1 = Rational, unsigned M_2, unsigned N_2, typename Field_2 = Rational>
bool operator==(const Matrix<M_1, N_1, Field_1> &matrix_1, const Matrix<M_2, N_2, Field_2> &matrix_2) {
  if (M_2 != M_1 || N_2 != N_1) {
    return false;
  }
  for (unsigned i = 0; i < M_1; ++i) {
    for (unsigned j = 0; j < N_1; ++j) {
      if (matrix_2[i][j] != matrix_1[i][j]) {
        return false;
      }
    }
  }
  return true;
}

template<unsigned M_1, unsigned N_1, typename Field_1 = Rational, unsigned M_2, unsigned N_2, typename Field_2 = Rational>
bool operator!=(const Matrix<M_1, N_1, Field_1> &matrix_1, const Matrix<M_2, N_2, Field_2> &matrix_2) {
  return !(matrix_1 == matrix_2);
}
