#include <vector>
#include <string>

using std::vector;
using std::string;


template<typename T>
class LightMatrix
{
public:
  LightMatrix() : nrows_(0), ncols_(0), data_(nullptr), is_data_owner_(true) {}
  LightMatrix(int nrows, int ncols) : nrows_(nrows), ncols_(ncols), is_data_owner_(true) { data_ = new T[nrows*ncols]; }
  LightMatrix(int nrows, int ncols, T* data, bool is_data_owner = false);
  LightMatrix(const LightMatrix<T>& other);
  ~LightMatrix() { if (is_data_owner_) delete[] data_; }

  int size() const { return ncols_*nrows_; }
  int cols() const { return ncols_; }
  void resize(int nrows, int ncols);
  void fill(const T& val);
  void set(const T* const data);
  void set_row(int i, const T& val);
  void set_col(int i, const T& val);
  void set_col(int i, const T* const data);
  void set_identity();
  void set_block(int row_off, int col_off, const LightMatrix<T>& block);
  void scale_col(int i, const T& val);
  void scale_row(int i, const T& val);
  const T* get_col(int i) const { return &((*this)(0, i)); };
  void transpose_in_place();
  T& operator()(int i_row, int i_col) { return data_[i_col*nrows_ + i_row]; }
  const T& operator()(int i_row, int i_col) const { return data_[i_col*nrows_ + i_row]; }
  LightMatrix<T>& operator=(const LightMatrix<T>& other);

  bool is_data_owner_;
  int nrows_;
  int ncols_;
  T *data_;
};

template<typename T, typename U, typename V>
void mat_mult(const LightMatrix<T>& lhs, const LightMatrix<U>& rhs, LightMatrix<V> *pout);

////////////////////////////////////////////////////////////
//////////////////// Definitions ///////////////////////////
////////////////////////////////////////////////////////////

template<typename T>
LightMatrix<T>& LightMatrix<T>::operator=(const LightMatrix<T>& other)
{ 
  resize(other.nrows_, other.ncols_);
  for (int i = 0; i < size(); i++)
    data_[i] = other.data_[i];
  is_data_owner_ = true;
  return *this; 
}

template<typename T>
LightMatrix<T>::LightMatrix(const LightMatrix<T>& other)
  : nrows_(0), ncols_(0), data_(nullptr), is_data_owner_(true)
{
  resize(other.nrows_, other.ncols_);
  for (int i = 0; i < size(); i++)
    data_[i] = other.data_[i];
  is_data_owner_ = true;
}

template<typename T>
LightMatrix<T>::LightMatrix(int nrows, int ncols, T* data, bool is_data_owner)
  : nrows_(nrows), ncols_(ncols), data_(data), is_data_owner_(is_data_owner) {}

template<typename T>
void LightMatrix<T>::set_block(int row_off, int col_off, const LightMatrix<T>& block)
{
  for (int i_col = 0; i_col < block.ncols_; i_col++)
    for (int i_row = 0; i_row < block.nrows_; i_row++)
      (*this)(i_row + row_off, i_col + col_off) = block(i_row, i_col);
}

template<typename T>
void LightMatrix<T>::transpose_in_place()
{
  for (int j = 0; j < ncols_; j++)
  {
    for (int i = j+1; i < nrows_; i++)
    {
      T tmp = (*this)(i, j);
      (*this)(i, j) = (*this)(j, i);
      (*this)(j, i) = tmp;
    }
  }
}

template<typename T>
void LightMatrix<T>::set_identity()
{
  for (int i_col = 0; i_col < ncols_; i_col++)
    for (int i_row = 0; i_row < nrows_; i_row++)
      if (i_col == i_row)
        (*this)(i_col, i_row) = T(1);
      else
        (*this)(i_col, i_row) = T(0);
}

template<typename T>
void LightMatrix<T>::scale_col(int i, const T& val)
{
  for (int j = 0; j < nrows_; j++)
    (*this)(j, i) *= val;
}

template<typename T>
void LightMatrix<T>::scale_row(int i, const T& val)
{
  for (int j = 0; j < ncols_; j++)
    (*this)(i, j) *= val;
}

template<typename T>
void LightMatrix<T>::set(const T* const data)
{
  for (int j = 0; j < size(); j++)
    data_[j] = data[j];
}

template<typename T>
void LightMatrix<T>::set_row(int i, const T& val)
{
  for (int j = 0; j < ncols_; j++)
    (*this)(i, j) = val;
}

template<typename T>
void LightMatrix<T>::set_col(int i, const T& val)
{
  for (int j = 0; j < nrows_; j++)
    (*this)(j, i) = val;
}

template<typename T>
void LightMatrix<T>::set_col(int i, const T* const data)
{
  for (int j = 0; j < nrows_; j++)
    (*this)(j, i) = data[j];
}

template<typename T>
void LightMatrix<T>::fill(const T& val)
{
  for (int i = 0; i < size(); i++)
    data_[i] = val;
}

template<typename T>
void LightMatrix<T>::resize(int nrows, int ncols)
{
  if (nrows_*ncols_ != nrows*ncols)
  {
    if (is_data_owner_)
      delete[] data_;
    if (nrows*ncols > 0)
      data_ = new T[ncols*nrows];
    else
      data_ = nullptr;
  }
  ncols_ = ncols;
  nrows_ = nrows;
}

template<typename T, typename U, typename V>
void mat_mult(const LightMatrix<T>& lhs, const LightMatrix<U>& rhs, LightMatrix<V> *pout)
{
  auto& out = *pout;
  out.resize(lhs.nrows_, rhs.ncols_);
  for (int i = 0; i < lhs.nrows_; i++)
  {
    for (int k = 0; k < rhs.ncols_; k++)
    {
      out(i, k) = lhs(i, 0) * rhs(0, k);
      for (int j = 1; j < lhs.ncols_; j++)
      {
        out(i, k) += lhs(i, j) * rhs(j, k);
      }
    }
  }
}

typedef struct { int verts[3]; } Triangle;
class HandModelLightMatrix
{
public:
  vector<string> bone_names;
  vector<int> parents; // assumimng that parent is earlier in the order of bones
  vector<LightMatrix<double>> base_relatives;
  vector<LightMatrix<double>> inverse_base_absolutes;
  LightMatrix<double> base_positions;
  LightMatrix<double> weights;
  vector<Triangle> triangles;
  bool is_mirrored;
};

class HandDataLightMatrix
{
public:
  HandModelLightMatrix model;
  vector<int> correspondences;
  LightMatrix<double> points;
};

