#pragma once
#include <iostream>
#include <map>
#include <boost/bimap.hpp>
#include <tuple>
#include <string>
#include <format>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <array>
#include <unsupported/Eigen/CXX11/Tensor>

class FDM3DSolver {
public:
  FDM3DSolver(int Nx, int Ny, int Nz);
  void solve();
  //const Eigen::Tensor<double, 3>& getPhi() const;
  Eigen::VectorXd getX() const;
private:
  const int Nx;
  const int Ny;
  const int Nz;
  const int npts;
  boost::bimap<std::tuple<int, int, int>, int> idx_tag_bimap;
  Eigen::SparseMatrix<double> A;
  Eigen::VectorXd b;
  Eigen::VectorXd x;
  
	
  const int& idx2tag(const std::tuple<int, int, int>& idx) const;
  const int& idx2tag(const int& i, const int& j, const int& k) const;
  const std::tuple<int, int, int>& tag2idx(const int& tag) const;
  void initBimap();
  void setMatrix();
  void initMatrix();
  void setEquation();
  void setBCs();
  void setDBCs();
  void solveMatrix();
  void setDBC(int tag, double value);
};