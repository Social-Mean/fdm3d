#include <fdm3d/fdm3dsolver.h>

FDM3DSolver::FDM3DSolver(int Nx, int Ny, int Nz)
  : Nx{ Nx }, Ny{ Ny }, Nz{ Nz }, npts{ Nx * Ny * Nz } {
  initBimap();
  for (int tag = 0; tag < Nx * Ny * Nz; ++tag) {
	const auto& [i, j, k] = tag2idx(tag);
	std::cout << std::format("({}, {}, {}) : {}", i, j, k, tag) << std::endl;
  }
  std::cout << std::endl;
}

const int& FDM3DSolver::idx2tag(const std::tuple<int, int, int>& idx) const {
  auto tag = idx_tag_bimap.left.at(idx);
  return tag;
}

const int& FDM3DSolver::idx2tag(const int& i, const int& j, const int& k) const {
  auto idx = std::make_tuple(i, j, k);
  auto tag = idx2tag(idx);
  std::cout << std::format("({}, {}, {}) : {}", i, j, k, tag) << std::endl;
  return tag;
}

const std::tuple<int, int, int>& FDM3DSolver::tag2idx(const int& tag) const {
  const auto& idx = idx_tag_bimap.right.at(tag);
  return idx;
}

void FDM3DSolver::initBimap() {
  for (int i = 0; i < Nx; ++i) {
	for (int j = 0; j < Ny; ++j) {
	  for (int k = 0; k < Nz; ++k) {
		auto idx = std::make_tuple(i, j, k);
		int tag = i * (Ny * Nz) + j * Nz + k;
		idx_tag_bimap.left.insert({ idx, tag });
	  }
	}
  }
 // for (int tag = 0; tag < npts; ++tag) {
	//auto [tmp, k] = std::div(tag, Nz);
	//auto [i, j] = std::div(tag, Ny * Nz);
	//auto idx = std::make_tuple(i, j, k);
	//idx_tag_bimap.right.insert({ tag, idx });
 // }
}

void FDM3DSolver::solve() {
  setMatrix();
  solveMatrix();
}

void FDM3DSolver::setMatrix() {
  initMatrix();
  setEquation();
  setBCs();
}

void FDM3DSolver::initMatrix() {
  A.resize(npts, npts);
  b.resize(npts);
  x.resize(npts);
}

void FDM3DSolver::setEquation() {
  using T = Eigen::Triplet<double>;
  auto listSize = 7 * npts;
  std::vector<T> tripletList;
  tripletList.reserve(listSize);
  for (auto i = 1; i < Nx-1; ++i) {
	for (auto j = 1; j < Ny-1; ++j) {
	  for (auto k = 1; k < Nz-1; ++k) {
		auto tag = idx2tag(i, j, k);
		tripletList.emplace_back(tag, tag, -6);
		tripletList.emplace_back(tag, idx2tag(i + 1, j, k), 1);
		tripletList.emplace_back(tag, idx2tag(i - 1, j, k), 1);
		tripletList.emplace_back(tag, idx2tag(i, j + 1, k), 1);
		tripletList.emplace_back(tag, idx2tag(i, j - 1, k), 1);
		tripletList.emplace_back(tag, idx2tag(i, j, k + 1), 1);
		tripletList.emplace_back(tag, idx2tag(i, j, k - 1), 1);
		b(tag) = 0;
	  }
	}
  }
  for (const auto& triplet : tripletList) {
	std::cout << std::format("({}, {}) : {}", triplet.row(), triplet.col(), triplet.value()) << std::endl;
  }
  A.setFromTriplets(tripletList.begin(), tripletList.end());
}

void FDM3DSolver::setBCs() {
  setDBCs();
}

void FDM3DSolver::setDBCs() {
  // TODO
  auto tag = idx2tag(std::make_tuple(Nx / 2, Ny / 2, Nz / 2));
  setDBC(tag, 0);

  for (int i = 0; i < Nx; ++i) {
	for (int j = 0; j < Ny; ++j) {
	  for (int k = 0; k < Nz; ++k) {
		if (i == 0 || i == Nx - 1 || j == 0 || j == Ny - 1 || k == 0 || k == Nz - 1) {
		  auto tag = idx2tag(std::make_tuple(i, j, k));
		  setDBC(tag, 1);
		}
	  }
	}
  }
}

void FDM3DSolver::setDBC(int tag, double value) {
  A.prune([tag](int row, int col, double value) { return row != tag && col != tag; });
  A.insert(tag, tag) = 1;
  b(tag) = value;
}

void FDM3DSolver::solveMatrix() {
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;
  cg.setTolerance(1e-3);
  cg.setMaxIterations(50);
  std::cout << A << std::endl;
  std::cout << b << std::endl;
  x = cg.compute(A).solve(b);
}

//const Eigen::Tensor<double, 3>& FDM3DSolver::getPhi() const {
  //Eigen::TensorMap<Eigen::Tensor<double, 3>> phi(x.data(), Nx, Ny, Nz);
  //return phi;
//}

Eigen::VectorXd FDM3DSolver::getX() const {
  return x;
}